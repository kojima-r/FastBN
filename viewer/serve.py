#!/usr/bin/env python3
# =============================================================================
# viewer/serve.py
#   cosmos.gl ビューア (viewer/index.html) 用のローカルサーバ。
#   Python 標準ライブラリのみ。npm / bundler / 外部通信は不要。
#
#   * 指定したディレクトリ以下から**推定済みネットワーク**を自動的に探す
#       学習網        : edges.tsv        + edges_named.tsv
#       コンセンサス網 : integ_edges2.tsv + integ_edges_named.tsv
#     付随ファイル (あれば自動で使う):
#       edge_importance.tsv / integ_edge_importance.tsv : エッジ重要度
#       integ_edges_score.tsv                           : ブートストラップ確率
#       eval_hc_edges.tsv / eval_bs_edges.tsv           : 正解構造との判定 (TP/FP/...)
#       ../data/var_map.tsv                             : 変数対応表 (変数総数の表示)
#       ../target_genes.txt                             : 注目遺伝子
#   * ブラウザには JSON API で渡す (/api/networks, /api/network?id=...)
#
# 使い方:
#   python3 viewer/serve.py                    # 既定: リポジトリ全体を探索して起動
#   python3 viewer/serve.py --root my_analysis # 特定の解析ディレクトリだけ
#   python3 viewer/serve.py --port 9000 --no-browser
#   python3 viewer/serve.py --list             # 見つかったネットワークを表示して終了
#   python3 viewer/serve.py --fetch-vendor     # cosmos.gl を取り直す (要ネットワーク)
#
# ノード番号は入力 TSV の列位置なので、edges.tsv と edges_named.tsv の**行対応**から
# 名前を復元する (var_map.tsv を作り直しても図がずれない、という repo の約束と同じ)。
# =============================================================================
import argparse
import html
import json
import os
import posixpath
import re
import socket
import sys
import threading
import urllib.parse
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

VIEWER_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(VIEWER_DIR)

# --- 探索対象のネットワーク定義 ---------------------------------------------
NET_SPECS = [
    {
        "kind": "learned",
        "label": "学習網 (Hill-Climb)",
        "edges": "edges.tsv",
        "named": "edges_named.tsv",
        "importance": "edge_importance.tsv",
        "prob": None,
        "eval": "eval_hc_edges.tsv",
    },
    {
        "kind": "consensus",
        "label": "コンセンサス網 (bootstrap)",
        "edges": "integ_edges2.tsv",
        "named": "integ_edges_named.tsv",
        "importance": "integ_edge_importance.tsv",
        "prob": "integ_edges_score.tsv",
        "eval": "eval_bs_edges.tsv",
    },
]

# edge_importance.tsv の列 (ヘッダは 6 列しか書かれていないので位置で読む)
IMP_COLUMNS = ["u", "v", "dlogL", "dBIC", "dK2", "dBDeu",
               "mean_dlogL_per_sample", "std_dlogL_per_sample"]
METRICS = IMP_COLUMNS[2:]

# 走査から外すディレクトリ (巨大 or 無関係)
PRUNE_DIRS = {
    ".git", ".hg", ".svn", "__pycache__", "node_modules", ".venv", "venv",
    ".claude", ".codex", ".agents", ".claude-plugin", ".codex-plugin",
    "viewer", "figures", "figures_bs", "subsets", "bs", "groups",
    "sample00", "sample00_bs", "logs", "source", "networks",
}
PRUNE_PREFIXES = ("data",)   # data/, data_bbknn_.../ など (var_map は直接参照する)

# 既定で最初に表示するネットワークの優先順 (前方一致)
DEFAULT_PREFERENCE = [
    "example_bulk/out::learned",
    "example_bulk/out::consensus",
    "example_sc/run_bbknn_bin100/out::learned",
    "example_sc/run_bbknn_binall/out::learned",
    "example_sachs/",
    "example_bnlearn/",
]


def log(*a):
    print("[viewer]", *a, file=sys.stderr, flush=True)


# --- ファイル読み込み --------------------------------------------------------

def count_data_lines(path):
    """空行・コメント行を除いた行数。"""
    n = 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fp:
            for line in fp:
                if line.strip() and not line.startswith("#"):
                    n += 1
    except OSError:
        return 0
    return n


def read_pairs(path):
    """`a<TAB>b` の 2 列を読む (ヘッダ行は自動判定して飛ばす)。"""
    out = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fp:
            for i, line in enumerate(fp):
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue
                a = line.split("\t") if "\t" in line else line.split()
                if len(a) < 2:
                    continue
                if i == 0 and a[0].lower() in ("u", "from", "source", "parent"):
                    continue
                out.append((a[0].strip(), a[1].strip()))
    except OSError:
        return []
    return out


def read_importance(path):
    """(u, v) -> {metric: |value|} を返す。値は visualize.py と同じく絶対値。"""
    imp = {}
    if not path or not os.path.isfile(path):
        return imp
    with open(path, "r", encoding="utf-8", errors="replace") as fp:
        for i, line in enumerate(fp):
            a = line.rstrip("\n").split("\t")
            if len(a) < 3:
                continue
            if i == 0 and not a[0].lstrip("-").isdigit():
                continue
            try:
                u, v = int(a[0]), int(a[1])
            except ValueError:
                continue
            vals = {}
            for j, m in enumerate(METRICS, start=2):
                if j < len(a):
                    try:
                        vals[m] = abs(float(a[j]))
                    except ValueError:
                        pass
            imp[(u, v)] = vals
    return imp


def read_prob(path):
    """integ_edges_score.tsv (u v count prob) -> (u, v) -> prob"""
    prob = {}
    if not path or not os.path.isfile(path):
        return prob
    with open(path, "r", encoding="utf-8", errors="replace") as fp:
        for line in fp:
            a = line.rstrip("\n").split("\t")
            if len(a) < 4:
                continue
            try:
                prob[(int(a[0]), int(a[1]))] = float(a[3])
            except ValueError:
                continue
    return prob


def read_eval(path):
    """eval_*_edges.tsv (u v status; 遺伝子名) -> (name_u, name_v) -> status"""
    ev = {}
    if not path or not os.path.isfile(path):
        return ev
    with open(path, "r", encoding="utf-8", errors="replace") as fp:
        for i, line in enumerate(fp):
            a = line.rstrip("\n").split("\t")
            if len(a) < 3:
                continue
            if i == 0 and a[2].strip().lower() == "status":
                continue
            ev[(a[0].strip(), a[1].strip())] = a[2].strip()
    return ev


def read_targets(paths):
    """注目遺伝子リスト (1 行 1 名前、# はコメント)。"""
    targets = set()
    for p in paths:
        if not p or not os.path.isfile(p):
            continue
        with open(p, "r", encoding="utf-8", errors="replace") as fp:
            for line in fp:
                line = line.strip()
                if line and not line.startswith("#"):
                    targets.add(line.split("\t")[0])
    return targets


def read_var_map(path):
    """var_map.tsv -> (総変数数, index -> column_name)"""
    if not path or not os.path.isfile(path):
        return 0, {}
    names = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fp:
        for i, line in enumerate(fp):
            a = line.rstrip("\n").split("\t")
            if i == 0 or len(a) < 2:
                continue
            try:
                names[int(a[0])] = a[1]
            except ValueError:
                continue
    return len(names), names


# --- ネットワークの探索 ------------------------------------------------------

def discover(roots, max_networks):
    """roots 以下からネットワークを探す。戻り値は catalog エントリのリスト。"""
    found = []
    seen_dirs = set()
    for root in roots:
        root = os.path.abspath(root)
        if not os.path.isdir(root):
            log(f"警告: ディレクトリがありません: {root}")
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = sorted(
                d for d in dirnames
                if d not in PRUNE_DIRS and not d.startswith(".")
                and not d.startswith(PRUNE_PREFIXES)
            )
            real = os.path.realpath(dirpath)
            if real in seen_dirs:
                dirnames[:] = []
                continue
            seen_dirs.add(real)
            fileset = set(filenames)
            for spec in NET_SPECS:
                if spec["edges"] not in fileset:
                    continue
                entry = build_entry(dirpath, spec, fileset)
                if entry is not None:
                    found.append(entry)
            if len(found) >= max_networks:
                log(f"警告: {max_networks} 件でネットワークの探索を打ち切りました "
                    "(--max-networks で変更)")
                return found
    return found


def build_entry(dirpath, spec, fileset):
    edges_path = os.path.join(dirpath, spec["edges"])
    n_edges = count_data_lines(edges_path)
    if n_edges == 0:
        return None
    named_path = os.path.join(dirpath, spec["named"])
    has_named = spec["named"] in fileset

    run_dir = os.path.dirname(os.path.abspath(dirpath))          # out/ の親
    rel = os.path.relpath(dirpath, REPO_DIR).replace(os.sep, "/")
    group = rel.split("/")[0] if "/" in rel else rel

    def opt(name):
        return os.path.join(dirpath, name) if name and name in fileset else None

    return {
        "id": f"{rel}::{spec['kind']}",
        "dir": rel,
        "kind": spec["kind"],
        "kind_label": spec["label"],
        "group": group,
        "label": f"{rel} — {spec['label']}",
        "n_edges": n_edges,
        "paths": {
            "edges": edges_path,
            "named": named_path if has_named else None,
            "importance": opt(spec["importance"]),
            "prob": opt(spec["prob"]),
            "eval": opt(spec["eval"]),
            "var_map": first_existing([
                os.path.join(run_dir, "data", "var_map.tsv"),
                os.path.join(dirpath, "var_map.tsv"),
            ]),
            "targets": [p for p in (
                os.path.join(run_dir, "target_genes.txt"),
                os.path.join(os.path.dirname(run_dir), "target_genes.txt"),
            ) if os.path.isfile(p)],
        },
    }


def first_existing(paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return None


def catalog_view(entries):
    """API に返す軽量版 (パスは外に出さない)。"""
    out = []
    for e in entries:
        p = e["paths"]
        out.append({
            "id": e["id"], "dir": e["dir"], "kind": e["kind"],
            "kind_label": e["kind_label"], "group": e["group"],
            "label": e["label"], "n_edges": e["n_edges"],
            "has_importance": bool(p["importance"]),
            "has_prob": bool(p["prob"]),
            "has_eval": bool(p["eval"]),
            "has_names": bool(p["named"]),
        })
    return out


def pick_default(entries):
    ids = [e["id"] for e in entries]
    for pref in DEFAULT_PREFERENCE:
        for i in ids:
            if i.startswith(pref):
                return i
    # 見つからなければ「重要度つき」で最小のものを選ぶ (最初の確認が軽く済む)
    scored = sorted(
        entries,
        key=lambda e: (not e["paths"]["importance"], e["n_edges"]),
    )
    return scored[0]["id"] if scored else None


# --- ネットワークの読み込み (JSON ペイロード) --------------------------------

def load_network(entry):
    p = entry["paths"]
    warnings = []
    idx_edges = read_pairs(p["edges"])
    named_edges = read_pairs(p["named"]) if p["named"] else []

    if named_edges and len(named_edges) != len(idx_edges):
        warnings.append(
            f"{os.path.basename(p['named'])} の行数 ({len(named_edges)}) が "
            f"{os.path.basename(p['edges'])} ({len(idx_edges)}) と一致しないため、"
            "名前は使わずノード番号で表示します")
        named_edges = []
    elif not named_edges:
        warnings.append("名前ファイルが無いのでノード番号で表示します")

    # ノード表 (元の列インデックス -> 連番)
    order = []
    compact = {}
    names = {}
    for k, (a, b) in enumerate(idx_edges):
        try:
            u, v = int(a), int(b)
        except ValueError:
            continue
        for node in (u, v):
            if node not in compact:
                compact[node] = len(order)
                order.append(node)
        if named_edges:
            nu, nv = named_edges[k]
            names.setdefault(u, nu)
            names.setdefault(v, nv)

    n_vars, var_names = read_var_map(p["var_map"])
    if var_names and names:
        sample = [i for i in order[:200] if i in var_names and i in names]
        mismatch = sum(1 for i in sample if var_names[i] != names[i])
        if sample and mismatch > len(sample) * 0.1:
            warnings.append(
                f"var_map.tsv の名前が edges_named.tsv と {mismatch}/{len(sample)} 件"
                "食い違っています (前処理をやり直した var_map かもしれません)。"
                "表示には edges_named.tsv の行対応を使います")

    imp = read_importance(p["importance"])
    prob = read_prob(p["prob"])
    ev = read_eval(p["eval"])
    targets = read_targets(p["targets"])

    nodes_in = [0] * len(order)
    nodes_out = [0] * len(order)
    links = []
    metric_values = {m: [] for m in METRICS}
    prob_values = []
    eval_values = []

    for k, (a, b) in enumerate(idx_edges):
        try:
            u, v = int(a), int(b)
        except ValueError:
            continue
        su, tv = compact[u], compact[v]
        links.append(su)
        links.append(tv)
        nodes_out[su] += 1
        nodes_in[tv] += 1
        vals = imp.get((u, v), {})
        for m in METRICS:
            metric_values[m].append(vals.get(m))
        if prob:
            prob_values.append(prob.get((u, v)))
        if ev:
            nu = names.get(u, str(u))
            nv = names.get(v, str(v))
            eval_values.append(ev.get((nu, nv)))

    present_metrics = [m for m in METRICS
                       if any(x is not None for x in metric_values[m])]

    nodes = []
    for compact_idx, orig in enumerate(order):
        name = names.get(orig) or var_names.get(orig) or f"#{orig}"
        nodes.append({
            "n": name,
            "i": orig,
            "in": nodes_in[compact_idx],
            "out": nodes_out[compact_idx],
            "t": 1 if name in targets else 0,
        })

    payload = {
        "id": entry["id"],
        "dir": entry["dir"],
        "kind": entry["kind"],
        "kind_label": entry["kind_label"],
        "label": entry["label"],
        "n_nodes": len(nodes),
        "n_edges": len(links) // 2,
        "n_vars": n_vars,
        "n_targets": sum(n["t"] for n in nodes),
        "metrics": present_metrics,
        "has_prob": bool(prob_values),
        "has_eval": bool(eval_values),
        "warnings": warnings,
        "nodes": nodes,
        "links": links,
        "link_metrics": {m: metric_values[m] for m in present_metrics},
        "link_prob": prob_values or None,
        "link_eval": eval_values or None,
        "files": {k: (os.path.relpath(v, REPO_DIR) if isinstance(v, str) else
                      [os.path.relpath(x, REPO_DIR) for x in v] if v else None)
                  for k, v in p.items()},
    }
    return payload


# --- HTTP -------------------------------------------------------------------

STATIC_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".map": "application/json; charset=utf-8",
    ".txt": "text/plain; charset=utf-8",
    ".md": "text/plain; charset=utf-8",
}


class State:
    def __init__(self, entries, default_id, quiet):
        self.entries = entries
        self.by_id = {e["id"]: e for e in entries}
        self.default_id = default_id
        self.quiet = quiet
        self.cache = {}
        self.lock = threading.Lock()

    def payload(self, net_id):
        with self.lock:
            if net_id in self.cache:
                return self.cache[net_id]
        entry = self.by_id.get(net_id)
        if entry is None:
            return None
        data = json.dumps(load_network(entry), ensure_ascii=False,
                          separators=(",", ":")).encode("utf-8")
        with self.lock:
            self.cache[net_id] = data
        return data


def make_handler(state):
    class Handler(BaseHTTPRequestHandler):
        server_version = "FastBNViewer/1.0"

        def log_message(self, fmt, *args):
            if not state.quiet:
                log(fmt % args)

        def _send(self, code, body, ctype, extra=None):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            for k, v in (extra or {}).items():
                self.send_header(k, v)
            self.end_headers()
            if self.command != "HEAD":
                self.wfile.write(body)

        def _json(self, obj, code=200):
            body = json.dumps(obj, ensure_ascii=False).encode("utf-8")
            self._send(code, body, "application/json; charset=utf-8")

        def _error(self, code, msg):
            body = (f"<!doctype html><meta charset=utf-8><title>{code}</title>"
                    f"<pre>{html.escape(msg)}</pre>").encode("utf-8")
            self._send(code, body, "text/html; charset=utf-8")

        def do_HEAD(self):
            self.do_GET()

        def do_GET(self):
            parsed = urllib.parse.urlsplit(self.path)
            path = urllib.parse.unquote(parsed.path)
            query = urllib.parse.parse_qs(parsed.query)

            if path == "/api/networks":
                self._json({
                    "networks": catalog_view(state.entries),
                    "default": state.default_id,
                    "repo": REPO_DIR,
                })
                return

            if path == "/api/network":
                net_id = (query.get("id") or [""])[0] or state.default_id
                body = state.payload(net_id)
                if body is None:
                    self._json({"error": f"unknown network id: {net_id}"}, 404)
                    return
                self._send(200, body, "application/json; charset=utf-8")
                return

            # 静的ファイル (viewer/ 以下だけ)
            rel = path.lstrip("/") or "index.html"
            rel = posixpath.normpath(rel)
            if rel.startswith("..") or os.path.isabs(rel):
                self._error(403, "forbidden")
                return
            full = os.path.join(VIEWER_DIR, rel.replace("/", os.sep))
            if os.path.isdir(full):
                full = os.path.join(full, "index.html")
            if not os.path.abspath(full).startswith(VIEWER_DIR + os.sep) or \
                    not os.path.isfile(full):
                self._error(404, f"not found: {path}")
                return
            ctype = STATIC_TYPES.get(os.path.splitext(full)[1].lower(),
                                     "application/octet-stream")
            with open(full, "rb") as fp:
                self._send(200, fp.read(), ctype)

    return Handler


def find_port(host, port, tries=20):
    for p in range(port, port + tries):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                s.bind((host, p))
                return p
            except OSError:
                continue
    raise SystemExit(f"[viewer] エラー: {port}〜{port + tries - 1} に空きポートがありません")


# --- cosmos.gl の取得 (任意) -------------------------------------------------

def fetch_vendor(version):
    """npm レジストリから UMD ビルドを取り直して viewer/vendor に置く。"""
    import hashlib
    import io
    import tarfile
    import urllib.request

    url = (f"https://registry.npmjs.org/@cosmos.gl/graph/-/graph-{version}.tgz")
    log(f"取得中: {url}")
    with urllib.request.urlopen(url, timeout=120) as resp:
        blob = resp.read()
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        js = tf.extractfile("package/dist/index.min.js").read()
        try:
            lic = tf.extractfile("package/LICENCE").read()
        except KeyError:
            lic = tf.extractfile("package/LICENSE").read()
    out_js = os.path.join(VIEWER_DIR, "vendor", f"cosmos.gl-{version}.min.js")
    os.makedirs(os.path.dirname(out_js), exist_ok=True)
    with open(out_js, "wb") as fp:
        fp.write(js)
    with open(os.path.join(VIEWER_DIR, "vendor", "LICENSE.cosmos.gl.txt"), "wb") as fp:
        fp.write(lic)
    log(f"書き出し: {os.path.relpath(out_js, REPO_DIR)} "
        f"({len(js) / 1024:.0f} KB, sha256={hashlib.sha256(js).hexdigest()[:16]}…)")

    # index.html の参照を張り替える
    index = os.path.join(VIEWER_DIR, "index.html")
    with open(index, "r", encoding="utf-8") as fp:
        src = fp.read()
    new = re.sub(r'vendor/cosmos\.gl-[0-9][^"\']*\.min\.js',
                 f"vendor/cosmos.gl-{version}.min.js", src)
    if new != src:
        with open(index, "w", encoding="utf-8") as fp:
            fp.write(new)
        log(f"index.html の script 参照を cosmos.gl-{version}.min.js に更新しました")


def check_vendor():
    vendor = os.path.join(VIEWER_DIR, "vendor")
    files = [f for f in os.listdir(vendor)] if os.path.isdir(vendor) else []
    if not any(f.startswith("cosmos.gl-") and f.endswith(".min.js") for f in files):
        log("警告: viewer/vendor に cosmos.gl のビルドがありません。"
            "`python3 viewer/serve.py --fetch-vendor` で取得してください。")


# --- main -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="推定済みネットワークを cosmos.gl で見るためのローカルサーバ")
    ap.add_argument("--root", action="append", default=None,
                    help="探索するディレクトリ (繰り返し指定可; 既定はリポジトリ全体)")
    ap.add_argument("--port", type=int, default=8765, help="ポート (既定 8765)")
    ap.add_argument("--host", default="127.0.0.1", help="バインドするホスト")
    ap.add_argument("--no-browser", action="store_true", help="ブラウザを開かない")
    ap.add_argument("--select", default=None, help="最初に表示するネットワーク ID")
    ap.add_argument("--max-networks", type=int, default=800,
                    help="探索を打ち切る件数 (既定 800)")
    ap.add_argument("--list", action="store_true",
                    help="見つかったネットワークを表示して終了")
    ap.add_argument("--fetch-vendor", action="store_true",
                    help="cosmos.gl の UMD ビルドを npm から取り直す")
    ap.add_argument("--vendor-version", default="3.4.1",
                    help="--fetch-vendor で取得するバージョン (既定 3.4.1)")
    ap.add_argument("--quiet", action="store_true", help="アクセスログを出さない")
    args = ap.parse_args()

    if args.fetch_vendor:
        fetch_vendor(args.vendor_version)
        return

    roots = args.root or [REPO_DIR]
    log(f"探索: {', '.join(os.path.relpath(r, REPO_DIR) if r != REPO_DIR else '.' for r in roots)}")
    entries = discover(roots, args.max_networks)
    if not entries:
        raise SystemExit(
            "[viewer] エラー: ネットワークが見つかりません。\n"
            "  edges.tsv (+ edges_named.tsv) か integ_edges2.tsv がある\n"
            "  ディレクトリを --root で指定してください。\n"
            "  まだ何も推定していない場合は example_bulk/run_all.sh を実行すると\n"
            "  1〜2 分で確認用のネットワークができます。")

    entries.sort(key=lambda e: (e["group"], e["dir"], e["kind"]))
    default_id = args.select or pick_default(entries)
    if default_id not in {e["id"] for e in entries}:
        log(f"警告: 指定された ID が見つかりません: {default_id}")
        default_id = pick_default(entries)

    groups = {}
    for e in entries:
        groups[e["group"]] = groups.get(e["group"], 0) + 1
    log(f"ネットワーク {len(entries)} 件: " +
        ", ".join(f"{g} {n}" for g, n in sorted(groups.items())))
    log(f"既定の表示: {default_id}")

    if args.list:
        for e in entries:
            p = e["paths"]
            flags = "".join([
                "I" if p["importance"] else "-",
                "P" if p["prob"] else "-",
                "E" if p["eval"] else "-",
                "N" if p["named"] else "-",
            ])
            print(f"{e['id']}\t{e['n_edges']} edges\t{flags}")
        print("\n凡例: I=重要度 P=ブートストラップ確率 E=正解との判定 N=遺伝子名",
              file=sys.stderr)
        return

    check_vendor()
    port = find_port(args.host, args.port)
    state = State(entries, default_id, args.quiet)
    httpd = ThreadingHTTPServer((args.host, port), make_handler(state))
    url = f"http://{args.host}:{port}/"
    log(f"起動しました: {url}  (Ctrl-C で停止)")
    if not args.no_browser:
        threading.Timer(0.4, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        log("停止しました")
    finally:
        httpd.server_close()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_report.py
==============
ベイジアンネットワーク解析の成果物 (図・重要度テーブル・データ要約) を
1 つの HTML レポートにまとめる汎用スクリプト。既定では図を相対リンクで
参照する (軽量)。--embed を付けると画像を base64 埋め込みし、出力 HTML
1 ファイルだけを共有・移動しても図が表示される自己完結型になる。
レポート上部のボタンでメトリクス (dlogL/dBIC/…) を切り替えて表示できる。

既定のパスは**カレントディレクトリ基準**。読み込むもの (存在すれば取り込み、
無ければそのセクションをスキップ):
  data/var_map.tsv                : 変数数などの要約
  data/expr_disc.tsv              : サンプル数・遺伝子数
  data/samples.tsv                : サンプルと群の一覧 (groups が無い場合の代替)
  groups/groups_manifest.tsv      : 実験群の一覧
  target_genes.txt                : 注目遺伝子 (ホワイトリスト)
  out/edges.tsv                   : 学習網のエッジ数
  out/integ_edges2.tsv            : コンセンサス網のエッジ数
  figures/                        : 学習網の図 (01〜05, subsets/)
  figures_bs/                     : コンセンサス網の図 (01〜06, subsets/)
  figures*/edge_importance_named_<metric>.tsv : 重要度降順エッジ表 (上位を表に)

出力:
  report.html (--out で変更可)
"""

import argparse
import base64
import html
import mimetypes
import os
import sys


def log(msg):
    print(f"[report] {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# 読み込みユーティリティ
# ---------------------------------------------------------------------------
def read_tsv(path, limit=None):
    """TSV を [header, row, row, ...] のリストで返す。無ければ None。"""
    if not path or not os.path.exists(path):
        return None
    rows = []
    with open(path, encoding="utf-8") as fp:
        for i, line in enumerate(fp):
            rows.append(line.rstrip("\n").split("\t"))
            if limit is not None and i >= limit:
                break
    return rows


def count_lines(path, skip_header=False):
    if not path or not os.path.exists(path):
        return None
    n = 0
    with open(path, encoding="utf-8") as fp:
        for _ in fp:
            n += 1
    return max(0, n - 1) if skip_header else n


def read_targets(path):
    if not path or not os.path.exists(path):
        return []
    out = []
    with open(path, encoding="utf-8") as fp:
        for line in fp:
            t = line.strip()
            if t and not t.startswith("#"):
                out.append(t)
    # 重複除去 (順序保持)
    seen, uniq = set(), []
    for t in out:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    return uniq


def img_data_uri(path):
    """画像を data URI (base64) にして返す。無ければ None。"""
    if not path or not os.path.exists(path):
        return None
    mime, _ = mimetypes.guess_type(path)
    mime = mime or "image/png"
    with open(path, "rb") as fp:
        b64 = base64.b64encode(fp.read()).decode("ascii")
    return f"data:{mime};base64,{b64}"


# ---------------------------------------------------------------------------
# HTML 生成ユーティリティ
# ---------------------------------------------------------------------------
def esc(x):
    return html.escape(str(x))


def fig_block(title, caption, path, embed, metric=None):
    """図 1 枚分の HTML。存在しなければ空文字。

    embed=True なら base64 埋め込み、embed が文字列なら出力 HTML の
    絶対パスとみなし、そこからの相対リンクを参照する。
    metric を渡すと data-metric 属性を付け、HTML 上のボタンで表示切替される
    (metric=None の図はメトリクスに依らず常に表示)。
    """
    if not path or not os.path.exists(path):
        return ""
    if embed is True:
        src = img_data_uri(path)
    else:
        src = os.path.relpath(path, os.path.dirname(embed))
    aid = esc(os.path.splitext(os.path.basename(path))[0])
    dm = f' data-metric="{esc(metric)}"' if metric else ""
    return f"""
    <figure class="fig" id="{aid}"{dm}>
      <a href="{src}" target="_blank" rel="noopener">
        <img src="{src}" alt="{esc(title)}" loading="lazy">
      </a>
      <figcaption><span class="figtitle">{esc(title)}</span>{(' — ' + esc(caption)) if caption else ''}</figcaption>
    </figure>"""


def table_block(rows, max_rows, numeric_cols=()):
    """TSV rows ([header,...]) を HTML テーブルに。max_rows でデータ行を制限。"""
    if not rows:
        return "<p class='muted'>(データなし)</p>"
    header, body = rows[0], rows[1:]
    total = len(body)
    body = body[:max_rows]
    th = "".join(f"<th>{esc(c)}</th>" for c in header)
    trs = []
    for r in body:
        tds = []
        for j, c in enumerate(r):
            cls = ' class="num"' if j in numeric_cols else ""
            tds.append(f"<td{cls}>{esc(c)}</td>")
        trs.append("<tr>" + "".join(tds) + "</tr>")
    note = ""
    if total > len(body):
        note = f"<p class='muted'>上位 {len(body)} 行を表示 (全 {total} 行)</p>"
    return f"""<table class="tbl"><thead><tr>{th}</tr></thead>
<tbody>{''.join(trs)}</tbody></table>{note}"""


# 重要度メトリクス (visualize.py の IMP_COLUMNS[2:] と一致)。図・テーブルの
# ファイル名接尾辞から検出し、HTML 上のボタンで切り替える。
KNOWN_METRICS = ["dlogL", "dBIC", "dK2", "dBDeu",
                 "mean_dlogL_per_sample", "std_dlogL_per_sample"]

# 図の説明 (visualize.py の出力に対応)
FIG_CAPTIONS = {
    "01_structure_full": ("全体構造", "ノードサイズ=次数、赤=注目遺伝子、ラベルはハブ+注目遺伝子"),
    "02_importance_full": ("エッジ重要度 (全体)", "エッジの色・太さ = 重要度 |metric|"),
    "03_importance_top": ("重要度 上位エッジ", "重要度上位の部分グラフ (全ノードにラベル)"),
    "04_targets_highlight": ("注目遺伝子の強調", "赤=注目遺伝子、緑=近傍、橙=関連エッジ"),
    "05_target_ego": ("注目遺伝子の近傍網", "注目遺伝子とその直接の親子のみ"),
    "06_bootstrap_prob": ("ブートストラップ確率", "エッジ出現頻度 = コンセンサス網の安定性"),
}

# メトリクスごとに切り替える (= metric 接尾辞を持つ) 図キー
METRIC_BEARING = {"02_importance_full", "03_importance_top", "05_target_ego"}


def strip_metric(stem):
    """末尾の `_<metric>` を剥がし (base, metric) を返す。metric が無ければ (stem, None)。"""
    for m in sorted(KNOWN_METRICS, key=len, reverse=True):
        if stem.endswith("_" + m):
            return stem[: -(len(m) + 1)], m
    return stem, None


def collect_main_figs(fig_dir, embed, allowed=None):
    """fig_dir 直下の 01〜06 系の図を収集。
    返り値: (blocks, metrics)  metrics は登場したメトリクス集合。
    メトリクス依存図 (02/03/05) は data-metric 付き、それ以外は常時表示。
    同一キーにメトリクス版があれば、接尾辞なしの旧ファイルは除外する。
    allowed が集合なら、そのメトリクスの図だけを含める。
    """
    if not fig_dir or not os.path.isdir(fig_dir):
        return [], set()
    files = sorted(f for f in os.listdir(fig_dir) if f.lower().endswith(".png"))
    # key -> list of (fn, metric)
    grouped = {}
    for fn in files:
        stem = os.path.splitext(fn)[0]
        key = next((k for k in FIG_CAPTIONS if stem == k or stem.startswith(k + "_")), None)
        if key is None:
            continue
        suffix = stem[len(key):].lstrip("_")
        metric = suffix if suffix in KNOWN_METRICS else None
        grouped.setdefault(key, []).append((fn, metric))

    blocks, metrics = [], set()
    for key in FIG_CAPTIONS:  # 図番号順
        entries = grouped.get(key, [])
        if not entries:
            continue
        title, cap = FIG_CAPTIONS[key]
        has_metric = any(m for _, m in entries)
        if has_metric:  # メトリクス版がある → 接尾辞なしの旧ファイルは無視 + allowed で絞る
            use = [(fn, m) for fn, m in entries
                   if m and (allowed is None or m in allowed)]
        else:  # 01/04/06 等は常時表示
            use = [(fn, None) for fn, _ in entries]
        for fn, m in sorted(use):
            if m:
                metrics.add(m)
            blocks.append(fig_block(title, cap, os.path.join(fig_dir, fn), embed, metric=m))
    return blocks, metrics


def collect_subset_figs(fig_dir, embed, allowed=None):
    """群別 (subsets/) の図を収集。返り値: (blocks, metrics)。全て metric 依存。
    allowed が集合なら、そのメトリクスの図だけを含める。"""
    sub = os.path.join(fig_dir, "subsets")
    if not os.path.isdir(sub):
        return [], set()

    def rank(fn):  # 統合図 (overlay/multichannel) → grid → 個別
        if "multichannel" in fn:
            return 0
        if "overlay" in fn:
            return 1
        if "grid" in fn:
            return 2
        return 3

    order = sorted(f for f in os.listdir(sub) if f.lower().endswith(".png"))
    order.sort(key=lambda fn: (rank(fn), fn))
    blocks, metrics = [], set()
    for fn in order:
        stem = os.path.splitext(fn)[0]
        base, metric = strip_metric(stem)
        if metric and allowed is not None and metric not in allowed:
            continue
        if metric:
            metrics.add(metric)
        if "multichannel" in fn:
            title = "群別 統合図 (3チャンネル)"
            cap = "色の種類=群、色の濃さ=重要度、線の太さ=ブートストラップ確率"
        elif "overlay" in fn:
            title = "群別 統合図"
            cap = "全群を1枚に重ねた図 (色の種類=群、色の濃さ=重要度)"
        elif "grid" in fn:
            title = "群別比較 (一覧)"
            cap = "全群を共通レイアウトで並べた図 (群ごとに Blues/Greens/Reds…)"
        else:
            label = base.replace("subset_", "")
            title, cap = f"群: {label}", "全体網を薄く背景表示し、この群で重要なエッジを強調"
        blocks.append(fig_block(title, cap, os.path.join(sub, fn), embed, metric=metric))
    return blocks, metrics


def collect_importance_tables(fig_dir, max_rows, allowed=None):
    """edge_importance_named[_<metric>].tsv を収集し、metric ごとの
    (metric, table_html) リストと metrics 集合を返す。
    allowed が集合なら、そのメトリクスのテーブルだけを含める。"""
    if not fig_dir or not os.path.isdir(fig_dir):
        return [], set()
    tables, metrics = [], set()
    files = sorted(f for f in os.listdir(fig_dir)
                   if f.startswith("edge_importance_named") and f.endswith(".tsv"))
    for fn in files:
        stem = fn[:-4]
        suffix = stem[len("edge_importance_named"):].lstrip("_")
        metric = suffix if suffix in KNOWN_METRICS else None
        if metric and allowed is not None and metric not in allowed:
            continue
        rows = read_tsv(os.path.join(fig_dir, fn))
        if not rows:
            continue
        if metric:
            metrics.add(metric)
        ncol = len(rows[0]) - 1 if rows[0] else 0
        tables.append((metric, table_block(rows, max_rows, numeric_cols=(3, 4, ncol))))
    return tables, metrics


# ---------------------------------------------------------------------------
# メイン
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="BN 解析結果の HTML レポートを生成")
    ap.add_argument("--base-dir", default=".",
                    help="解析ディレクトリ (既定: カレントディレクトリ)")
    ap.add_argument("--out", default="report.html",
                    help="出力 HTML パス (既定: report.html)")
    ap.add_argument("--figures", default="figures",
                    help="学習網の図ディレクトリ (--base-dir からの相対 or 絶対)")
    ap.add_argument("--figures-bs", default="figures_bs",
                    help="コンセンサス網の図ディレクトリ")
    ap.add_argument("--input-tsv", default="data/expr_disc.tsv",
                    help="離散化済み入力 (サンプル数・変数数の要約に使用)")
    ap.add_argument("--var-map", default="data/var_map.tsv",
                    help="変数対応表 var_map.tsv")
    ap.add_argument("--samples", default="data/samples.tsv",
                    help="サンプル表 (row_index / sample_id / group)")
    ap.add_argument("--groups-manifest", default="groups/groups_manifest.tsv",
                    help="群のマニフェスト (make_groups.py の出力)")
    ap.add_argument("--target-file", default="target_genes.txt",
                    help="注目遺伝子 (ホワイトリスト) ファイル")
    ap.add_argument("--edges", default="out/edges.tsv",
                    help="学習網のエッジファイル (エッジ数の要約に使用)")
    ap.add_argument("--integ-edges", default="out/integ_edges2.tsv",
                    help="コンセンサス網のエッジファイル")
    ap.add_argument("--top-edges", type=int, default=40,
                    help="重要度テーブルに載せる上位エッジ数 (既定 40)")
    ap.add_argument("--metrics", default="",
                    help="レポートに含めるメトリクスをカンマ区切りで限定 (既定: 検出した全て)。"
                         f"選択肢: {','.join(KNOWN_METRICS)}。埋め込み版のサイズ削減に有用")
    ap.add_argument("--embed", action="store_true",
                    help="画像を base64 で HTML に埋め込む (単一ファイルで自己完結。"
                         "既定は相対リンク参照で軽量)")
    ap.add_argument("--link-images", action="store_true",
                    help="(既定の挙動) 画像を相対リンク参照にする。互換のため残置")
    ap.add_argument("--title", default="ベイジアンネットワーク解析レポート",
                    help="レポートのタイトル")
    ap.add_argument("--subtitle", default="FastBN / fast_bn パイプライン成果物レポート",
                    help="タイトル下に表示する副題")
    args = ap.parse_args()

    base = args.base_dir
    # 既定は相対リンク参照 (軽量)。--embed 指定時のみ base64 埋め込み。
    embed = args.embed
    # link モードのとき fig_block は出力 HTML の場所を基準に相対パスを作る
    embed_ref = True if embed else os.path.abspath(args.out)

    def p(path):
        """--base-dir 基準でパスを解決する (絶対パスはそのまま)。"""
        return path if os.path.isabs(path) else os.path.join(base, path)

    # --- データ要約 ---------------------------------------------------------
    n_genes = count_lines(p(args.var_map), skip_header=True)
    expr = read_tsv(p(args.input_tsv), limit=0)
    n_samples = count_lines(p(args.input_tsv), skip_header=True)
    # 離散化済み入力はヘッダ=変数名のみ (サンプル名列は持たない)
    n_feat = len(expr[0]) if expr and expr[0] else None
    groups = read_tsv(p(args.groups_manifest))
    if not groups:
        groups = read_tsv(p(args.samples))
    targets = read_targets(p(args.target_file))
    n_edges_hc = count_lines(p(args.edges))
    n_edges_bs = count_lines(p(args.integ_edges))

    summary_rows = []
    if n_samples is not None:
        summary_rows.append(("サンプル数", n_samples))
    if n_feat is not None:
        summary_rows.append(("解析に使用した遺伝子数 (離散化後)", n_feat))
    if n_genes is not None:
        summary_rows.append(("var_map の遺伝子数", n_genes))
    if n_edges_hc is not None:
        summary_rows.append(("学習網のエッジ数 (Hill-Climb)", n_edges_hc))
    if n_edges_bs is not None:
        summary_rows.append(("コンセンサス網のエッジ数 (Bootstrap)", n_edges_bs))
    if targets:
        summary_rows.append(("注目遺伝子 (ホワイトリスト) 数", len(targets)))

    summary_html = "".join(
        f"<tr><th>{esc(k)}</th><td class='num'>{esc(v)}</td></tr>" for k, v in summary_rows
    ) or "<tr><td class='muted'>(要約データが見つかりません)</td></tr>"

    groups_html = ""
    if groups:
        groups_html = "<h3>実験群 / サンプル</h3>" + table_block(groups, 100)

    targets_html = ""
    if targets:
        chips = " ".join(f"<span class='chip'>{esc(t)}</span>" for t in targets)
        targets_html = f"<h3>注目遺伝子 (ホワイトリスト)</h3><div class='chips'>{chips}</div>"

    # --- ネットワーク別セクション ------------------------------------------
    allowed = {m.strip() for m in args.metrics.split(",") if m.strip()} or None
    if allowed:
        bad = allowed - set(KNOWN_METRICS)
        if bad:
            log(f"警告: 未知のメトリクス {sorted(bad)} は無視 (選択肢: {KNOWN_METRICS})")
        allowed = allowed & set(KNOWN_METRICS) or None
    all_metrics = set()

    def network_section(sec_id, name, desc, fig_dir):
        main_figs, m1 = collect_main_figs(fig_dir, embed_ref, allowed)
        subset_figs, m2 = collect_subset_figs(fig_dir, embed_ref, allowed)
        tables, m3 = collect_importance_tables(fig_dir, args.top_edges, allowed)
        all_metrics.update(m1 | m2 | m3)
        if not main_figs and not subset_figs and not tables:
            return "", None
        parts = [f"<section id='{sec_id}'><h2>{esc(name)}</h2><p>{esc(desc)}</p>"]
        if main_figs:
            parts.append("<div class='grid'>" + "".join(main_figs) + "</div>")
        if subset_figs:
            parts.append("<h3>実験群 (サブセット) 別の重要度比較</h3>")
            parts.append("<div class='grid'>" + "".join(subset_figs) + "</div>")
        if tables:
            parts.append(f"<h3>エッジ重要度 上位 {args.top_edges}</h3>")
            for metric, tbl in tables:
                dm = f' data-metric="{esc(metric)}"' if metric else ""
                parts.append(f"<div class='tbl-wrap'{dm}>{tbl}</div>")
        parts.append("</section>")
        return "".join(parts), name

    sections = []
    nav = []
    sec, nm = network_section(
        "hc", "Hill-Climb 学習網",
        "全サンプルで Hill-Climb により推定した DAG とそのエッジ重要度。",
        p(args.figures))
    if sec:
        sections.append(sec)
        nav.append(("hc", nm))
    sec, nm = network_section(
        "bs", "Bootstrap コンセンサス網",
        "ブートストラップ・リサンプリングで安定なエッジを採用したコンセンサス網。"
        "06 図はエッジの出現頻度 (安定性) を表す。",
        p(args.figures_bs))
    if sec:
        sections.append(sec)
        nav.append(("bs", nm))

    if not sections:
        log("図・重要度テーブルが見つかりませんでした。先に viz.sh などを実行してください。")

    # メトリクス切替ボタン (登場したメトリクスを既知の順で)
    metrics_ordered = [m for m in KNOWN_METRICS if m in all_metrics]
    default_metric = metrics_ordered[0] if metrics_ordered else ""
    if metrics_ordered:
        btns = "".join(
            f"<button class='metric-btn' type='button' data-mbtn='{esc(m)}'>{esc(m)}</button>"
            for m in metrics_ordered)
        metricbar_html = (f"<div class='metricbar'><span class='metriclabel'>"
                          f"表示メトリクス:</span>{btns}</div>")
        metric_script = f"""<script>
(function() {{
  function apply(m) {{
    document.querySelectorAll('[data-metric]').forEach(function(el) {{
      el.style.display = (el.getAttribute('data-metric') === m) ? '' : 'none';
    }});
    document.querySelectorAll('.metric-btn').forEach(function(b) {{
      b.classList.toggle('active', b.getAttribute('data-mbtn') === m);
    }});
  }}
  document.querySelectorAll('.metric-btn').forEach(function(b) {{
    b.addEventListener('click', function() {{ apply(b.getAttribute('data-mbtn')); }});
  }});
  apply({default_metric!r});
}})();
</script>"""
        log(f"メトリクス切替ボタン: {metrics_ordered} (既定 {default_metric})")
    else:
        metricbar_html = ""
        metric_script = ""

    nav_html = "".join(
        f"<a href='#{sid}'>{esc(nm)}</a>" for sid, nm in
        [("summary", "概要")] + nav)

    # --- HTML 組み立て ------------------------------------------------------
    doc = f"""<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(args.title)}</title>
<style>
  :root {{
    --bg:#ffffff; --fg:#1a1d21; --muted:#6b7280; --line:#e3e6ea;
    --card:#f7f8fa; --accent:#2b6cb0; --chip:#e6f0fa;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg:#15181c; --fg:#e6e8eb; --muted:#9aa3ad; --line:#2a2f36;
      --card:#1c2026; --accent:#7fb2e6; --chip:#233245;
    }}
  }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; font-family:system-ui,-apple-system,"Segoe UI",
          "Hiragino Sans","Noto Sans JP",sans-serif; color:var(--fg);
          background:var(--bg); line-height:1.6; }}
  header {{ padding:2rem 1.5rem 1.2rem; border-bottom:1px solid var(--line); }}
  header h1 {{ margin:0 0 .3rem; font-size:1.6rem; }}
  header .meta {{ color:var(--muted); font-size:.9rem; }}
  .topbar {{ position:sticky; top:0; background:var(--bg);
             border-bottom:1px solid var(--line); z-index:5; }}
  nav {{ padding:.6rem 1.5rem; display:flex; gap:1rem; flex-wrap:wrap; }}
  nav a {{ color:var(--accent); text-decoration:none; font-size:.92rem; }}
  nav a:hover {{ text-decoration:underline; }}
  .metricbar {{ padding:.5rem 1.5rem; display:flex; gap:.5rem; flex-wrap:wrap;
                align-items:center; border-top:1px solid var(--line); }}
  .metriclabel {{ color:var(--muted); font-size:.85rem; margin-right:.2rem; }}
  .metric-btn {{ cursor:pointer; border:1px solid var(--line); background:var(--card);
                 color:var(--fg); border-radius:999px; padding:.25rem .95rem;
                 font-size:.85rem; font-variant-numeric:tabular-nums; }}
  .metric-btn:hover {{ border-color:var(--accent); }}
  .metric-btn.active {{ background:var(--accent); color:#fff; border-color:var(--accent); }}
  main {{ max-width:1200px; margin:0 auto; padding:1.5rem; }}
  section {{ margin:0 0 3rem; }}
  h2 {{ font-size:1.3rem; border-bottom:2px solid var(--accent); padding-bottom:.3rem;
        margin-top:2.5rem; }}
  h3 {{ font-size:1.05rem; color:var(--muted); margin:1.8rem 0 .8rem; }}
  .grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(340px,1fr));
           gap:1.2rem; }}
  figure.fig {{ margin:0; background:var(--card); border:1px solid var(--line);
                border-radius:10px; padding:.6rem; overflow:hidden; }}
  figure.fig img {{ width:100%; height:auto; border-radius:6px; display:block;
                    background:#fff; }}
  figcaption {{ font-size:.85rem; color:var(--muted); margin-top:.5rem; }}
  .figtitle {{ color:var(--fg); font-weight:600; }}
  table.tbl {{ border-collapse:collapse; width:100%; font-size:.85rem; margin:.5rem 0; }}
  table.tbl th, table.tbl td {{ border:1px solid var(--line); padding:.35rem .55rem;
                                text-align:left; }}
  table.tbl thead th {{ background:var(--card); position:sticky; top:0; }}
  .num {{ text-align:right; font-variant-numeric:tabular-nums; }}
  .tbl-wrap {{ overflow-x:auto; }}
  .kv {{ border-collapse:collapse; }}
  .kv th {{ text-align:left; padding:.35rem .8rem .35rem 0; color:var(--muted);
            font-weight:500; }}
  .kv td {{ padding:.35rem 0; }}
  .chips {{ display:flex; flex-wrap:wrap; gap:.4rem; }}
  .chip {{ background:var(--chip); color:var(--fg); border-radius:999px;
           padding:.15rem .7rem; font-size:.82rem; }}
  .muted {{ color:var(--muted); font-size:.85rem; }}
  footer {{ border-top:1px solid var(--line); padding:1.2rem 1.5rem; color:var(--muted);
            font-size:.82rem; text-align:center; }}
</style>
</head>
<body>
<header>
  <h1>{esc(args.title)}</h1>
  <div class="meta">{esc(args.subtitle)}</div>
</header>
<div class="topbar">
  <nav>{nav_html}</nav>
  {metricbar_html}
</div>
<main>
  <section id="summary">
    <h2>概要</h2>
    <p>バルク RNA 発現量データを離散化し、ベイジアンネットワーク構造を推定・
       エッジ重要度を評価した結果をまとめています。学習網 (Hill-Climb) と
       ブートストラップ・コンセンサス網の 2 種のネットワークを掲載します。</p>
    <table class="kv">{summary_html}</table>
    {groups_html}
    {targets_html}
  </section>
  {"".join(f'<div class="tbl-wrap">{s}</div>' if False else s for s in sections)}
</main>
<footer>Generated by make_report.py — FastBN analysis pipeline</footer>
{metric_script}
</body>
</html>"""

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fp:
        fp.write(doc)
    size_mb = os.path.getsize(args.out) / 1e6
    log(f"出力: {args.out} ({size_mb:.1f} MB, "
        f"{'画像埋め込み' if embed else '画像リンク参照'})")
    log(f"セクション数: {len(sections)}")


if __name__ == "__main__":
    main()

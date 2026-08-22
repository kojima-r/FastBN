#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_benchmark_report.py
========================
構造学習ベンチマークの成果物 (集約表・指標グラフ・正解 vs 学習ネットワークの
比較図) を 1 つの HTML にまとめる。

    python3 make_benchmark_report.py \
        --benchmark results/benchmark.tsv \
        --summary results/summary.tsv \
        --plot results/summary.png \
        --compare-dir results/figures \
        --out results/report.html

既定では画像を相対リンクで参照する (軽量)。--embed を付けると base64 で埋め込み、
HTML 1 ファイルだけで共有できるようになる。

比較図は <compare-dir>/<network>/<何か>.png という配置を想定し、ネットワークごとに
セクションを分けて並べる (フラットに置いてあってもまとめて 1 セクションになる)。
"""

import argparse
import base64
import html
import mimetypes
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def log(*args):
    print("[report]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[report] エラー: {msg}")


def read_tsv(path):
    if not path or not os.path.isfile(path):
        return None, []
    with open(path, encoding="utf-8") as fp:
        rows = [line.rstrip("\n").split("\t") for line in fp if line.strip()]
    if not rows:
        return None, []
    return rows[0], rows[1:]


def esc(text):
    return html.escape(str(text), quote=True)


def img_src(path, out_html, embed):
    if embed:
        mime = mimetypes.guess_type(path)[0] or "image/png"
        with open(path, "rb") as fp:
            data = base64.b64encode(fp.read()).decode("ascii")
        return f"data:{mime};base64,{data}"
    return esc(os.path.relpath(path, os.path.dirname(os.path.abspath(out_html))))


# 表示名 (列名 -> 見出し)
PRETTY = {
    "shd": "SHD",
    "precision_directed": "Precision (dir)",
    "recall_directed": "Recall (dir)",
    "f1_directed": "F1 (dir)",
    "precision_skeleton": "Precision (skel)",
    "recall_skeleton": "Recall (skel)",
    "f1_skeleton": "F1 (skel)",
    "sid": "SID",
    "sid_normalized": "SID (正規化)",
    "kl_divergence": "KL",
    "n_runs": "runs",
    "network": "ネットワーク",
    "n": "サンプル数",
    "score": "スコア",
    "rep": "反復",
}


def pretty(col):
    return PRETTY.get(col, col)


def summary_table_html(header, rows):
    """summarize_benchmark.py の <metric>_mean / _sd を「平均 ± SD」に畳む。"""
    if not header:
        return "<p>集約表がありません。</p>"
    metrics, keys = [], []
    for col in header:
        if col.endswith("_mean"):
            metrics.append(col[:-5])
        elif not col.endswith("_sd") and not col.endswith("_missing"):
            keys.append(col)
    idx = {c: i for i, c in enumerate(header)}

    out = ["<table><thead><tr>"]
    for k in keys:
        out.append(f"<th>{esc(pretty(k))}</th>")
    for m in metrics:
        out.append(f"<th>{esc(pretty(m))}</th>")
    out.append("</tr></thead><tbody>")
    for r in rows:
        out.append("<tr>")
        for k in keys:
            out.append(f"<td>{esc(r[idx[k]])}</td>")
        for m in metrics:
            mean, sd = r[idx[f"{m}_mean"]], r[idx[f"{m}_sd"]]
            cell = "NA" if mean == "NA" else f"{mean} <span class='sd'>± {sd}</span>"
            out.append(f"<td class='num'>{cell}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def collect_compare_figures(compare_dir):
    """<compare-dir> 以下の PNG を (ネットワーク名, [パス]) にまとめる。"""
    groups = {}
    if not compare_dir or not os.path.isdir(compare_dir):
        return groups
    for root, _dirs, files in os.walk(compare_dir):
        pngs = sorted(os.path.join(root, f) for f in files if f.endswith(".png"))
        if not pngs:
            continue
        rel = os.path.relpath(root, compare_dir)
        label = "全体" if rel == "." else rel
        groups.setdefault(label, []).extend(pngs)
    return dict(sorted(groups.items()))


CSS = """
:root { color-scheme: light dark; }
* { box-sizing: border-box; }
body { margin: 0; padding: 2rem 1.25rem 4rem; font-family: system-ui, -apple-system,
       "Hiragino Sans", "Noto Sans JP", Meiryo, sans-serif; line-height: 1.65;
       max-width: 1200px; margin-inline: auto; }
h1 { font-size: 1.6rem; margin: 0 0 .25rem; }
h2 { font-size: 1.2rem; margin: 2.5rem 0 .75rem; padding-bottom: .3rem;
     border-bottom: 2px solid currentColor; }
h3 { font-size: 1rem; margin: 1.5rem 0 .5rem; }
p.lead { margin: 0 0 1.5rem; opacity: .75; }
table { border-collapse: collapse; width: 100%; font-size: .85rem; margin: .5rem 0 1rem; }
th, td { border: 1px solid rgba(128,128,128,.35); padding: .3rem .5rem; text-align: left; }
th { background: rgba(128,128,128,.12); font-weight: 600; white-space: nowrap; }
td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
.sd { opacity: .55; font-size: .85em; }
.scroll { overflow-x: auto; }
figure { margin: 0 0 1.5rem; }
figure img { width: 100%; height: auto; border: 1px solid rgba(128,128,128,.3);
             border-radius: 4px; }
figcaption { font-size: .8rem; opacity: .7; margin-top: .3rem; }
.legend { display: flex; flex-wrap: wrap; gap: 1rem; font-size: .85rem; margin: .5rem 0 1rem; }
.legend span { display: inline-flex; align-items: center; gap: .4rem; }
.swatch { width: 1.6rem; height: 3px; border-radius: 2px; display: inline-block; }
.meta { font-size: .85rem; opacity: .75; }
"""

LEGEND = """
<div class="legend">
  <span><i class="swatch" style="background:#2e7d32"></i>matched (向きまで一致)</span>
  <span><i class="swatch" style="background:#ef6c00"></i>reversed (向きが逆)</span>
  <span><i class="swatch" style="background:#c62828"></i>extra (余分)</span>
  <span><i class="swatch" style="background:#9e9e9e"></i>missing (見落とし)</span>
</div>
"""


def main():
    ap = argparse.ArgumentParser(description="ベンチマーク結果の HTML レポート")
    ap.add_argument("--benchmark", default=None, help="実行ごとの生指標 TSV")
    ap.add_argument("--summary", default=None, help="集約 TSV")
    ap.add_argument("--summary-overall", default=None,
                    help="サンプル数だけで集約した TSV (任意)")
    ap.add_argument("--networks", default=None, help="ネットワーク一覧の TSV (任意)")
    ap.add_argument("--plot", default=None, help="指標グラフ PNG")
    ap.add_argument("--compare-dir", default=None,
                    help="正解 vs 学習の比較図が入ったディレクトリ")
    ap.add_argument("--out", required=True, help="出力 HTML")
    ap.add_argument("--title", default="構造学習ベンチマーク")
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--embed", action="store_true", help="画像を base64 埋め込み")
    args = ap.parse_args()

    parts = [f"<h1>{esc(args.title)}</h1>"]
    if args.subtitle:
        parts.append(f"<p class='lead'>{esc(args.subtitle)}</p>")

    # --- ネットワーク一覧 ---------------------------------------------------
    nh, nr = read_tsv(args.networks)
    if nh:
        parts.append("<h2>対象ネットワーク</h2><div class='scroll'><table><thead><tr>")
        parts += [f"<th>{esc(pretty(c))}</th>" for c in nh]
        parts.append("</tr></thead><tbody>")
        for r in nr:
            parts.append("<tr>" + "".join(f"<td>{esc(c)}</td>" for c in r) + "</tr>")
        parts.append("</tbody></table></div>")

    # --- 集約表 -------------------------------------------------------------
    sh, sr = read_tsv(args.summary_overall)
    if sh:
        parts.append("<h2>サンプル数ごとの精度 (全ネットワーク・全スコアの平均)</h2>")
        parts.append("<div class='scroll'>" + summary_table_html(sh, sr) + "</div>")

    sh, sr = read_tsv(args.summary)
    if sh:
        parts.append("<h2>条件ごとの精度</h2>")
        parts.append("<p class='meta'>各セルは繰り返しにわたる平均 ± 標準偏差。"
                     "SHD・SID・KL は小さいほど、Precision / Recall / F1 は"
                     "大きいほど良い。</p>")
        parts.append("<div class='scroll'>" + summary_table_html(sh, sr) + "</div>")

    # --- グラフ -------------------------------------------------------------
    if args.plot and os.path.isfile(args.plot):
        parts.append("<h2>サンプル数に対する推移</h2>")
        parts.append(f"<figure><img src='{img_src(args.plot, args.out, args.embed)}' "
                     "alt='metrics vs sample size'>"
                     "<figcaption>横軸 = サンプル数 (対数)、線 = スコア関数、"
                     "列 = ネットワーク。</figcaption></figure>")

    # --- 比較図 -------------------------------------------------------------
    groups = collect_compare_figures(args.compare_dir)
    if groups:
        parts.append("<h2>正解ネットワークと学習ネットワークの比較</h2>")
        parts.append("<p class='meta'>左が正解、右が学習結果。ノードの配置は正解 DAG "
                     "から決めているので、両パネルおよび条件をまたいで同じ位置に"
                     "描かれる。</p>")
        parts.append(LEGEND)
        for label, pngs in groups.items():
            parts.append(f"<h3>{esc(label)}</h3>")
            for png in pngs:
                cap = os.path.splitext(os.path.basename(png))[0]
                parts.append(
                    f"<figure><img src='{img_src(png, args.out, args.embed)}' "
                    f"alt='{esc(cap)}'><figcaption>{esc(cap)}</figcaption></figure>")

    # --- 生データへの案内 ---------------------------------------------------
    bh, br = read_tsv(args.benchmark)
    if bh:
        parts.append("<h2>生データ</h2>")
        rel = os.path.relpath(os.path.abspath(args.benchmark),
                              os.path.dirname(os.path.abspath(args.out)))
        parts.append(f"<p class='meta'>実行ごとの全指標 ({len(br)} 行): "
                     f"<code>{esc(rel)}</code></p>")

    doc = ("<!doctype html><html lang='ja'><head><meta charset='utf-8'>"
           "<meta name='viewport' content='width=device-width,initial-scale=1'>"
           f"<title>{esc(args.title)}</title><style>{CSS}</style></head><body>"
           + "".join(parts) + "</body></html>")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fp:
        fp.write(doc)
    size = os.path.getsize(args.out) / 1e6
    log(f"出力: {args.out} ({size:.1f} MB, "
        f"{'埋め込み' if args.embed else '画像リンク参照'})")


if __name__ == "__main__":
    main()

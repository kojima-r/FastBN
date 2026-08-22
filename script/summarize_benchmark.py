#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
summarize_benchmark.py
======================
evaluate_structure.py が 1 行ずつ追記したベンチマーク結果 (TSV) を、指定した
キーで集約して平均 ± 標準偏差の表にする。Markdown 表と、必要なら折れ線グラフも
出力する。

    python3 summarize_benchmark.py \
        --input results/benchmark.tsv \
        --group-by network,n,score \
        --metrics shd,precision_directed,recall_directed,f1_directed,kl_divergence \
        --plot-metrics shd,f1_directed,kl_divergence \
        --out results/summary.tsv --markdown results/summary.md \
        --plot results/summary.png --plot-x n --plot-facet network --plot-series score

数値に変換できない値 (NA / inf) は平均から除外し、除外件数を n_missing 列に出す。
"""

import argparse
import math
import os
import sys

import numpy as np


def log(*args):
    print("[summary]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[summary] エラー: {msg}")


def read_tsv(path):
    if not os.path.isfile(path):
        die(f"{path} がありません")
    with open(path, encoding="utf-8") as fp:
        rows = [line.rstrip("\n").split("\t") for line in fp if line.strip()]
    if len(rows) < 2:
        die(f"{path} にデータ行がありません")
    header, body = rows[0], rows[1:]
    return header, [dict(zip(header, r)) for r in body]


def to_float(text):
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(v) or math.isinf(v)) else v


def sort_key(values):
    """数値に見えるものは数値として並べる。"""
    out = []
    for v in values:
        try:
            out.append((0, float(v), ""))
        except (TypeError, ValueError):
            out.append((1, 0.0, str(v)))
    return tuple(out)


def main():
    ap = argparse.ArgumentParser(description="ベンチマーク結果の集約")
    ap.add_argument("--input", required=True, help="evaluate_structure.py の出力 TSV")
    ap.add_argument("--group-by", default="network,n,score",
                    help="集約キー (カンマ区切り)")
    ap.add_argument("--metrics",
                    default="shd,precision_directed,recall_directed,f1_directed,"
                            "precision_skeleton,recall_skeleton,f1_skeleton,"
                            "sid_normalized,kl_divergence",
                    help="集約する指標 (カンマ区切り)")
    ap.add_argument("--plot-metrics", default=None,
                    help="グラフに描く指標 (カンマ区切り)。既定は --metrics と同じ。"
                         "表には P/R/F1 を全部入れつつ、図は絞りたいときに使う")
    ap.add_argument("--out", default=None, help="集約結果の TSV")
    ap.add_argument("--markdown", default=None, help="集約結果の Markdown 表")
    ap.add_argument("--plot", default=None, help="折れ線グラフの PNG")
    ap.add_argument("--plot-x", default="n", help="グラフの横軸に使う列")
    ap.add_argument("--plot-facet", default="network", help="サブプロットを分ける列")
    ap.add_argument("--plot-series", default="score", help="線を分ける列")
    ap.add_argument("--precision", type=int, default=3)
    args = ap.parse_args()

    header, rows = read_tsv(args.input)
    keys = [k for k in args.group_by.split(",") if k]
    metrics = [m for m in args.metrics.split(",") if m]
    for col in keys + metrics:
        if col not in header:
            die(f"列 '{col}' が {args.input} にありません (列: {', '.join(header)})")

    groups = {}
    for r in rows:
        gk = tuple(r[k] for k in keys)
        groups.setdefault(gk, []).append(r)

    out_header = keys + ["n_runs"]
    for m in metrics:
        out_header += [f"{m}_mean", f"{m}_sd", f"{m}_missing"]

    out_rows = []
    for gk in sorted(groups, key=sort_key):
        members = groups[gk]
        row = list(gk) + [str(len(members))]
        for m in metrics:
            vals = [to_float(r.get(m)) for r in members]
            good = [v for v in vals if v is not None]
            missing = len(vals) - len(good)
            if good:
                mean = float(np.mean(good))
                sd = float(np.std(good, ddof=1)) if len(good) > 1 else 0.0
                row += [f"{mean:.{args.precision}f}", f"{sd:.{args.precision}f}",
                        str(missing)]
            else:
                row += ["NA", "NA", str(missing)]
        out_rows.append(row)

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fp:
            fp.write("\t".join(out_header) + "\n")
            for r in out_rows:
                fp.write("\t".join(r) + "\n")
        log(f"集約結果: {args.out} ({len(out_rows)} 行)")

    md_header = keys + ["runs"] + metrics
    md_lines = ["| " + " | ".join(md_header) + " |",
                "| " + " | ".join(["---"] * len(md_header)) + " |"]
    for r in out_rows:
        cells = list(r[:len(keys) + 1])
        for i, _ in enumerate(metrics):
            mean, sd = r[len(keys) + 1 + 3 * i], r[len(keys) + 2 + 3 * i]
            cells.append(f"{mean} ± {sd}" if mean != "NA" else "NA")
        md_lines.append("| " + " | ".join(cells) + " |")
    md = "\n".join(md_lines)
    if args.markdown:
        os.makedirs(os.path.dirname(os.path.abspath(args.markdown)) or ".",
                    exist_ok=True)
        with open(args.markdown, "w", encoding="utf-8") as fp:
            fp.write(md + "\n")
        log(f"Markdown 表: {args.markdown}")
    if not args.out and not args.markdown:
        print(md)

    if args.plot:
        plot_metrics = ([m for m in args.plot_metrics.split(",") if m]
                        if args.plot_metrics else metrics)
        missing = [m for m in plot_metrics if m not in header]
        if missing:
            die(f"--plot-metrics の列が {args.input} にありません: {', '.join(missing)}")
        make_plot(args, rows, plot_metrics)


def make_plot(args, rows, metrics):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log("matplotlib が無いのでグラフはスキップします")
        return

    xcol, fcol, scol = args.plot_x, args.plot_facet, args.plot_series
    for col in (xcol, fcol, scol):
        if col not in rows[0]:
            log(f"列 '{col}' が無いのでグラフはスキップします")
            return

    facets = sorted({r[fcol] for r in rows}, key=lambda v: sort_key([v]))
    series = sorted({r[scol] for r in rows}, key=lambda v: sort_key([v]))
    n_row, n_col = len(metrics), len(facets)
    fig, axes = plt.subplots(n_row, n_col, figsize=(3.0 * n_col, 2.5 * n_row),
                             squeeze=False)

    for mi, metric in enumerate(metrics):
        for fi, facet in enumerate(facets):
            ax = axes[mi][fi]
            for s in series:
                pts = {}
                for r in rows:
                    if r[fcol] != facet or r[scol] != s:
                        continue
                    v = to_float(r.get(metric))
                    if v is None:
                        continue
                    pts.setdefault(r[xcol], []).append(v)
                if not pts:
                    continue
                xs = sorted(pts, key=lambda v: sort_key([v]))
                ys = [float(np.mean(pts[x])) for x in xs]
                ax.plot([float(x) for x in xs], ys, marker="o", markersize=3,
                        linewidth=1.2, label=s)
            ax.set_xscale("log")
            if mi == 0:
                ax.set_title(facet, fontsize=10)
            if fi == 0:
                ax.set_ylabel(metric, fontsize=8)
            if mi == n_row - 1:
                ax.set_xlabel(xcol, fontsize=8)
            ax.tick_params(labelsize=7)
            ax.grid(alpha=0.3, linewidth=0.5)
    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=len(labels),
                   fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    os.makedirs(os.path.dirname(os.path.abspath(args.plot)) or ".", exist_ok=True)
    fig.savefig(args.plot, dpi=150)
    plt.close(fig)
    log(f"グラフ: {args.plot}")


if __name__ == "__main__":
    main()

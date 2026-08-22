#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_dag_comparison.py
======================
正解 DAG と学習された DAG を**同じノード配置**で並べて描き、どのエッジが当たって
いてどれを間違えたかを一目で分かるようにする。

左: 正解ネットワーク / 右: 学習ネットワーク。右のエッジは判定で色分けする。

    緑 実線   正解と向きまで一致          (matched)
    橙 実線   骨格は合っているが向きが逆  (reversed)
    赤 実線   正解に無い余分なエッジ      (extra / false positive)
    灰 破線   学習が見落としたエッジ      (missing / false negative)

図中の文字は ../script/visualize.py と同じく英語にしてある (日本語フォントが
無い環境でも文字化けしないため)。

ノードが多すぎる場合は --max-nodes で次数上位のノードだけを取り出した部分グラフを
描く (副題の指標はネットワーク全体に対する値のまま)。

ノード配置は**正解 DAG から決める**ので、同じネットワークの図はスコア関数や
サンプル数が違っても同じ配置になり、直接見比べられる。配置は正解 DAG の
トポロジカルな階層 (親から子へ上から下) を基本とし、階層化できない場合は
バネモデルにフォールバックする。

使い方:
    python3 plot_dag_comparison.py \
        --true-bif asia.bif \
        --pred-edges out/edges.tsv \
        --input data.tsv \
        --out figures/asia_bic.png \
        --title "asia / bic / n=1000"

正解は --true-bif (BIF) でも --true-edges (エッジ表) でもよい。--input を与えると
図の副題に SHD / F1 / SID / KL などの指標も入る (KL は --true-bif が必要)。
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from evaluate_structure import (build_parents, fmt, is_dag, kl_divergence,
                                read_data, read_edge_file, structural_metrics,
                                structural_intervention_distance)

COLOR_CORRECT = "#2e7d32"   # 緑: 向きまで一致
COLOR_REVERSED = "#ef6c00"  # 橙: 向きが逆
COLOR_EXTRA = "#c62828"     # 赤: 余分
COLOR_MISSING = "#9e9e9e"   # 灰: 見落とし
COLOR_TRUE = "#37474f"      # 正解パネルのエッジ


def log(*args):
    print("[plot_dag]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[plot_dag] エラー: {msg}")


def hierarchy_layout(nodes, parents, seed=0):
    """正解 DAG のトポロジカルな深さで階層配置する (親が上、子が下)。"""
    depth = {}

    def calc(v, guard):
        if v in depth:
            return depth[v]
        if v in guard:
            return 0
        guard.add(v)
        d = 0 if not parents[v] else 1 + max(calc(p, guard) for p in parents[v])
        guard.discard(v)
        depth[v] = d
        return d

    for v in nodes:
        calc(v, set())

    levels = {}
    for v in nodes:
        levels.setdefault(depth[v], []).append(v)
    max_w = max(len(vs) for vs in levels.values())
    max_d = max(levels)

    pos = {}
    for d, vs in levels.items():
        vs = sorted(vs)
        for i, v in enumerate(vs):
            # 各段を中央揃えで横に並べる
            x = (i + 0.5) / len(vs) * max_w
            y = -d
            pos[v] = (x, y)
    # 単一段しか無い場合はバネ配置に任せる
    if max_d == 0:
        g = nx.Graph()
        g.add_nodes_from(nodes)
        for v in nodes:
            for p in parents[v]:
                g.add_edge(p, v)
        return nx.spring_layout(g, seed=seed)
    return pos


def classify_edges(true_edges, pred_edges):
    T, P = set(true_edges), set(pred_edges)
    skelT = {frozenset(e) for e in T}
    correct = sorted(T & P)
    reversed_ = sorted(e for e in P if e not in T and (e[1], e[0]) in T)
    extra = sorted(e for e in P if frozenset(e) not in skelT)
    missing = sorted(e for e in T if frozenset(e) not in {frozenset(x) for x in P})
    return correct, reversed_, extra, missing


def draw_panel(ax, nodes, pos, groups, title, node_size, font_size):
    g = nx.DiGraph()
    g.add_nodes_from(nodes)
    nx.draw_networkx_nodes(g, pos, ax=ax, node_size=node_size,
                           node_color="#ffffff", edgecolors="#455a64",
                           linewidths=1.2)
    nx.draw_networkx_labels(g, pos, ax=ax, font_size=font_size)
    for edges, color, style, width in groups:
        if not edges:
            continue
        h = nx.DiGraph()
        h.add_nodes_from(nodes)
        h.add_edges_from(edges)
        nx.draw_networkx_edges(
            h, pos, ax=ax, edgelist=edges, edge_color=color, style=style,
            width=width, arrowsize=11, arrowstyle="-|>",
            node_size=node_size, connectionstyle="arc3,rad=0.06")
    ax.set_title(title, fontsize=11)
    ax.set_axis_off()
    ax.margins(0.12)


def main():
    ap = argparse.ArgumentParser(
        description="正解 DAG と学習 DAG を同じ配置で並べて描く")
    ap.add_argument("--true-bif", default=None, help="正解ネットワーク (BIF)")
    ap.add_argument("--true-edges", default=None, help="正解エッジ表")
    ap.add_argument("--pred-edges", required=True, help="学習エッジ表")
    ap.add_argument("--input", default=None,
                    help="学習に使ったデータ (ノード名の復元と KL に使う)")
    ap.add_argument("--out", required=True, help="出力 PNG")
    ap.add_argument("--title", default="", help="図全体のタイトル")
    ap.add_argument("--alpha", type=float, default=1.0, help="KL 用の平滑化")
    ap.add_argument("--max-states", type=int, default=2_000_000)
    ap.add_argument("--no-metrics", action="store_true",
                    help="副題に指標を入れない")
    ap.add_argument("--max-nodes", type=int, default=0,
                    help="描画するノード数の上限 (0 = 制限なし)。これを超える場合は "
                         "次数の大きいノードだけを取り出した部分グラフを描く。"
                         "副題の指標は常にネットワーク全体に対する値")
    ap.add_argument("--figsize", default="11,5", help="図のサイズ (幅,高さ)")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    if not args.true_bif and not args.true_edges:
        die("--true-bif か --true-edges が必要です")

    true_bn, node_names, X = None, None, None
    if args.true_bif:
        from bif_io import read_bif
        true_bn = read_bif(args.true_bif)
        node_names = list(true_bn.names)
    if args.input:
        header, X = read_data(args.input, node_names)
        if node_names is None:
            node_names = header

    true_edges = (true_bn.edges() if true_bn
                  else read_edge_file(args.true_edges, node_names))
    if node_names is None:
        node_names = sorted({v for e in true_edges for v in e})
    pred_edges = read_edge_file(args.pred_edges, node_names)

    nodes = list(node_names)
    true_parents = build_parents(nodes, true_edges)
    pred_parents = build_parents(nodes, pred_edges)

    correct, reversed_, extra, missing = classify_edges(true_edges, pred_edges)

    # --- 指標 (副題用) ------------------------------------------------------
    subtitle = ""
    if not args.no_metrics:
        m = structural_metrics(nodes, true_edges, pred_edges)
        parts = [f"SHD={m['shd']}",
                 f"F1(dir)={m['f1_directed']:.2f}",
                 f"F1(skel)={m['f1_skeleton']:.2f}"]
        if is_dag(nodes, pred_parents) and is_dag(nodes, true_parents):
            sid = structural_intervention_distance(nodes, true_parents, pred_parents)
            parts.append(f"SID={sid}/{len(nodes) * (len(nodes) - 1)}")
        if true_bn is not None and X is not None \
                and true_bn.state_space_size() <= args.max_states:
            parts.append(f"KL={fmt(kl_divergence(true_bn, pred_parents, X, args.alpha))}")
        subtitle = "   ".join(parts)

    # --- 大きすぎる場合はハブ部分グラフに絞る -------------------------------
    # 指標は上でネットワーク全体に対して計算済み。ここで削るのは描画だけ。
    note = ""
    if args.max_nodes and len(nodes) > args.max_nodes:
        deg = {v: 0 for v in nodes}
        for u, v in list(true_edges) + list(pred_edges):
            deg[u] += 1
            deg[v] += 1
        keep = set(sorted(nodes, key=lambda v: (-deg[v], v))[:args.max_nodes])
        sub = lambda es: [e for e in es if e[0] in keep and e[1] in keep]
        full_n, full_true = len(nodes), len(true_edges)
        nodes = [v for v in nodes if v in keep]
        true_edges = sub(true_edges)
        correct, reversed_, extra, missing = (sub(correct), sub(reversed_),
                                              sub(extra), sub(missing))
        true_parents = build_parents(nodes, true_edges)
        pred_edges = correct + reversed_ + extra
        note = (f"top-{len(nodes)} hub subgraph of {full_n} nodes "
                f"(metrics above are for the full network)")
        log(f"ノードが多いので次数上位 {len(nodes)} / {full_n} の部分グラフを描画します")

    # --- 描画 ---------------------------------------------------------------
    pos = hierarchy_layout(nodes, true_parents)
    n = len(nodes)
    node_size = max(300, min(1400, int(9000 / max(1, n))))
    font_size = max(6, min(10, int(70 / max(1, n)) + 5))
    w, h = (float(x) for x in args.figsize.split(","))
    fig, axes = plt.subplots(1, 2, figsize=(w, h))

    draw_panel(axes[0], nodes, pos,
               [(sorted(true_edges), COLOR_TRUE, "solid", 1.8)],
               f"true network ({len(true_edges)} edges)", node_size, font_size)
    draw_panel(axes[1], nodes, pos,
               [(missing, COLOR_MISSING, "dashed", 1.2),
                (correct, COLOR_CORRECT, "solid", 2.0),
                (reversed_, COLOR_REVERSED, "solid", 2.0),
                (extra, COLOR_EXTRA, "solid", 1.6)],
               f"learned network ({len(pred_edges)} edges)", node_size, font_size)

    handles = [
        Line2D([], [], color=COLOR_CORRECT, lw=2.0,
               label=f"matched ({len(correct)})"),
        Line2D([], [], color=COLOR_REVERSED, lw=2.0,
               label=f"reversed ({len(reversed_)})"),
        Line2D([], [], color=COLOR_EXTRA, lw=1.6,
               label=f"extra ({len(extra)})"),
        Line2D([], [], color=COLOR_MISSING, lw=1.2, ls="--",
               label=f"missing ({len(missing)})"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9,
               frameon=False)

    if args.title or subtitle or note:
        head = args.title
        if subtitle:
            head = f"{head}\n{subtitle}" if head else subtitle
        if note:
            head = f"{head}\n{note}" if head else note
        fig.suptitle(head, fontsize=12)
    fig.tight_layout(rect=(0, 0.07, 1, 0.93 if (args.title or subtitle) else 1.0))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)
    log(f"{args.out} (一致 {len(correct)} / 逆 {len(reversed_)} / "
        f"余分 {len(extra)} / 見落とし {len(missing)})")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
viz_subsets.py
==============
サブセット (実験群) ごとのエッジ重要度を、
  ・全体グラフを背景に薄く表示
  ・そのサブセットで重要なエッジ (とノード) のみ強調
という形で描画する。**ノード配置 (レイアウト) は全サブセットで共通**にするため、
全体ネットワークから 1 度だけレイアウトを計算して全パネルで使い回す。

**配色**: 群ごとに単色系カラーマップ (Blues / Greens / Reds / Purples / Oranges)
を割り当て、色の濃さで重要度を表す。全サンプル参照パネルは Greys。

対象は `<prefix>_g{N}_<label>.tsv` 形式の群別重要度ファイル群
(importance_groups.sh の出力)。全サンプル版 `<prefix>.tsv` があれば
参照パネルとして先頭に含める。

既定のパスは**カレントディレクトリ基準** (./out, ./figures/subsets,
./data/var_map.tsv, ./target_genes.txt) なので、任意の解析ディレクトリで
実行できる。

出力 (--fig-dir):
  subset_<label>_<metric>.png         : サブセットごとの個別図 (共通レイアウト, 群色)
  subsets_grid_<metric>.png           : 全サブセットを並べた比較図 (共通レイアウト)
  subsets_overlay_<metric>.png        : 全群を 1 枚に重ねた統合図
                                        (色の種類=群, 色の濃さ=重要度)
  subsets_multichannel_<metric>.png   : 追加の統合図
                                        (色の種類=群, 色の濃さ=重要度,
                                         線の太さ=ブートストラップ確率)

使い方の簡便化のため viz_subsets.sh / viz_bs_subsets.sh を用意している。
"""

import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

# 同ディレクトリの visualize.py からヘルパーと配色を再利用
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import visualize as V  # noqa: E402


# 群に割り当てる単色系カラーマップ (色の種類で群を区別)
GROUP_PALETTE = ["Blues", "Greens", "Reds", "Purples", "Oranges", "BuPu", "YlOrBr"]
ALL_CMAP = "Greys"  # 全サンプル参照パネル用

# 線幅の倍率: ノード数に応じた倍率 (V.STYLE["edge"]) をそのまま掛けると
# 強調エッジが太すぎて重なるため、控えめ (0.6 乗) に効かせる。
EDGE_W = 1.0


def set_edge_width_scale():
    global EDGE_W
    EDGE_W = float(V.STYLE["edge"]) ** 0.6


def discover_subsets(out_dir, prefix, include_all):
    """`<prefix>_g{N}_<label>.tsv` を探し、(label, path) のリストを返す。
    include_all が真で `<prefix>.tsv` があれば先頭に ("ALL_samples", path) を追加。
    """
    subsets = []
    all_path = os.path.join(out_dir, f"{prefix}.tsv")
    if include_all and os.path.exists(all_path):
        subsets.append(("ALL_samples", all_path))
    pat = re.compile(re.escape(prefix) + r"_g(\d+)_(.+)\.tsv$")
    found = []
    for p in glob.glob(os.path.join(out_dir, f"{prefix}_g*_*.tsv")):
        m = pat.search(os.path.basename(p))
        if m:
            found.append((int(m.group(1)), m.group(2), p))
    found.sort()
    subsets.extend([(label, p) for _, label, p in found])
    return subsets


def assign_cmaps(labels):
    """群ラベル -> 単色系カラーマップ名 の割り当て。
    ALL_samples は Greys、それ以外は GROUP_PALETTE を順に循環させる。"""
    out = {}
    gi = 0
    for lab in labels:
        if lab == "ALL_samples":
            out[lab] = ALL_CMAP
        else:
            out[lab] = GROUP_PALETTE[gi % len(GROUP_PALETTE)]
            gi += 1
    return out


def subset_edge_values(G, imp):
    """グラフの各エッジに対するサブセット重要度値の配列を返す (順序は G.edges())。"""
    return np.array([imp.get(G[u][v].get("_uv", (None, None)), 0.0)
                     for u, v in G.edges()])


def edge_prob_values(G, prob):
    """グラフの各エッジに対するブートストラップ確率の配列 (順序は G.edges())。"""
    return np.array([prob.get(G[u][v].get("_uv", (None, None)), 0.0)
                     for u, v in G.edges()])


def draw_background(ax, G, pos):
    """全体グラフを薄く背景に描く (全図で共通の見た目)。"""
    ax.set_axis_off()
    V.nx.draw_networkx_edges(G, pos, ax=ax, edgelist=list(G.edges()),
                             edge_color="#b8bcc2", width=0.3 * V.STYLE["edge"],
                             arrowsize=3, alpha=0.12)
    V.nx.draw_networkx_nodes(G, pos, ax=ax, node_size=8 * V.STYLE["node"],
                             node_color=V.COL_BG, linewidths=0, alpha=0.35)


def top_indices(vals, top_n):
    order = np.argsort(vals)[::-1]
    return [i for i in order if vals[i] > 0][:top_n]


def draw_subset(ax, G, pos, vals, targets, top_n, vmax, cmap, title):
    """1 サブセット分を ax に描画。背景=全体グラフ薄く、強調=重要上位 top_n エッジ。"""
    edges = list(G.edges())
    draw_background(ax, G, pos)
    ax.set_title(title, fontsize=12)

    top_idx = top_indices(vals, top_n)
    if not top_idx:
        return
    hi_edges = [edges[i] for i in top_idx]
    hi_vals = np.array([vals[i] for i in top_idx])
    wn = hi_vals / vmax if vmax > 0 else hi_vals
    V.nx.draw_networkx_edges(G, pos, ax=ax, edgelist=hi_edges,
                             edge_color=hi_vals, edge_cmap=cmap,
                             width=(1.0 + 4.0 * wn) * EDGE_W,
                             arrowsize=9 * V.STYLE["arrow"], alpha=0.85,
                             edge_vmin=0.0, edge_vmax=vmax)

    # 強調エッジの端点ノードのみ前面に + ラベル
    hi_nodes = set()
    for u, v in hi_edges:
        hi_nodes.add(u); hi_nodes.add(v)
    tgt = [n for n in hi_nodes if V.is_target(n, targets)]
    oth = [n for n in hi_nodes if not V.is_target(n, targets)]
    V.nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=oth,
                             node_size=120 * V.STYLE["node"],
                             node_color=V.COL_NODE, linewidths=0, alpha=V.NODE_ALPHA)
    if tgt:
        V.nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=tgt,
                                 node_size=220 * V.STYLE["node"],
                                 node_color=V.COL_TARGET, linewidths=0,
                                 alpha=V.NODE_ALPHA)
    V.nx.draw_networkx_labels(G, pos, ax=ax,
                              labels={n: n for n in hi_nodes},
                              font_size=V.STYLE["label"], bbox=V.LABEL_BBOX)


def draw_overlay(ax, G, pos, group_vals, cmaps, targets, top_n, vmax,
                 title, prob_vals=None):
    """全群の重要エッジを 1 枚に重ねた統合図。
      色の種類 = 群 (cmaps)、色の濃さ = 重要度。
    prob_vals (idx 対応のエッジ確率配列) を渡すと、線の太さ = ブートストラップ確率
    とする (multichannel 図)。渡さなければ太さは重要度に応じる。
    group_vals: [(label, vals_array), ...] (ALL_samples は除いておくこと)
    """
    edges = list(G.edges())
    draw_background(ax, G, pos)
    ax.set_title(title, fontsize=13)

    labeled = set()
    for label, vals in group_vals:
        cmap = V.make_cmap(cmaps[label])
        top_idx = top_indices(vals, top_n)
        if not top_idx:
            continue
        hi_edges = [edges[i] for i in top_idx]
        hi_vals = np.array([vals[i] for i in top_idx])
        if prob_vals is not None:
            # 線の太さ = ブートストラップ確率 (0..1)
            pw = np.array([prob_vals[i] for i in top_idx])
            width = 0.6 + 4.5 * pw
        else:
            width = 1.0 + 3.5 * (hi_vals / vmax if vmax > 0 else hi_vals)
        V.nx.draw_networkx_edges(G, pos, ax=ax, edgelist=hi_edges,
                                 edge_color=hi_vals, edge_cmap=cmap,
                                 width=width * EDGE_W,
                                 arrowsize=8 * V.STYLE["arrow"], alpha=0.75,
                                 edge_vmin=0.0, edge_vmax=vmax)
        labeled |= {n for e in hi_edges for n in e}

    # ノード: 強調エッジ端点を淡色で、注目遺伝子は赤で強調
    tgt = [n for n in labeled if V.is_target(n, targets)]
    oth = [n for n in labeled if not V.is_target(n, targets)]
    V.nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=oth,
                             node_size=90 * V.STYLE["node"],
                             node_color=V.COL_NODE, linewidths=0, alpha=0.6)
    if tgt:
        V.nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=tgt,
                                 node_size=200 * V.STYLE["node"],
                                 node_color=V.COL_TARGET, linewidths=0,
                                 alpha=V.NODE_ALPHA)
    # 強調エッジの端点ノードすべてに遺伝子名ラベルを表示
    V.nx.draw_networkx_labels(G, pos, ax=ax, labels={n: n for n in labeled},
                              font_size=V.STYLE["label"], bbox=V.LABEL_BBOX)

    # 凡例: 群 -> 代表色
    handles = [Patch(facecolor=V.make_cmap(cmaps[label])(0.85), edgecolor="none",
                     label=label) for label, _ in group_vals]
    ax.legend(handles=handles, loc="upper right", fontsize=9,
              title="group (hue)", framealpha=0.85)


def main():
    ap = argparse.ArgumentParser(
        description="サブセット別 重要度可視化 (背景薄表示+重要部分強調, 共通レイアウト, 群色)")
    ap.add_argument("--out-dir", default="./out")
    ap.add_argument("--edges", default="edges.tsv",
                    help="背景/レイアウト用エッジ (インデックス)")
    ap.add_argument("--edges-named", default="edges_named.tsv",
                    help="背景/レイアウト用エッジ (遺伝子名, --edges と行対応)")
    ap.add_argument("--prefix", default="edge_importance",
                    help="群別重要度ファイルの接頭辞 (例 edge_importance / integ_edge_importance)")
    ap.add_argument("--edge-prob", default=None,
                    help="ブートストラップ確率ファイル (u v count prob; 例 integ_edges_score.tsv)。"
                         "指定すると multichannel 図の線の太さに用いる")
    ap.add_argument("--include-all", action="store_true",
                    help="全サンプル版 <prefix>.tsv も参照パネルとして含める")
    ap.add_argument("--metrics", default="dlogL,dBIC,dK2,dBDeu",
                    help="可視化するスコア (カンマ区切り)。"
                         f"選択肢: {','.join(V.IMP_COLUMNS[2:])} (既定 dlogL,dBIC,dK2,dBDeu)")
    ap.add_argument("--top-n", type=int, default=40,
                    help="各サブセットで強調する上位エッジ数 (既定 40)")
    ap.add_argument("--fig-dir", default=os.path.join("figures", "subsets"))
    ap.add_argument("--target-file", default="./target_genes.txt")
    ap.add_argument("--var-map", dest="var_map", default="./data/var_map.tsv")
    ap.add_argument("--layout", default="spring",
                    choices=["spring", "kamada", "circular"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--style-scale", type=float, default=None,
                    help="ノード・線・文字の大きさの倍率 (既定はノード数から自動決定)")
    args = ap.parse_args()

    def resolve(p):
        return p if os.path.isabs(p) else os.path.join(args.out_dir, p)

    edges_path = resolve(args.edges)
    named_path = resolve(args.edges_named)
    for p in (edges_path, named_path):
        if not os.path.exists(p):
            sys.exit(f"[viz_subsets] 必須ファイルがありません: {p}")

    idx2name = V.load_node_names(edges_path, named_path)
    edges = V.load_edges(edges_path)
    targets = V.load_targets(args.target_file, args.var_map)
    G = V.build_graph(edges, idx2name, imp={})  # _uv を各エッジに保持
    V.set_style(G.number_of_nodes(), args.style_scale)
    set_edge_width_scale()

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    bad = [m for m in metrics if m not in V.IMP_COLUMNS[2:]]
    if bad:
        sys.exit(f"[viz_subsets] 未知のスコア: {bad} (選択肢: {V.IMP_COLUMNS[2:]})")

    subsets = discover_subsets(args.out_dir, args.prefix, args.include_all)
    if not subsets:
        sys.exit(f"[viz_subsets] {args.prefix}_g*_*.tsv が {args.out_dir} に見つかりません。"
                 " 先に run_importance_groups.sh を実行してください。")
    labels = [l for l, _ in subsets]
    cmaps = assign_cmaps(labels)
    V.log(f"サブセット {len(subsets)} 件: {labels}, スコア={metrics}")
    V.log(f"群色: " + ", ".join(f"{l}={cmaps[l]}" for l in labels))

    # ブートストラップ確率 (multichannel 図の線の太さ用)
    prob = {}
    if args.edge_prob:
        pp = resolve(args.edge_prob)
        if os.path.exists(pp):
            prob = V.load_edge_prob(pp)
            V.log(f"ブートストラップ確率 {len(prob)} 件を {pp} から読込")
        else:
            V.log(f"警告: --edge-prob {pp} が見つかりません (太さは重要度で代用)")
    prob_edge_vals = edge_prob_values(G, prob) if prob else None

    os.makedirs(args.fig_dir, exist_ok=True)

    # --- 共通レイアウトを 1 度だけ計算 (全サブセット・全スコアで使い回す) ---
    V.log(f"共通レイアウト計算中 ({args.layout}) ...")
    if args.layout == "spring":
        pos = V.nx.spring_layout(G, seed=args.seed,
                                 k=1.0 / np.sqrt(max(1, G.number_of_nodes())))
    elif args.layout == "kamada":
        pos = V.nx.kamada_kawai_layout(G)
    else:
        pos = V.nx.circular_layout(G)

    # --- スコアごとに描画 (ファイル名に _<metric> を付与) ---
    for metric in metrics:
        # 各サブセットの重要度をロードし、グラフエッジ順の値配列に変換
        sub_vals = [(label, subset_edge_values(G, V.load_importance(path, metric)))
                    for label, path in subsets]
        # 共通カラースケール (このスコアの全サブセット共通の vmax) で比較可能に
        vmax = max((v.max() for _, v in sub_vals if v.size), default=0.0)
        vmax = float(vmax) if vmax > 0 else 1.0
        V.log(f"スコア |{metric}|: vmax={vmax:.3f}, 強調 top-{args.top_n}")

        # 個別図 (共通レイアウト・群ごとの色)
        for label, vals in sub_vals:
            cmap = V.make_cmap(cmaps[label])
            fig, ax = plt.subplots(figsize=(15, 13))
            draw_subset(ax, G, pos, vals, targets, args.top_n, vmax, cmap,
                        f"{label}: top {args.top_n} important edges (|{metric}|)")
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax))
            sm.set_array([])
            fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.01, label=f"|{metric}|")
            out = os.path.join(args.fig_dir, f"subset_{label}_{metric}.png")
            fig.tight_layout()
            fig.savefig(out, dpi=200, bbox_inches="tight")
            plt.close(fig)
            V.log(f"出力: {out}")

        # 並べた比較図 (共通レイアウト・共通スケール・群ごとの色)
        n = len(sub_vals)
        ncols = 2 if n > 1 else 1
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(13 * ncols, 11 * nrows))
        axes = np.atleast_1d(axes).ravel()
        for ax, (label, vals) in zip(axes, sub_vals):
            cmap = V.make_cmap(cmaps[label])
            draw_subset(ax, G, pos, vals, targets, args.top_n, vmax, cmap, label)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax))
            sm.set_array([])
            fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.01, label=f"|{metric}|")
        for ax in axes[n:]:
            ax.set_axis_off()
        fig.suptitle(f"Per-subset important edges |{metric}| "
                     f"(top {args.top_n}, shared layout; hue = group)", fontsize=15)
        grid = os.path.join(args.fig_dir, f"subsets_grid_{metric}.png")
        fig.savefig(grid, dpi=170, bbox_inches="tight")
        plt.close(fig)
        V.log(f"出力: {grid}")

        # 統合図に載せる群 (全サンプル参照パネルは除く)
        group_vals = [(l, v) for l, v in sub_vals if l != "ALL_samples"]

        # 統合オーバーレイ図 (色の種類=群, 色の濃さ=重要度)
        fig, ax = plt.subplots(figsize=(17, 15))
        draw_overlay(ax, G, pos, group_vals, cmaps, targets, args.top_n, vmax,
                     f"Integrated per-group important edges |{metric}| "
                     f"(top {args.top_n}; hue = group, shade = importance)")
        overlay = os.path.join(args.fig_dir, f"subsets_overlay_{metric}.png")
        fig.tight_layout()
        fig.savefig(overlay, dpi=200, bbox_inches="tight")
        plt.close(fig)
        V.log(f"出力: {overlay}")

        # 3 チャンネル統合図 (色の種類=群, 色の濃さ=重要度, 線の太さ=ブートストラップ確率)
        width_note = ("width = bootstrap prob" if prob_edge_vals is not None
                      else "width = importance (no --edge-prob)")
        fig, ax = plt.subplots(figsize=(17, 15))
        draw_overlay(ax, G, pos, group_vals, cmaps, targets, args.top_n, vmax,
                     f"Integrated |{metric}| (top {args.top_n}; hue = group, "
                     f"shade = importance, {width_note})",
                     prob_vals=prob_edge_vals)
        multi = os.path.join(args.fig_dir, f"subsets_multichannel_{metric}.png")
        fig.tight_layout()
        fig.savefig(multi, dpi=200, bbox_inches="tight")
        plt.close(fig)
        V.log(f"出力: {multi}")
    V.log("完了")


if __name__ == "__main__":
    main()

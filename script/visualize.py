#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualize.py
============
構造学習結果 (out/ 以下) を読み込み、networkx + matplotlib で
ベイジアンネットワークを描画する汎用スクリプト。表示方法を切り替えた
複数バージョンの図を figures/ に出力する。

既定のパスは**カレントディレクトリ基準** (./out, ./figures, ./data/var_map.tsv,
./target_genes.txt) なので、任意の解析ディレクトリで実行できる。
ラッパスクリプト viz.sh / viz_bs.sh / viz_subsets.sh / viz_bs_subsets.sh も参照。

読み込むファイル (--out-dir, 既定 ./out):
  edges.tsv         : エッジ (ノード=列インデックス: u v)         [必須]
  edges_named.tsv   : 同じエッジを遺伝子名で表記 (行は edges.tsv と対応) [必須]
  edge_importance.tsv : エッジ重要度 (u v ΔlogL ΔBIC ΔK2 ΔBDeu ...) [任意]

ノードのインデックス→遺伝子名の対応は edges.tsv と edges_named.tsv が
行単位で 1:1 対応していることを利用して復元する (var_map.tsv の世代ずれに
依存しない)。

出力する図 (figures/):
  01_structure_full.png   : 全構造 (ノードサイズ=次数, ラベルはハブ/注目遺伝子のみ)
  02_importance_full.png  : 全構造 (エッジの色/太さ=重要度)
  03_importance_top.png   : 重要度上位 N エッジのみの部分グラフ (全ノードにラベル)
  04_targets_highlight.png: 注目遺伝子 (ホワイトリスト) を強調、関連エッジを着色
  05_target_ego.png       : 注目遺伝子とその近傍のみの部分グラフ
さらに edge_importance_named.tsv (重要度順の名前付きエッジ表) を出力する。
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")  # 画面表示なしで PNG 出力
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
import numpy as np


# 重要度メトリクスの列 (edge_importance.tsv のデータ行は 8 列)
IMP_COLUMNS = ["u", "v", "dlogL", "dBIC", "dK2", "dBDeu",
               "mean_dlogL_per_sample", "std_dlogL_per_sample"]

# 配色: 濃い色を避け、半透明にしてラベル (文字) を見やすくする
COL_NODE = "#add0ec"      # 通常ノード (淡い青)
COL_TARGET = "#f4a6a6"    # 注目遺伝子 (淡い赤)
COL_NEIGHBOR = "#b4e0bd"  # 近傍 (淡い緑)
COL_BG = "#e9e9e9"        # 背景ノード (淡いグレー)
NODE_ALPHA = 0.6          # ノードの不透明度 (半透明)
EDGE_ALPHA = 0.3          # 単色エッジの不透明度 (半透明)
EDGE_ALPHA_C = 0.5        # カラーマップエッジの不透明度 (半透明だが色は判別可)
# ラベルの背景 (白の半透明) を敷いて文字を読みやすくする
LABEL_BBOX = dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.55)

# 既定のエッジ配色: 重要度・ブートストラップ確率にかかわらず単色グラデーション
# (Blues) を用いる。--cmap で任意の matplotlib カラーマップ名に変更可能。
DEFAULT_CMAP = "Blues"

# 描画スケール: 既定値は数百ノードの大きな網に合わせてある。ノード数が少ない
# 網ではノード・エッジ・文字が小さすぎて見えなくなるため、set_style() で
# ノード数に応じて倍率を決め、各描画関数がこれを掛ける。
REF_NODES = 300.0        # この規模のとき倍率 1.0 (= 従来の見た目)
STYLE = {"node": 1.0, "edge": 1.0, "font": 7.0, "label": 6.0,
         "edge_alpha": EDGE_ALPHA, "edge_alpha_c": EDGE_ALPHA_C, "arrow": 1.0}


def set_style(n_nodes, override=None):
    """ノード数に応じた描画倍率を決める (override があればそれを使う)。"""
    k = override if override else float(np.clip(np.sqrt(REF_NODES / max(1, n_nodes)),
                                                1.0, 5.0))
    STYLE["node"] = k ** 2          # ノードは面積指定なので 2 乗
    STYLE["edge"] = k
    STYLE["arrow"] = k
    STYLE["font"] = float(np.clip(7.0 * np.sqrt(k), 7.0, 13.0))
    STYLE["label"] = float(np.clip(6.0 * np.sqrt(k), 6.0, 12.0))
    STYLE["edge_alpha"] = float(np.clip(EDGE_ALPHA * np.sqrt(k), EDGE_ALPHA, 0.75))
    STYLE["edge_alpha_c"] = float(np.clip(EDGE_ALPHA_C * np.sqrt(k),
                                          EDGE_ALPHA_C, 0.9))
    log(f"描画スケール: ノード数 {n_nodes} -> 倍率 {k:.2f} "
        f"(node x{STYLE['node']:.1f}, font {STYLE['font']:.1f}pt)")
    return k


def make_cmap(name=DEFAULT_CMAP, lo=0.20, hi=1.0):
    """単色系カラーマップを、最小値 (=0) でも白飛びせず視認できるよう
    下端を切り詰めて返す。未知の名前は Blues にフォールバック。"""
    try:
        base = plt.get_cmap(name)
    except ValueError:
        log(f"未知のカラーマップ '{name}' -> Blues を使用")
        base = plt.get_cmap(DEFAULT_CMAP)
    return LinearSegmentedColormap.from_list(
        f"{name}_trunc", base(np.linspace(lo, hi, 256)))


def log(msg):
    print(f"[viz] {msg}", file=sys.stderr)


def load_node_names(edges_path, named_path):
    """edges.tsv と edges_named.tsv の行対応から idx->遺伝子名 を復元。"""
    idx2name = {}
    with open(edges_path) as fe, open(named_path) as fn:
        for le, ln in zip(fe, fn):
            ea = le.strip().split("\t")
            na = ln.strip().split("\t")
            if len(ea) < 2 or len(na) < 2:
                continue
            u, v = int(ea[0]), int(ea[1])
            idx2name[u] = na[0]
            idx2name[v] = na[1]
    return idx2name


def load_edges(edges_path):
    edges = []
    with open(edges_path) as fp:
        for line in fp:
            a = line.strip().split("\t")
            if len(a) >= 2:
                edges.append((int(a[0]), int(a[1])))
    return edges


def load_importance(path, metric):
    """edge_importance.tsv を読み、(u,v)->重要度 の dict を返す。

    重要度は選択メトリクスの絶対値 (|ΔlogL| など)。値が大きいほど重要。
    """
    imp = {}
    col = IMP_COLUMNS.index(metric)
    with open(path) as fp:
        first = True
        for line in fp:
            a = line.strip().split("\t")
            if first:
                first = False
                # ヘッダ行 (u v ...) はスキップ
                if a and a[0] == "u":
                    continue
            if len(a) < col + 1:
                continue
            try:
                u, v = int(a[0]), int(a[1])
                val = abs(float(a[col]))
            except ValueError:
                continue
            imp[(u, v)] = val
    return imp


def load_edge_prob(path):
    """ブートストラップ確率ファイル (u v count prob, ヘッダ無し) を読み込み、
    (u,v)->prob の dict を返す。compute_bs_prob.py の --out 出力 (integ_edges_score.tsv)。
    """
    prob = {}
    with open(path) as fp:
        for line in fp:
            a = line.rstrip("\n").split("\t")
            if len(a) < 4:
                continue
            try:
                u, v = int(a[0]), int(a[1])
                prob[(u, v)] = float(a[3])
            except ValueError:
                continue  # ヘッダ行等はスキップ
    return prob


def load_targets(target_file, var_map):
    """注目遺伝子 (ホワイトリスト) 名の集合を返す。

    target_genes.txt があればそれを優先。無ければ var_map.tsv の
    whitelisted==1 列を使用する。いずれも無ければ空集合。
    """
    gene_map = var_map
    targets = set()
    if target_file and os.path.exists(target_file):
        with open(target_file, encoding="utf-8") as fp:
            for line in fp:
                t = line.strip()
                if t and not t.startswith("#"):
                    targets.add(t)
        log(f"注目遺伝子を {target_file} から読込 ({len(targets)} 件)")
    elif gene_map and os.path.exists(gene_map):
        with open(gene_map, encoding="utf-8") as fp:
            header = fp.readline().rstrip("\n").split("\t")
            try:
                ni = header.index("gene_name")
                wi = header.index("whitelisted")
            except ValueError:
                ni = wi = None
            if ni is not None:
                for line in fp:
                    a = line.rstrip("\n").split("\t")
                    if len(a) > wi and a[wi] == "1":
                        targets.add(a[ni])
        log(f"注目遺伝子を {gene_map} (whitelisted) から読込 ({len(targets)} 件)")
    return targets


def base_name(name):
    """重複回避で付与された "__ENSG..." 接尾辞を除いた素の遺伝子名。"""
    return name.split("__")[0]


def is_target(name, targets):
    return base_name(name) in targets or name in targets


def build_graph(edges, idx2name, imp):
    G = nx.DiGraph()
    for u, v in edges:
        nu = idx2name.get(u, str(u))
        nv = idx2name.get(v, str(v))
        w = imp.get((u, v), 0.0) if imp else 0.0
        G.add_edge(nu, nv, importance=w, _uv=(u, v))
    return G


def set_importance(G, imp):
    """各エッジの importance 属性を、指定スコアの重要度 dict で更新する。"""
    for u, v in G.edges():
        G[u][v]["importance"] = imp.get(G[u][v].get("_uv", (None, None)), 0.0)


# ---------------------------------------------------------------------------
# 描画ユーティリティ
# ---------------------------------------------------------------------------
def _node_sizes(G, nodes, scale=18, base=20):
    m = STYLE["node"]
    return [m * (base + scale * G.degree(n)) for n in nodes]


def _w(width):
    """線幅 (スカラー or 配列) に描画倍率を掛ける。"""
    return width * STYLE["edge"]


def _a(arrowsize):
    return max(3.0, arrowsize * STYLE["arrow"])


def save_fig(path, title):
    plt.title(title, fontsize=13)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    log(f"出力: {path}")


def draw_structure_full(G, pos, targets, fig_dir, hub_labels=15):
    plt.figure(figsize=(16, 16))
    nodes = list(G.nodes())
    tgt = [n for n in nodes if is_target(n, targets)]
    oth = [n for n in nodes if not is_target(n, targets)]
    nx.draw_networkx_edges(G, pos, edge_color="#9aa3ad", width=_w(0.5),
                           arrowsize=_a(6), alpha=STYLE["edge_alpha"])
    nx.draw_networkx_nodes(G, pos, nodelist=oth, node_size=_node_sizes(G, oth),
                           node_color=COL_NODE, linewidths=0, alpha=NODE_ALPHA)
    if tgt:
        nx.draw_networkx_nodes(G, pos, nodelist=tgt,
                               node_size=_node_sizes(G, tgt, scale=30, base=120),
                               node_color=COL_TARGET, linewidths=0, alpha=NODE_ALPHA)
    # ラベルは高次数ハブ + 注目遺伝子のみ (混雑回避)
    hubs = sorted(G.degree, key=lambda x: x[1], reverse=True)[:max(0, hub_labels)]
    label_nodes = {n for n, _ in hubs} | set(tgt)
    nx.draw_networkx_labels(G, pos, labels={n: n for n in label_nodes},
                            font_size=STYLE["font"], bbox=LABEL_BBOX)
    save_fig(os.path.join(fig_dir, "01_structure_full.png"),
             f"BN structure (nodes={G.number_of_nodes()}, "
             f"edges={G.number_of_edges()}); red=target, label=hubs+targets")


def draw_importance_full(G, pos, fig_dir, metric, cmap, suffix=""):
    plt.figure(figsize=(16, 16))
    edges = list(G.edges())
    weights = np.array([G[u][v]["importance"] for u, v in edges])
    if weights.max() <= 0:
        log(f"重要度が全て 0/未取得のため 02_importance_full{suffix} はスキップ")
        plt.close()
        return
    wn = weights / weights.max()
    nx.draw_networkx_nodes(G, pos, node_size=_node_sizes(G, list(G.nodes())),
                           node_color=COL_NODE, linewidths=0, alpha=NODE_ALPHA)
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=weights,
                           edge_cmap=cmap, width=_w(0.5 + 3.5 * wn),
                           arrowsize=_a(7), alpha=STYLE["edge_alpha_c"],
                           edge_vmin=0, edge_vmax=weights.max())
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=0, vmax=weights.max()))
    sm.set_array([])
    plt.colorbar(sm, ax=plt.gca(), fraction=0.03, pad=0.01,
                 label=f"edge importance |{metric}|")
    save_fig(os.path.join(fig_dir, f"02_importance_full{suffix}.png"),
             f"Edge importance (|{metric}|): color/width = importance")


def draw_importance_top(G, targets, fig_dir, metric, top_n, seed, cmap, suffix=""):
    edges = [(u, v, d["importance"]) for u, v, d in G.edges(data=True)]
    edges = [e for e in edges if e[2] > 0]
    if not edges:
        log(f"重要度付きエッジが無いため 03_importance_top{suffix} はスキップ")
        return
    edges.sort(key=lambda x: x[2], reverse=True)
    top = edges[:top_n]
    H = nx.DiGraph()
    for u, v, w in top:
        H.add_edge(u, v, importance=w)
    pos = nx.spring_layout(H, seed=seed, k=0.6)
    plt.figure(figsize=(14, 12))
    weights = np.array([w for _, _, w in top])
    wn = weights / weights.max()
    tgt = [n for n in H.nodes() if is_target(n, targets)]
    oth = [n for n in H.nodes() if not is_target(n, targets)]
    nx.draw_networkx_nodes(H, pos, nodelist=oth, node_size=300,
                           node_color=COL_NODE, linewidths=0, alpha=NODE_ALPHA)
    if tgt:
        nx.draw_networkx_nodes(H, pos, nodelist=tgt, node_size=450,
                               node_color=COL_TARGET, linewidths=0, alpha=NODE_ALPHA)
    nx.draw_networkx_edges(H, pos, edge_color=weights, edge_cmap=cmap,
                           width=1.0 + 4.0 * wn, arrowsize=12,
                           alpha=STYLE["edge_alpha_c"],
                           edge_vmin=0, edge_vmax=weights.max())
    nx.draw_networkx_labels(H, pos, font_size=8, bbox=LABEL_BBOX)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=0, vmax=weights.max()))
    sm.set_array([])
    plt.colorbar(sm, ax=plt.gca(), fraction=0.03, pad=0.01,
                 label=f"|{metric}|")
    save_fig(os.path.join(fig_dir, f"03_importance_top{suffix}.png"),
             f"Top {len(top)} edges by importance |{metric}|")


def draw_targets_highlight(G, pos, targets, fig_dir):
    if not targets:
        log("注目遺伝子が無いため 04_targets_highlight はスキップ")
        return
    tgt_nodes = [n for n in G.nodes() if is_target(n, targets)]
    if not tgt_nodes:
        log("注目遺伝子がグラフ中に存在しないため 04_targets_highlight はスキップ")
        return
    tgt_set = set(tgt_nodes)
    inc_edges = [(u, v) for u, v in G.edges() if u in tgt_set or v in tgt_set]
    plt.figure(figsize=(16, 16))
    # 背景: 全グラフを淡色で
    nx.draw_networkx_edges(G, pos, edge_color="#dddddd", width=_w(0.4),
                           arrowsize=_a(4), alpha=STYLE["edge_alpha"])
    nx.draw_networkx_nodes(G, pos, node_size=15 * STYLE["node"], node_color=COL_BG,
                           linewidths=0, alpha=NODE_ALPHA)
    # 強調: 注目遺伝子に接続するエッジと、注目遺伝子の近傍
    nbrs = set()
    for u, v in inc_edges:
        nbrs.add(u)
        nbrs.add(v)
    nbr_only = [n for n in nbrs if n not in tgt_set]
    nx.draw_networkx_edges(G, pos, edgelist=inc_edges, edge_color="#f0a64d",
                           width=_w(1.5), arrowsize=_a(10), alpha=0.55)
    nx.draw_networkx_nodes(G, pos, nodelist=nbr_only, node_size=120 * STYLE["node"],
                           node_color=COL_NEIGHBOR, linewidths=0, alpha=NODE_ALPHA)
    nx.draw_networkx_nodes(G, pos, nodelist=tgt_nodes, node_size=300 * STYLE["node"],
                           node_color=COL_TARGET, linewidths=0, alpha=NODE_ALPHA)
    label_nodes = {n: n for n in (set(tgt_nodes) | set(nbr_only))}
    nx.draw_networkx_labels(G, pos, labels=label_nodes, font_size=STYLE["font"],
                            bbox=LABEL_BBOX)
    save_fig(os.path.join(fig_dir, "04_targets_highlight.png"),
             f"Target genes (red) + neighbors (green); "
             f"{len(tgt_nodes)} targets, {len(inc_edges)} incident edges")


def draw_target_ego(G, targets, fig_dir, metric, seed, cmap, suffix=""):
    if not targets:
        log("注目遺伝子が無いため 05_target_ego はスキップ")
        return
    tgt_nodes = [n for n in G.nodes() if is_target(n, targets)]
    if not tgt_nodes:
        log("注目遺伝子がグラフ中に存在しないため 05_target_ego はスキップ")
        return
    keep = set(tgt_nodes)
    for n in tgt_nodes:
        keep |= set(G.successors(n)) | set(G.predecessors(n))
    H = G.subgraph(keep).copy()
    pos = nx.spring_layout(H, seed=seed, k=0.5)
    plt.figure(figsize=(15, 13))
    weights = np.array([H[u][v]["importance"] for u, v in H.edges()])
    has_imp = weights.size > 0 and weights.max() > 0
    if has_imp:
        wn = weights / weights.max()
        nx.draw_networkx_edges(H, pos, edge_color=weights,
                               edge_cmap=cmap,
                               width=0.8 + 3.0 * wn, arrowsize=10, alpha=EDGE_ALPHA_C,
                               edge_vmin=0, edge_vmax=weights.max())
        sm = plt.cm.ScalarMappable(cmap=cmap,
                                   norm=plt.Normalize(0, weights.max()))
        sm.set_array([])
        plt.colorbar(sm, ax=plt.gca(), fraction=0.03, pad=0.01,
                     label=f"|{metric}|")
    else:
        nx.draw_networkx_edges(H, pos, edge_color="#9aa3ad", width=0.8,
                               arrowsize=10, alpha=EDGE_ALPHA)
    tgt = [n for n in H.nodes() if is_target(n, targets)]
    oth = [n for n in H.nodes() if not is_target(n, targets)]
    nx.draw_networkx_nodes(H, pos, nodelist=oth, node_size=200,
                           node_color=COL_NODE, linewidths=0, alpha=NODE_ALPHA)
    nx.draw_networkx_nodes(H, pos, nodelist=tgt, node_size=400,
                           node_color=COL_TARGET, linewidths=0, alpha=NODE_ALPHA)
    nx.draw_networkx_labels(H, pos, font_size=7, bbox=LABEL_BBOX)
    save_fig(os.path.join(fig_dir, f"05_target_ego{suffix}.png"),
             f"Ego network of target genes (|{metric}|; "
             f"nodes={H.number_of_nodes()}, edges={H.number_of_edges()})")


def draw_bootstrap_prob(G, pos, prob, fig_dir, cmap):
    """エッジのブートストラップ確率 (出現頻度) で色と太さを表現した全体図。
    コンセンサス網の安定性を可視化する (確率が高いほど安定なエッジ)。
    重要度と同じ単色グラデーション (既定 Blues) を用いる。
    """
    edges = list(G.edges())
    vals = np.array([prob.get(G[u][v].get("_uv", (None, None)), 0.0)
                     for u, v in edges])
    if vals.max() <= 0:
        log("ブートストラップ確率が取得できないため 06_bootstrap_prob はスキップ")
        return
    plt.figure(figsize=(16, 16))
    nx.draw_networkx_nodes(G, pos, node_size=_node_sizes(G, list(G.nodes())),
                           node_color=COL_NODE, linewidths=0, alpha=NODE_ALPHA)
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=vals,
                           edge_cmap=cmap, width=_w(0.5 + 3.5 * vals),
                           arrowsize=_a(7), alpha=STYLE["edge_alpha_c"],
                           edge_vmin=0.0, edge_vmax=1.0)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.0, vmax=1.0))
    sm.set_array([])
    plt.colorbar(sm, ax=plt.gca(), fraction=0.03, pad=0.01,
                 label="bootstrap probability")
    save_fig(os.path.join(fig_dir, "06_bootstrap_prob.png"),
             f"Bootstrap edge probability (consensus stability); "
             f"mean={vals.mean():.2f}, edges={len(edges)}")


def write_named_importance(G, idx2name, imp, path, metric):
    """重要度順の名前付きエッジ表を出力 (読みやすさ用)。"""
    if not imp:
        return
    rows = []
    for (u, v), w in imp.items():
        rows.append((w, idx2name.get(u, str(u)), idx2name.get(v, str(v)), u, v))
    rows.sort(reverse=True)
    with open(path, "w", encoding="utf-8") as fp:
        fp.write(f"rank\tu_name\tv_name\tu_idx\tv_idx\timportance_abs_{metric}\n")
        for i, (w, un, vn, u, v) in enumerate(rows, 1):
            fp.write(f"{i}\t{un}\t{vn}\t{u}\t{v}\t{w:.6g}\n")
    log(f"出力: {path} (重要度順エッジ表)")


def main():
    ap = argparse.ArgumentParser(description="BN 結果を networkx で可視化")
    ap.add_argument("--out-dir", default="./out",
                    help="学習結果ディレクトリ (edges.tsv 等; 既定 ./out)")
    ap.add_argument("--var-map", dest="var_map", default="./data/var_map.tsv",
                    help="変数対応表 var_map.tsv (注目遺伝子フォールバック用)")
    ap.add_argument("--target-file", default="./target_genes.txt",
                    help="ホワイトリスト (注目遺伝子) ファイル")
    ap.add_argument("--fig-dir", default="./figures",
                    help="図の出力先 (既定 ./figures)")
    # 入力ファイル名 (--out-dir からの相対 or 絶対パス)。コンセンサス網を描く場合は
    # integ_edges2.tsv / integ_edges_named.tsv / integ_edge_importance.tsv を指定。
    ap.add_argument("--edges", default="edges.tsv",
                    help="インデックス表記のエッジファイル (既定 edges.tsv)")
    ap.add_argument("--edges-named", default="edges_named.tsv",
                    help="遺伝子名表記のエッジファイル (--edges と行対応; 既定 edges_named.tsv)")
    ap.add_argument("--importance", default="edge_importance.tsv",
                    help="エッジ重要度ファイル (既定 edge_importance.tsv)")
    ap.add_argument("--edge-prob", default=None,
                    help="ブートストラップ確率ファイル (u v count prob; 例 integ_edges_score.tsv)。"
                         "指定すると確率で着色した図も出力する")
    ap.add_argument("--metrics", default="dlogL,dBIC,dK2,dBDeu",
                    help="重要度可視化を行うスコア (カンマ区切り)。"
                         f"選択肢: {','.join(IMP_COLUMNS[2:])} (既定 dlogL,dBIC,dK2,dBDeu)")
    ap.add_argument("--top-n", type=int, default=60,
                    help="重要度上位エッジ数 (03 図, 既定 60)")
    ap.add_argument("--hub-labels", type=int, default=15,
                    help="全体図 (01) でラベルを付けるハブ (高次数ノード) の数")
    ap.add_argument("--style-scale", type=float, default=None,
                    help="ノード・線・文字の大きさの倍率。既定はノード数から自動決定 "
                         f"(ノード数 {int(REF_NODES)} で 1.0、少ないほど大きく描く)")
    ap.add_argument("--cmap", default=DEFAULT_CMAP,
                    help="エッジの色に用いる matplotlib カラーマップ名。重要度・"
                         f"ブートストラップ確率とも共通 (既定 {DEFAULT_CMAP})。"
                         "例: Blues / Greens / Reds / viridis / plasma")
    ap.add_argument("--layout", default="spring",
                    choices=["spring", "kamada", "circular"],
                    help="全体図のレイアウト (既定 spring)")
    ap.add_argument("--seed", type=int, default=42, help="レイアウト乱数シード")
    args = ap.parse_args()

    def resolve(p):
        return p if os.path.isabs(p) else os.path.join(args.out_dir, p)

    edges_path = resolve(args.edges)
    named_path = resolve(args.edges_named)
    imp_path = resolve(args.importance)
    prob_path = resolve(args.edge_prob) if args.edge_prob else None

    for p in (edges_path, named_path):
        if not os.path.exists(p):
            sys.exit(f"[viz] 必須ファイルがありません: {p}")

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    bad = [m for m in metrics if m not in IMP_COLUMNS[2:]]
    if bad:
        sys.exit(f"[viz] 未知のスコア: {bad} (選択肢: {IMP_COLUMNS[2:]})")

    idx2name = load_node_names(edges_path, named_path)
    edges = load_edges(edges_path)
    if not edges:
        sys.exit(f"[viz] エッジが 0 件です: {edges_path}")
    has_imp = os.path.exists(imp_path)
    prob = load_edge_prob(prob_path) if prob_path and os.path.exists(prob_path) else {}
    log(f"ノード {len(idx2name)}, エッジ {len(edges)}, "
        f"ブートストラップ確率 {len(prob)} 件, スコア={metrics if has_imp else 'なし'}")
    targets = load_targets(args.target_file, args.var_map)

    G = build_graph(edges, idx2name, {})
    os.makedirs(args.fig_dir, exist_ok=True)
    cmap = make_cmap(args.cmap)
    set_style(G.number_of_nodes(), args.style_scale)

    # 全体図共通レイアウト (スコアに依らず 1 度だけ計算)
    log(f"レイアウト計算中 ({args.layout}) ...")
    if args.layout == "spring":
        pos = nx.spring_layout(G, seed=args.seed, k=1.0 / np.sqrt(max(1, G.number_of_nodes())))
    elif args.layout == "kamada":
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.circular_layout(G)

    # スコアに依らない図は 1 度だけ
    draw_structure_full(G, pos, targets, args.fig_dir, args.hub_labels)
    draw_targets_highlight(G, pos, targets, args.fig_dir)
    if prob:
        draw_bootstrap_prob(G, pos, prob, args.fig_dir, cmap)

    # 重要度ネットワーク図はスコアごとに出力 (ファイル名に _<metric> を付与)
    if has_imp:
        for metric in metrics:
            imp = load_importance(imp_path, metric)
            set_importance(G, imp)
            suffix = f"_{metric}"
            log(f"スコア |{metric}|: {len(imp)} エッジ -> 図を出力")
            draw_importance_full(G, pos, args.fig_dir, metric, cmap, suffix)
            draw_importance_top(G, targets, args.fig_dir, metric, args.top_n,
                                args.seed, cmap, suffix)
            draw_target_ego(G, targets, args.fig_dir, metric, args.seed, cmap, suffix)
            write_named_importance(
                G, idx2name, imp,
                os.path.join(args.fig_dir, f"edge_importance_named{suffix}.tsv"),
                metric)
    else:
        log("重要度ファイルが無いため重要度ネットワーク図はスキップ")

    log("完了")


if __name__ == "__main__":
    main()

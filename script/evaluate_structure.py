#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
evaluate_structure.py
=====================
正解ネットワーク (真の DAG) が分かっている場合に、学習された構造の精度を
以下の指標で評価する。

  1. Structural Hamming Distance (SHD)
  2. Directed Edge  Precision / Recall / F1   (向きまで一致するか)
  3. Skeleton Edge  Precision / Recall / F1   (向きを無視した骨格の一致)
  4. Structural Intervention Distance (SID)
  5. KL divergence  KL(P_true || P_learned)

--- 各指標の定義 -------------------------------------------------------------

SHD
  「学習した DAG を正解 DAG に変換するのに必要な編集回数」。本スクリプトでは
  よく使われる次の数え方を採用する (bnlearn の shd() と同じ):
      SHD = (欠損エッジ数) + (余分なエッジ数) + (向きだけ違うエッジ数)
  骨格が違えば 1、骨格は合っていて向きだけ違えば 1 と数える (反転を 2 とは
  数えない)。小さいほど良い。

Directed / Skeleton Precision, Recall, F1
  Directed は順序つきペア (u -> v)、Skeleton は無向ペア {u, v} を集合とみなして
      Precision = TP / (学習エッジ数),  Recall = TP / (正解エッジ数)
  を計算する。Skeleton は向きの誤りに寛容なので、Directed より必ず高くなる。
  観測データだけではマルコフ同値な DAG を区別できないため、向きの一致 (Directed)
  は原理的に低く出やすい。

SID (Peters & Bühlmann, 2015; "Structural Intervention Distance for Evaluating
Causal Graphs", Neural Computation 27(3))
  「学習した DAG を信じて親調整 (parent adjustment) で介入分布を計算したとき、
  何組の (i, j) で誤った答えを得るか」を数える因果的な距離。
  順序つきペア (i, j), i != j それぞれについて、学習 DAG H での i の親集合
  S = pa_H(i) が、真の DAG G において (i, j) の妥当な調整集合 (バックドア基準を
  満たす集合) であるかを判定し、満たさないペアを 1 つの誤りとして数える。
      妥当 <=> j not in S
               かつ S が G における i の子孫を含まない
               かつ G から i の出て行く辺を除いたグラフで S が i と j を d 分離する
  値域は 0 〜 p(p-1)。SID(G, G) = 0。SHD が小さくても SID が大きい (= 因果効果の
  推定を誤る) ことがあるため、構造の「使い道」に近い評価になる。

KL divergence
  KL(P_true || P_learned) を全状態空間の厳密な列挙で計算する (単位: nat)。
  P_learned は「学習した構造 + 学習に使ったデータから推定した CPT」の同時分布で、
  CPT は Dirichlet 平滑化 (--alpha) つきで推定する:
      P(v = k | pa = j) = (n_jk + alpha) / (n_j + alpha * r_v)
  alpha > 0 なら P_learned はどの状態にも正の確率を与えるので KL は有限になる。
  変数の基数は正解ネットワーク (--true-bif) のものを使うため、データに出現
  しなかった状態も正しく扱われる。状態空間が --max-states を超える場合は
  計算を省略する。

--- 使い方 -------------------------------------------------------------------

    # 正解が BIF で与えられる場合 (5 指標すべて)
    python3 evaluate_structure.py \
        --true-bif asia.bif \
        --pred-edges out/edges.tsv \
        --input data/asia_n1000.tsv \
        --out eval.tsv

    # 正解がエッジ表だけの場合 (KL 以外)
    python3 evaluate_structure.py \
        --true-edges true_edges.tsv \
        --pred-edges out/edges.tsv --input data.tsv \
        --out eval.tsv

    # ベンチマークの 1 行として追記する
    python3 evaluate_structure.py ... --append benchmark.tsv \
        --extra network=asia --extra n=1000 --extra score=bic

エッジファイルは 1 行 1 エッジの "u<TAB>v" (u -> v)。整数なら列インデックス、
文字列なら変数名として解釈する (自動判定)。インデックス表記の場合はノード名を
--input のヘッダか --true-bif から復元する。
"""

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def log(*args):
    print("[eval]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[eval] エラー: {msg}")


# ---------------------------------------------------------------------------
# 入出力
# ---------------------------------------------------------------------------

def read_edge_file(path, node_names):
    """エッジ表を (u, v) 名前ペアの集合として読む。整数なら列位置とみなす。"""
    if not os.path.isfile(path):
        die(f"{path} がありません")
    edges = []
    with open(path, encoding="utf-8") as fp:
        for lineno, line in enumerate(fp, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", "\t").split()
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            if a.lstrip("-").isdigit() and b.lstrip("-").isdigit():
                if node_names is None:
                    die(f"{path} はインデックス表記ですが、ノード名が不明です "
                        "(--input か --true-bif を指定してください)")
                ia, ib = int(a), int(b)
                if not (0 <= ia < len(node_names) and 0 <= ib < len(node_names)):
                    die(f"{path}:{lineno} のノード番号が範囲外です: {a} {b} "
                        f"(変数は {len(node_names)} 個)")
                a, b = node_names[ia], node_names[ib]
            edges.append((a, b))
    return edges


def read_data(path, node_names):
    """離散化済み TSV を (ヘッダ, int 配列) として読む。"""
    if not os.path.isfile(path):
        die(f"{path} がありません")
    with open(path, encoding="utf-8") as fp:
        header = fp.readline().rstrip("\n").split("\t")
    X = np.loadtxt(path, dtype=np.int64, delimiter="\t", skiprows=1, ndmin=2)
    if X.shape[1] != len(header):
        die(f"{path}: ヘッダ {len(header)} 列に対しデータ {X.shape[1]} 列")
    if node_names is not None and header != list(node_names):
        die(f"{path} のヘッダが正解ネットワークの変数と一致しません。\n"
            f"  データ: {header[:5]}...\n  正解  : {list(node_names)[:5]}...")
    return header, X


# ---------------------------------------------------------------------------
# グラフのユーティリティ
# ---------------------------------------------------------------------------

def build_parents(nodes, edges):
    parents = {v: set() for v in nodes}
    for u, v in edges:
        if u not in parents or v not in parents:
            die(f"エッジ ({u}, {v}) に未知のノードが含まれます")
        if u == v:
            die(f"自己ループがあります: {u}")
        parents[v].add(u)
    return parents


def children_of(nodes, parents):
    ch = {v: set() for v in nodes}
    for v in nodes:
        for p in parents[v]:
            ch[p].add(v)
    return ch


def is_dag(nodes, parents):
    color = {v: 0 for v in nodes}

    def visit(v):
        color[v] = 1
        for p in parents[v]:
            if color[p] == 1:
                return False
            if color[p] == 0 and not visit(p):
                return False
        color[v] = 2
        return True

    return all(color[v] != 0 or visit(v) for v in nodes)


def descendants(nodes, parents, start):
    """start の真の子孫集合 (start 自身は含まない)。"""
    ch = children_of(nodes, parents)
    seen, stack = set(), [start]
    while stack:
        v = stack.pop()
        for c in ch[v]:
            if c not in seen:
                seen.add(c)
                stack.append(c)
    seen.discard(start)
    return seen


def ancestors_of_set(parents, targets):
    """targets 自身を含む祖先集合。"""
    seen, stack = set(targets), list(targets)
    while stack:
        v = stack.pop()
        for p in parents[v]:
            if p not in seen:
                seen.add(p)
                stack.append(p)
    return seen


def d_separated(nodes, parents, x, y, z):
    """DAG 上で x と y が z によって d 分離されるか (モラル化祖先グラフ法)。"""
    z = set(z)
    keep = ancestors_of_set(parents, {x, y} | z)
    adj = {v: set() for v in keep}
    for v in keep:
        pa = [p for p in parents[v] if p in keep]
        for p in pa:
            adj[v].add(p)
            adj[p].add(v)
        for a in range(len(pa)):
            for b in range(a + 1, len(pa)):
                adj[pa[a]].add(pa[b])
                adj[pa[b]].add(pa[a])
    if x in z or y in z:
        return True                      # 端点が条件集合にある場合は分離扱い
    seen, stack = {x}, [x]
    while stack:
        v = stack.pop()
        if v == y:
            return False
        for w in adj[v]:
            if w not in seen and w not in z:
                seen.add(w)
                stack.append(w)
    return y not in seen


# ---------------------------------------------------------------------------
# 指標
# ---------------------------------------------------------------------------

def prf(tp, n_pred, n_true):
    precision = tp / n_pred if n_pred else (1.0 if n_true == 0 else 0.0)
    recall = tp / n_true if n_true else (1.0 if n_pred == 0 else 0.0)
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return precision, recall, f1


def structural_metrics(nodes, true_edges, pred_edges):
    """SHD と directed / skeleton の precision, recall, F1。"""
    T = set(true_edges)
    P = set(pred_edges)
    skelT = {frozenset(e) for e in T}
    skelP = {frozenset(e) for e in P}

    tp_dir = len(T & P)
    common_skel = skelT & skelP
    n_reversed = len({frozenset(e) for e in P if e not in T and (e[1], e[0]) in T})
    missing = len(skelT - skelP)
    extra = len(skelP - skelT)

    p_dir, r_dir, f_dir = prf(tp_dir, len(P), len(T))
    p_skel, r_skel, f_skel = prf(len(common_skel), len(skelP), len(skelT))

    return {
        "n_true_edges": len(T),
        "n_pred_edges": len(P),
        "true_positive_directed": tp_dir,
        "false_positive_directed": len(P) - tp_dir,
        "false_negative_directed": len(T) - tp_dir,
        "reversed_edges": n_reversed,
        "missing_edges": missing,
        "extra_edges": extra,
        "shd": missing + extra + n_reversed,
        "precision_directed": p_dir,
        "recall_directed": r_dir,
        "f1_directed": f_dir,
        "precision_skeleton": p_skel,
        "recall_skeleton": r_skel,
        "f1_skeleton": f_skel,
    }


def adjustment_is_valid(nodes, true_parents, desc, ch_true, i, j, S):
    """S = pa_H(i) が真の DAG で (i, j) の妥当な調整集合かを判定する。

    Perkovic らの (generalized) adjustment criterion を用いる。これは
    Peters & Bühlmann の SID 実装 (R パッケージ SID) が用いる判定と等価:

      * i から j への因果的経路 (proper causal path) 上の子 c を集める
      * 禁止集合 = それらの c とその子孫。S がこれと交われば無効
      * それらの c について辺 i -> c を除いたグラフ (proper back-door graph)
        で S が i と j を d 分離しなければ無効 (= 非因果的経路が開いている)
    """
    if j in desc[i]:
        causal_children = [c for c in ch_true[i] if c == j or j in desc[c]]
        forbidden = set()
        for c in causal_children:
            forbidden.add(c)
            forbidden |= desc[c]
        if S & forbidden:
            return False
        pbd_parents = {v: set(true_parents[v]) for v in nodes}
        for c in causal_children:
            pbd_parents[c].discard(i)
    else:
        pbd_parents = true_parents
    return d_separated(nodes, pbd_parents, i, j, S)


def structural_intervention_distance(nodes, true_parents, pred_parents):
    """SID: 親調整で介入分布を誤って計算する順序つきペア (i, j) の数。"""
    ch_true = children_of(nodes, true_parents)
    desc = {v: descendants(nodes, true_parents, v) for v in nodes}
    mistakes = 0
    for i in nodes:
        S = set(pred_parents[i])
        if S == set(true_parents[i]):
            continue                     # 親集合が真と一致すれば必ず正しい
        for j in nodes:
            if j == i:
                continue
            if not adjustment_is_valid(nodes, true_parents, desc, ch_true, i, j, S):
                mistakes += 1
    return mistakes


def fit_cpts(names, cards, parents, X, alpha):
    """学習構造 + データから CPT を Dirichlet 平滑化つきで推定する。"""
    idx = {v: i for i, v in enumerate(names)}
    cpt = {}
    n = X.shape[0]
    for v in names:
        pa = sorted(parents[v], key=lambda p: idx[p])
        k = cards[idx[v]]
        shape = tuple(cards[idx[p]] for p in pa) + (k,)
        counts = np.zeros(shape, dtype=np.float64)
        flat = counts.reshape(-1, k)
        if pa:
            ridx = np.zeros(n, dtype=np.int64)
            for p in pa:
                ridx = ridx * cards[idx[p]] + X[:, idx[p]]
        else:
            ridx = np.zeros(n, dtype=np.int64)
        np.add.at(flat, (ridx, X[:, idx[v]]), 1.0)
        counts = flat.reshape(shape)
        total = counts.sum(axis=-1, keepdims=True)
        denom = total + alpha * k
        with np.errstate(invalid="ignore", divide="ignore"):
            table = (counts + alpha) / denom
        # 観測ゼロかつ alpha=0 の親設定は一様分布で埋める
        bad = (denom <= 0)
        if np.any(bad):
            table = np.where(bad, 1.0 / k, table)
        cpt[v] = table
        parents[v] = pa            # 軸の順序を確定させる
    return cpt


def kl_divergence(true_bn, pred_parents, X, alpha):
    """KL(P_true || P_learned) を厳密列挙で計算する (nat)。"""
    from bif_io import BayesNet

    names = true_bn.names
    cards = true_bn.cards
    parents = {v: sorted(pred_parents[v], key=lambda p: true_bn.index[p])
               for v in names}
    cpt = fit_cpts(names, cards, parents, X, alpha)
    learned = BayesNet(names, true_bn.states, parents, cpt)

    P = true_bn.joint()
    Q = learned.joint()
    mask = P > 0
    if np.any(Q[mask] <= 0):
        return float("inf")
    return float(np.sum(P[mask] * (np.log(P[mask]) - np.log(Q[mask]))))


# ---------------------------------------------------------------------------
# メイン
# ---------------------------------------------------------------------------

METRIC_ORDER = [
    "n_nodes", "n_true_edges", "n_pred_edges",
    "true_positive_directed", "false_positive_directed",
    "false_negative_directed", "reversed_edges", "missing_edges", "extra_edges",
    "shd",
    "precision_directed", "recall_directed", "f1_directed",
    "precision_skeleton", "recall_skeleton", "f1_skeleton",
    "sid", "sid_max", "sid_normalized",
    "kl_divergence",
    "true_is_dag", "pred_is_dag",
    "n_true_edges_raw", "n_pred_edges_raw", "n_pred_not_evaluable",
]


def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isinf(value):
            return "inf"
        if math.isnan(value):
            return "NA"
        return f"{value:.6g}"
    return str(value)


def main():
    ap = argparse.ArgumentParser(
        description="正解 DAG に対する構造学習結果の評価 (SHD / P-R-F1 / SID / KL)")
    ap.add_argument("--true-bif", default=None,
                    help="正解ネットワーク (BIF)。KL の計算に必要")
    ap.add_argument("--true-edges", default=None,
                    help="正解エッジ表 (--true-bif の代わり; KL は計算しない)")
    ap.add_argument("--pred-edges", required=True,
                    help="学習エッジ表 (fast_bn の edges.tsv など)")
    ap.add_argument("--input", default=None,
                    help="学習に使ったデータ TSV。ノード名の復元と KL 用の "
                         "パラメータ推定に使う")
    ap.add_argument("--alpha", type=float, default=1.0,
                    help="KL 用 CPT 推定の Dirichlet 平滑化 (default=1.0)")
    ap.add_argument("--max-states", type=int, default=2_000_000,
                    help="KL を計算する状態空間の上限 (default=2000000)")
    ap.add_argument("--eval-pairs", default=None,
                    help="評価対象の (u, v) ペアを列挙したファイル。DREAM5 のように "
                         "gold standard が一部のペア (TF x 遺伝子) しか判定して "
                         "いない場合に、そのペアだけで指標を計算する")
    ap.add_argument("--skip-sid", action="store_true",
                    help="SID を計算しない (大規模ネットワークで時間がかかる場合)")
    ap.add_argument("--max-sid-nodes", type=int, default=500,
                    help="SID を計算するノード数の上限 (default=500)。"
                         "SID は O(p^2) 回の d 分離判定を行うため大規模だと重い")
    ap.add_argument("--out", default=None, help="指標を書き出す TSV (2 行)")
    ap.add_argument("--append", default=None,
                    help="指標を 1 行として追記する TSV (無ければヘッダも書く)")
    ap.add_argument("--json", dest="json_out", default=None, help="JSON 出力")
    ap.add_argument("--extra", action="append", default=[], metavar="KEY=VALUE",
                    help="出力に付け加える列 (複数指定可)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    if not args.true_bif and not args.true_edges:
        die("--true-bif か --true-edges のどちらかが必要です")

    # --- ノード名の確定 -----------------------------------------------------
    true_bn = None
    node_names = None
    if args.true_bif:
        from bif_io import read_bif
        true_bn = read_bif(args.true_bif)
        node_names = list(true_bn.names)
    X = None
    if args.input:
        header, X = read_data(args.input, node_names)
        if node_names is None:
            node_names = header

    # --- エッジの読み込み ---------------------------------------------------
    if args.true_bif:
        true_edges = true_bn.edges()
    else:
        true_edges = read_edge_file(args.true_edges, node_names)
        if node_names is None:
            node_names = sorted({v for e in true_edges for v in e})
    pred_edges = read_edge_file(args.pred_edges, node_names)

    nodes = list(node_names)
    known = set(nodes)
    unknown = {v for e in true_edges + pred_edges for v in e} - known
    if unknown:
        die(f"ノード集合に無い名前がエッジに含まれます: {sorted(unknown)[:5]}")

    true_parents = build_parents(nodes, true_edges)
    pred_parents = build_parents(nodes, pred_edges)

    # --- 評価対象ペアの制限 --------------------------------------------------
    # SID は「グラフ全体」に対する指標なので、フィルタ前のグラフで計算する。
    n_excluded = 0
    metric_true, metric_pred = true_edges, pred_edges
    if args.eval_pairs:
        evaluable = set(read_edge_file(args.eval_pairs, node_names))
        skel_evaluable = {frozenset(e) for e in evaluable}
        metric_true = [e for e in true_edges if e in evaluable]
        metric_pred = [e for e in pred_edges if frozenset(e) in skel_evaluable]
        n_excluded = len(pred_edges) - len(metric_pred)
        if not args.quiet:
            log(f"評価対象ペアで制限: 正解 {len(true_edges)} -> {len(metric_true)}, "
                f"学習 {len(pred_edges)} -> {len(metric_pred)} "
                f"(判定対象外の学習エッジ {n_excluded} 本を除外)")

    # --- 指標 ---------------------------------------------------------------
    results = {"n_nodes": len(nodes)}
    results.update(structural_metrics(nodes, metric_true, metric_pred))
    # n_true_edges / n_pred_edges は --eval-pairs で絞った**後**の本数なので、
    # 絞る前の本数も別列で残す (両者が違うと解釈を誤りやすい)
    results["n_true_edges_raw"] = len(set(true_edges))
    results["n_pred_edges_raw"] = len(set(pred_edges))
    results["n_pred_not_evaluable"] = n_excluded

    # SID は DAG 同士でしか定義されない。生物ネットワークの正解構造 (DREAM,
    # Sachs など) はフィードバックループを含むことがあるので、その場合は
    # SID だけを NA にして他の指標は計算する。
    sid = None
    sid_max = len(nodes) * (len(nodes) - 1)
    true_is_dag = is_dag(nodes, true_parents)
    pred_is_dag = is_dag(nodes, pred_parents)
    if args.skip_sid:
        if not args.quiet:
            log("SID: --skip-sid が指定されたのでスキップします")
    elif not true_is_dag:
        if not args.quiet:
            log("SID: 正解グラフが DAG ではない (閉路がある) のでスキップします")
    elif not pred_is_dag:
        if not args.quiet:
            log("SID: 学習グラフが DAG ではないのでスキップします")
    elif sid_max > args.max_sid_nodes * (args.max_sid_nodes - 1):
        if not args.quiet:
            log(f"SID: ノード数 {len(nodes)} が --max-sid-nodes "
                f"{args.max_sid_nodes} を超えるのでスキップします")
    else:
        sid = structural_intervention_distance(nodes, true_parents, pred_parents)
    results["sid"] = sid
    results["sid_max"] = sid_max
    results["sid_normalized"] = (sid / sid_max) if (sid is not None and sid_max) else None
    results["true_is_dag"] = int(true_is_dag)
    results["pred_is_dag"] = int(pred_is_dag)

    kl = None
    if true_bn is None:
        if not args.quiet:
            log("KL: --true-bif が無いのでスキップします")
    elif X is None:
        if not args.quiet:
            log("KL: --input が無いのでスキップします (学習側の CPT を推定できません)")
    elif true_bn.state_space_size() > args.max_states:
        if not args.quiet:
            log(f"KL: 状態空間 {true_bn.state_space_size()} が --max-states "
                f"{args.max_states} を超えるのでスキップします")
    else:
        bad = [(v, int(X[:, i].max())) for i, v in enumerate(nodes)
               if X[:, i].max() >= true_bn.card(v) or X[:, i].min() < 0]
        if bad:
            die(f"データの値が正解ネットワークの状態数を超えています: {bad[:3]}")
        kl = kl_divergence(true_bn, pred_parents, X, args.alpha)
    results["kl_divergence"] = kl

    # --- 出力 ---------------------------------------------------------------
    extra = {}
    for item in args.extra:
        if "=" not in item:
            die(f"--extra は KEY=VALUE 形式です: {item}")
        k, v = item.split("=", 1)
        extra[k] = v

    columns = list(extra.keys()) + METRIC_ORDER
    values = [extra[k] for k in extra] + [fmt(results.get(k)) for k in METRIC_ORDER]

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fp:
            fp.write("\t".join(columns) + "\n")
            fp.write("\t".join(values) + "\n")
    if args.append:
        os.makedirs(os.path.dirname(os.path.abspath(args.append)) or ".",
                    exist_ok=True)
        need_header = (not os.path.exists(args.append)
                       or os.path.getsize(args.append) == 0)
        with open(args.append, "a", encoding="utf-8") as fp:
            if need_header:
                fp.write("\t".join(columns) + "\n")
            fp.write("\t".join(values) + "\n")
    if args.json_out:
        payload = dict(extra)
        payload.update(results)
        with open(args.json_out, "w", encoding="utf-8") as fp:
            json.dump(payload, fp, ensure_ascii=False, indent=2)

    if not args.quiet:
        label = " ".join(f"{k}={v}" for k, v in extra.items())
        log(f"{label + ' | ' if label else ''}"
            f"SHD={results['shd']} "
            f"F1(dir)={results['f1_directed']:.3f} "
            f"F1(skel)={results['f1_skeleton']:.3f} "
            f"SID={fmt(sid)}/{sid_max} "
            f"KL={fmt(kl)}")
    if not args.out and not args.append and not args.json_out:
        print("\t".join(columns))
        print("\t".join(values))


if __name__ == "__main__":
    main()

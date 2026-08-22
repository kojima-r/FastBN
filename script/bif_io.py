#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bif_io.py
=========
BIF (Bayesian Interchange Format) 形式の離散ベイジアンネットワークを読み込み、

  * 構造 (エッジ) の取り出し
  * 祖先サンプリングによるデータ生成
  * 同時分布の厳密な展開 (小規模ネットワーク用)

を行うライブラリ兼コマンドラインツール。bnlearn Bayesian Network Repository
(https://www.bnlearn.com/bnrepository/) が配布する ``*.bif`` をそのまま扱える
(``.gz`` も透過的に読む)。

ライブラリとして:
    from bif_io import read_bif
    bn = read_bif("asia.bif")
    X  = bn.sample(1000, seed=1)      # (1000, D) の整数コード配列
    P  = bn.joint()                   # 同時分布 (各変数が 1 軸)

コマンドラインとして:
    # ネットワークの要約
    python3 bif_io.py info --bif asia.bif

    # データ生成 (fast_bn の入力形式 TSV) と正解エッジの書き出し
    python3 bif_io.py sample --bif asia.bif --n 1000 --seed 1 \
        --out data.tsv --out-edges true_edges.tsv

状態は BIF の宣言順に 0, 1, 2, ... と整数コード化する。出力 TSV のヘッダは
変数名 (BIF の宣言順) で、fast_bn のノード番号 = この列位置になる。
"""

import argparse
import gzip
import os
import re
import sys

import numpy as np


def log(*args):
    print("[bif_io]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[bif_io] エラー: {msg}")


# ---------------------------------------------------------------------------
# ネットワーク表現
# ---------------------------------------------------------------------------

class BayesNet:
    """離散ベイジアンネットワーク。

    names   : 変数名 (BIF の宣言順)
    states  : {変数名: [状態ラベル, ...]}
    parents : {変数名: [親の変数名, ...]}  (BIF の probability 宣言での順序)
    cpt     : {変数名: ndarray}
              shape = tuple(card[p] for p in parents) + (card[v],)
              最後の軸が自分の状態。親の軸は parents の順。
    """

    def __init__(self, names, states, parents, cpt):
        self.names = list(names)
        self.states = states
        self.parents = parents
        self.cpt = cpt
        self.index = {v: i for i, v in enumerate(self.names)}

    # --- 基本情報 ----------------------------------------------------------
    @property
    def n_vars(self):
        return len(self.names)

    def card(self, v):
        return len(self.states[v])

    @property
    def cards(self):
        return [self.card(v) for v in self.names]

    def edges(self):
        """(親, 子) の名前ペアを返す。"""
        out = []
        for v in self.names:
            for p in self.parents[v]:
                out.append((p, v))
        return out

    def n_params(self):
        """自由パラメータ数 (bnlearn の 'Number of parameters' と同じ数え方)。"""
        total = 0
        for v in self.names:
            q = 1
            for p in self.parents[v]:
                q *= self.card(p)
            total += q * (self.card(v) - 1)
        return total

    def topological_order(self):
        indeg = {v: len(self.parents[v]) for v in self.names}
        ready = [v for v in self.names if indeg[v] == 0]
        children = {v: [] for v in self.names}
        for v in self.names:
            for p in self.parents[v]:
                children[p].append(v)
        order = []
        while ready:
            v = ready.pop(0)
            order.append(v)
            for c in children[v]:
                indeg[c] -= 1
                if indeg[c] == 0:
                    ready.append(c)
        if len(order) != self.n_vars:
            die("ネットワークが DAG ではありません (閉路があります)")
        return order

    # --- サンプリング ------------------------------------------------------
    def sample(self, n, seed=None, rng=None):
        """祖先サンプリング。(n, D) の int32 配列 (列は self.names の順) を返す。"""
        if rng is None:
            rng = np.random.default_rng(seed)
        X = np.zeros((n, self.n_vars), dtype=np.int32)
        for v in self.topological_order():
            pa = self.parents[v]
            k = self.card(v)
            flat = self.cpt[v].reshape(-1, k)
            if pa:
                # 親設定のフラット添字 (親の宣言順・先頭が上位桁 = C 順)
                ridx = np.zeros(n, dtype=np.int64)
                for p in pa:
                    ridx = ridx * self.card(p) + X[:, self.index[p]]
            else:
                ridx = np.zeros(n, dtype=np.int64)
            probs = flat[ridx]                       # (n, k)
            cum = np.cumsum(probs, axis=1)
            cum = cum / cum[:, -1:]                  # 数値誤差の保険
            u = rng.random(n)[:, None]
            X[:, self.index[v]] = (u > cum).sum(axis=1).clip(0, k - 1)
        return X

    # --- 同時分布 ----------------------------------------------------------
    def state_space_size(self):
        size = 1
        for c in self.cards:
            size *= c
        return size

    def joint(self):
        """同時分布を ndarray (shape = cards, 軸の順 = self.names) で返す。"""
        shape = self.cards
        out = np.ones(shape, dtype=np.float64)
        for v in self.names:
            out = out * self._broadcast_cpt(v, self.parents[v], self.cpt[v])
        return out

    def _broadcast_cpt(self, v, parents, cpt):
        """CPT を全変数分の軸を持つ形に整形して返す (掛け算用)。"""
        dims = [self.index[p] for p in parents] + [self.index[v]]
        order = np.argsort(dims)
        arr = np.transpose(cpt, axes=tuple(order))
        shape = [1] * self.n_vars
        for d in dims:
            shape[d] = self.cards[d]
        return arr.reshape(shape)


# ---------------------------------------------------------------------------
# BIF パーサ
# ---------------------------------------------------------------------------

_COMMENT_RE = re.compile(r"//[^\n]*")
_VAR_RE = re.compile(r"\bvariable\s+(\"[^\"]+\"|[^\s{]+)\s*\{")
_PROB_RE = re.compile(r"\bprobability\s*\(([^)]*)\)\s*\{")
_TYPE_RE = re.compile(r"type\s+discrete\s*\[\s*(\d+)\s*\]\s*\{([^}]*)\}\s*;", re.S)
_NUM_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def _strip_quotes(tok):
    tok = tok.strip()
    if len(tok) >= 2 and tok[0] == tok[-1] and tok[0] in "\"'":
        return tok[1:-1]
    return tok


def _matching_block(text, open_pos):
    """text[open_pos] == '{' として、対応する '}' の直前までの中身を返す。"""
    depth = 0
    for i in range(open_pos, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[open_pos + 1:i], i
    die("BIF の括弧が閉じていません")


def _open_text(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fp:
            return fp.read()
    with open(path, encoding="utf-8", errors="replace") as fp:
        return fp.read()


def read_bif(path):
    """BIF ファイルを読み込んで BayesNet を返す。"""
    if not os.path.isfile(path):
        die(f"{path} がありません")
    text = _COMMENT_RE.sub("", _open_text(path))

    # --- variable ブロック --------------------------------------------------
    names, states = [], {}
    for m in _VAR_RE.finditer(text):
        name = _strip_quotes(m.group(1))
        body, _ = _matching_block(text, m.end() - 1)
        tm = _TYPE_RE.search(body)
        if not tm:
            log(f"警告: 変数 {name} に 'type discrete' がありません (スキップ)")
            continue
        labels = [_strip_quotes(s) for s in tm.group(2).split(",") if s.strip()]
        declared = int(tm.group(1))
        if declared != len(labels):
            log(f"警告: 変数 {name} の宣言状態数 {declared} と "
                f"ラベル数 {len(labels)} が違います (ラベル数を採用)")
        names.append(name)
        states[name] = labels
    if not names:
        die(f"{path} に variable 宣言が見つかりません")

    # --- probability ブロック -----------------------------------------------
    parents, cpt = {v: [] for v in names}, {}
    for m in _PROB_RE.finditer(text):
        head = [_strip_quotes(t) for t in re.split(r"[|,]", m.group(1)) if t.strip()]
        child, pa = head[0], head[1:]
        if child not in states:
            log(f"警告: 未宣言の変数 {child} の probability をスキップ")
            continue
        for p in pa:
            if p not in states:
                die(f"変数 {child} の親 {p} が宣言されていません")
        body, _ = _matching_block(text, m.end() - 1)
        parents[child] = pa
        cpt[child] = _parse_cpt(child, pa, states, body)

    missing = [v for v in names if v not in cpt]
    if missing:
        die(f"probability 宣言が無い変数があります: {', '.join(missing)}")

    return BayesNet(names, states, parents, cpt)


def _parse_cpt(child, pa, states, body):
    """probability ブロックの中身を CPT 配列に変換する。"""
    k = len(states[child])
    shape = tuple(len(states[p]) for p in pa) + (k,)
    table = np.zeros(shape, dtype=np.float64)
    filled = np.zeros(shape[:-1], dtype=bool) if pa else np.zeros((), dtype=bool)

    for stmt in body.split(";"):
        stmt = stmt.strip()
        if not stmt:
            continue
        if stmt.startswith("table"):
            vals = [float(x) for x in _NUM_RE.findall(stmt[len("table"):])]
            if pa:
                # 親つきの table 形式: 値は親設定を並べたもの (列優先) として扱う
                if len(vals) != table.size:
                    die(f"{child}: table の値数 {len(vals)} が "
                        f"期待値 {table.size} と一致しません")
                table[...] = np.asarray(vals).reshape(shape[:-1] + (k,), order="F")
                filled[...] = True
            else:
                if len(vals) != k:
                    die(f"{child}: table の値数 {len(vals)} が状態数 {k} と違います")
                table[...] = vals
                filled[...] = True
            continue
        if stmt.startswith("default"):
            log(f"警告: {child} の 'default' 指定は無視します")
            continue
        m = re.match(r"\(([^)]*)\)(.*)", stmt, re.S)
        if not m:
            continue
        labels = [_strip_quotes(t) for t in m.group(1).split(",") if t.strip()]
        if len(labels) != len(pa):
            die(f"{child}: 親設定 ({m.group(1)}) の要素数が親の数 {len(pa)} と違います")
        try:
            idx = tuple(states[p].index(lab) for p, lab in zip(pa, labels))
        except ValueError:
            die(f"{child}: 親設定 ({m.group(1)}) に未知の状態ラベルがあります")
        vals = [float(x) for x in _NUM_RE.findall(m.group(2))]
        if len(vals) != k:
            die(f"{child}: 親設定 ({m.group(1)}) の値数 {len(vals)} が "
                f"状態数 {k} と違います")
        table[idx] = vals
        filled[idx] = True

    if not np.all(filled):
        die(f"{child}: CPT に未指定の親設定があります")

    # 正規化 (BIF の丸め誤差を吸収)
    total = table.sum(axis=-1, keepdims=True)
    if np.any(total <= 0):
        die(f"{child}: 確率の合計が 0 の親設定があります")
    if np.max(np.abs(total - 1.0)) > 1e-3:
        log(f"警告: {child} の CPT の合計が 1 から離れています "
            f"(最大 {float(np.max(np.abs(total - 1.0))):.3g}) — 正規化します")
    return table / total


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def write_edges(path, edges):
    with open(path, "w", encoding="utf-8") as fp:
        for u, v in edges:
            fp.write(f"{u}\t{v}\n")


def cmd_info(args):
    bn = read_bif(args.bif)
    print(f"file\t{args.bif}")
    print(f"nodes\t{bn.n_vars}")
    print(f"arcs\t{len(bn.edges())}")
    print(f"parameters\t{bn.n_params()}")
    print(f"max_in_degree\t{max(len(bn.parents[v]) for v in bn.names)}")
    print(f"state_space\t{bn.state_space_size()}")
    print("variables\t" + ", ".join(f"{v}({bn.card(v)})" for v in bn.names))


def cmd_sample(args):
    bn = read_bif(args.bif)
    X = bn.sample(args.n, seed=args.seed)
    header = "\t".join(bn.names)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    np.savetxt(args.out, X, fmt="%d", delimiter="\t", header=header, comments="")
    log(f"{args.n} サンプル x {bn.n_vars} 変数 -> {args.out} (seed={args.seed})")
    if args.out_edges:
        write_edges(args.out_edges, bn.edges())
        log(f"正解エッジ {len(bn.edges())} 本 -> {args.out_edges}")
    # 出現しなかった状態があると fast_bn 側の基数が真の基数より小さくなる
    for j, v in enumerate(bn.names):
        seen = len(np.unique(X[:, j]))
        if seen < bn.card(v):
            log(f"注意: 変数 {v} は {bn.card(v)} 状態中 {seen} 状態しか出現しません "
                f"(n={args.n})")


def main():
    ap = argparse.ArgumentParser(description="BIF ネットワークの読み込みとサンプリング")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_info = sub.add_parser("info", help="ネットワークの要約を表示")
    p_info.add_argument("--bif", required=True)
    p_info.set_defaults(func=cmd_info)

    p_smp = sub.add_parser("sample", help="祖先サンプリングでデータを生成")
    p_smp.add_argument("--bif", required=True)
    p_smp.add_argument("--n", type=int, required=True, help="サンプル数")
    p_smp.add_argument("--seed", type=int, default=0, help="乱数シード")
    p_smp.add_argument("--out", required=True, help="出力 TSV (fast_bn 入力形式)")
    p_smp.add_argument("--out-edges", default=None, help="正解エッジの出力 TSV")
    p_smp.set_defaults(func=cmd_sample)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

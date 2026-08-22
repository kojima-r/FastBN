#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_sachs.py
================
Sachs らのフローサイトメトリー・タンパク質シグナル伝達データ
(https://zenodo.org/records/7681811, CC-BY-4.0) を fast_bn の入力形式に整える。

元データは 14 の実験条件ごとの CSV (行 = 単一細胞, 列 = 11 タンパク質の測定値) と、
既知のシグナル伝達経路 (GroundTruth.csv) から成る。本スクリプトは

  1. 条件を選んで結合し (--conditions)
  2. タンパク質ごとに離散化して (--bins / --method)
  3. fast_bn 入力 TSV + サンプル表 + 正解エッジ表

を書き出す。列順は全データセットで共通 (GroundTruth に現れる順ではなく、
CSV のヘッダ順) なので、fast_bn のノード番号 = 列位置 の対応が保たれる。

データセットの選び方 (--preset):
  obs   : 一般刺激のみ (cd3cd28)。介入なしの観測データに相当する古典的な設定
  int   : 阻害剤などの介入条件のみ
  all   : 14 条件すべてを結合 (最大サンプル数)

正解構造の注意: GroundTruth.csv は PKA <-> PIP3 のような相互作用を含むため
**DAG ではありません**。SID は DAG 同士でしか定義できないので、評価スクリプトは
SID を NA にして他の指標を計算します。
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "script"))
from discretize_matrix import discretize_column  # noqa: E402


def log(*args):
    print("[sachs]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[sachs] エラー: {msg}")


# 一般刺激のみ (介入なしとみなす) / 阻害剤などの介入あり
OBSERVATIONAL = ["cd3cd28"]
INTERVENTIONAL = [
    "cd3cd28_aktinhib", "cd3cd28_g0076", "cd3cd28_ly", "cd3cd28_psitect",
    "cd3cd28_u0126", "cd3cd28icam2_aktinhib", "cd3cd28icam2_g0076",
    "cd3cd28icam2_ly", "cd3cd28icam2_psit", "cd3cd28icam2_u0126",
    "cd3cd28_icam2", "pma", "b2camp",
]
PRESETS = {
    "obs": OBSERVATIONAL,
    "int": INTERVENTIONAL,
    "all": OBSERVATIONAL + INTERVENTIONAL,
}


def read_csv(path):
    with open(path, encoding="utf-8") as fp:
        header = [h.strip().strip('"') for h in fp.readline().rstrip("\n").split(",")]
        rows = []
        for line in fp:
            line = line.strip()
            if not line:
                continue
            rows.append([float(c) for c in line.split(",")])
    return header, np.asarray(rows, dtype=np.float64)


def main():
    ap = argparse.ArgumentParser(description="Sachs データを fast_bn 入力に整える")
    ap.add_argument("--data-dir", required=True,
                    help="展開した sachs.zip の 'Data Files' ディレクトリ")
    ap.add_argument("--ground-truth", required=True, help="GroundTruth.csv")
    ap.add_argument("--preset", choices=sorted(PRESETS), default="all")
    ap.add_argument("--conditions", default=None,
                    help="条件名をカンマ区切りで直接指定 (--preset より優先)")
    ap.add_argument("--bins", type=int, default=3)
    ap.add_argument("--method", choices=["quantile", "uniform"], default="quantile")
    ap.add_argument("--log2", action="store_true",
                    help="log2(x+1) 変換してから離散化する (蛍光強度は裾が長い)")
    ap.add_argument("--out", required=True, help="出力 TSV (fast_bn 入力)")
    ap.add_argument("--out-samples", default=None,
                    help="row_index / sample_id / group (= 条件名) の表")
    ap.add_argument("--out-edges", default=None, help="正解エッジ表の出力")
    ap.add_argument("--out-map", default=None, help="変数対応表の出力")
    args = ap.parse_args()

    conds = ([c.strip() for c in args.conditions.split(",") if c.strip()]
             if args.conditions else PRESETS[args.preset])

    # --- 条件ごとの CSV を結合 ---------------------------------------------
    names, blocks, groups = None, [], []
    for c in conds:
        path = os.path.join(args.data_dir, f"{c}.csv")
        if not os.path.isfile(path):
            die(f"{path} がありません (先に ./00download.sh を実行してください)")
        h, X = read_csv(path)
        if names is None:
            names = h
        elif h != names:
            die(f"{c}.csv の列順が他と違います: {h} != {names}")
        blocks.append(X)
        groups += [c] * X.shape[0]
    X = np.vstack(blocks)
    log(f"条件 {len(conds)} 件を結合: {X.shape[0]} 細胞 x {X.shape[1]} タンパク質")

    if args.log2:
        X = np.log2(np.clip(X, 0, None) + 1.0)

    # --- 離散化 (タンパク質ごと) --------------------------------------------
    codes = np.zeros(X.shape, dtype=np.int64)
    used = np.zeros(X.shape[1], dtype=int)
    for j in range(X.shape[1]):
        codes[:, j], used[j] = discretize_column(X[:, j], args.bins, args.method)
    if np.any(used < args.bins):
        log(f"注意: {int((used < args.bins).sum())} 変数が {args.bins} 段階に"
            f"分かれませんでした (最小 {int(used.min())})")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fp:
        fp.write("\t".join(names) + "\n")
        for row in codes:
            fp.write("\t".join(map(str, row)) + "\n")
    log(f"出力: {args.out}")

    if args.out_samples:
        with open(args.out_samples, "w", encoding="utf-8") as fp:
            fp.write("row_index\tsample_id\tgroup\n")
            counter = {}
            for i, g in enumerate(groups):
                counter[g] = counter.get(g, 0) + 1
                fp.write(f"{i}\t{g}_{counter[g]}\t{g}\n")
        log(f"サンプル表: {args.out_samples}")

    if args.out_map:
        with open(args.out_map, "w", encoding="utf-8") as fp:
            fp.write("index\tcolumn_name\tgene_id\tgene_name\tvariance\t"
                     "detected_frac\tused_levels\twhitelisted\n")
            for j, nm in enumerate(names):
                fp.write(f"{j}\t{nm}\t{nm}\t{nm}\t{np.var(X[:, j]):.6g}\t"
                         f"{float(np.mean(codes[:, j] > 0)):.6g}\t{used[j]}\t0\n")
        log(f"対応表: {args.out_map}")

    # --- 正解エッジ ---------------------------------------------------------
    if args.out_edges:
        known = set(names)
        edges, skipped = [], []
        with open(args.ground_truth, encoding="utf-8") as fp:
            head = fp.readline()
            for line in fp:
                line = line.strip()
                if not line:
                    continue
                parts = [p.strip().strip('"') for p in line.split(",")]
                if len(parts) < 2:
                    continue
                u, v = parts[0], parts[1]
                (edges if (u in known and v in known) else skipped).append((u, v))
        if skipped:
            log(f"警告: データに無いタンパク質のエッジを {len(skipped)} 本除外: "
                f"{skipped[:3]}")
        with open(args.out_edges, "w", encoding="utf-8") as fp:
            for u, v in edges:
                fp.write(f"{u}\t{v}\n")
        # 双方向のペア (= 閉路) があるか報告する
        st = set(edges)
        cyc = sorted({tuple(sorted(e)) for e in st if (e[1], e[0]) in st})
        log(f"正解エッジ {len(edges)} 本 -> {args.out_edges}"
            + (f" (相互作用 {len(cyc)} 組を含むため DAG ではありません: {cyc})"
               if cyc else ""))


if __name__ == "__main__":
    main()

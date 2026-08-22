#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
discretize_matrix.py
====================
連続値の数値行列 (行 = サンプル, 列 = 変数, 1 行目 = 変数名) を fast_bn の入力形式
(整数コードの TSV) に変換する汎用ツール。

`preprocess_expr.py` がバルク RNA 用に正規化・フィルタまで面倒を見るのに対し、
こちらは「すでに使える形の数値行列を離散化するだけ」の薄いツールで、
ベンチマーク用データ (DREAM, Sachs など) に使う。

    python3 discretize_matrix.py --input expr.tsv --out disc.tsv \
        --bins 3 --method quantile --out-map var_map.tsv

主なオプション:
  --bins N          離散化の段階数 (既定 3)
  --method          quantile (等頻度; 既定) / uniform (等幅)
  --log2            log2(x + --pseudocount) 変換してから離散化する
  --max-vars N      分散上位 N 変数だけ残す (0 = 制限なし)
  --keep-vars FILE  1 行 1 変数名。--max-vars で削られても必ず残す
  --drop-constant   分散 0 (値が 1 種類) の列を落とす
  --transpose       入力が「行 = 変数, 列 = サンプル」の場合に指定

等頻度離散化では、同じ値が多い変数 (ゼロが多いなど) で分位点が重複し、実際の
段階数が --bins より少なくなることがある。その場合は警告し、使われた段階数を
var_map の used_levels 列に出す。fast_bn 側の基数は「観測された最大値 + 1」なので、
コードは 0 から詰めて振り直す。
"""

import argparse
import os
import sys

import numpy as np


def log(*args):
    print("[discretize]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[discretize] エラー: {msg}")


def read_matrix(path, transpose=False):
    """数値行列を (変数名リスト, ndarray[サンプル, 変数]) で返す。"""
    if not os.path.isfile(path):
        die(f"{path} がありません")
    sep = "," if os.path.splitext(path)[1].lower() == ".csv" else "\t"
    with open(path, encoding="utf-8") as fp:
        header = [h.strip().strip('"') for h in
                  fp.readline().rstrip("\n").split(sep)]
        rows = []
        for lineno, line in enumerate(fp, start=2):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cells = line.split(sep)
            if len(cells) == len(header) + 1:
                cells = cells[1:]          # 先頭が行名の場合は落とす
            if len(cells) != len(header):
                die(f"{path}:{lineno} の列数 {len(cells)} がヘッダ {len(header)} と違います")
            try:
                rows.append([float(c) if c.strip() not in ("", "NA", "NaN")
                             else np.nan for c in cells])
            except ValueError:
                die(f"{path}:{lineno} に数値でない値があります")
    if not rows:
        die(f"{path} にデータ行がありません")
    X = np.asarray(rows, dtype=np.float64)
    if transpose:
        X = X.T
        header = [f"S{i + 1}" for i in range(X.shape[1])]
        die("--transpose は変数名の対応が取れないため未対応です "
            "(行=サンプル・列=変数に整形してから渡してください)")
    return header, X


def discretize_column(x, bins, method):
    """1 列を 0..k-1 の整数コードにする。実際に使われた段階数も返す。"""
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return np.zeros(x.shape, dtype=np.int64), 1
    if method == "uniform":
        lo, hi = float(finite.min()), float(finite.max())
        if hi <= lo:
            return np.zeros(x.shape, dtype=np.int64), 1
        edges = np.linspace(lo, hi, bins + 1)[1:-1]
    else:  # quantile
        qs = np.linspace(0, 100, bins + 1)[1:-1]
        edges = np.percentile(finite, qs)
        edges = np.unique(edges)           # 分位点が重複したら段階数が減る
    codes = np.digitize(np.nan_to_num(x, nan=float(np.median(finite))), edges,
                        right=False).astype(np.int64)
    # 出現したコードを 0 から詰め直す (fast_bn の基数 = 最大値+1 のため)
    uniq = np.unique(codes)
    remap = {c: i for i, c in enumerate(uniq)}
    codes = np.vectorize(remap.get)(codes).astype(np.int64)
    return codes, len(uniq)


def main():
    ap = argparse.ArgumentParser(description="連続値行列を離散化して fast_bn 入力にする")
    ap.add_argument("--input", required=True, help="数値行列 (行=サンプル, 列=変数)")
    ap.add_argument("--out", required=True, help="出力 TSV (整数コード)")
    ap.add_argument("--out-map", default=None, help="変数対応表 TSV")
    ap.add_argument("--bins", type=int, default=3, help="離散化の段階数 (既定 3)")
    ap.add_argument("--method", choices=["quantile", "uniform"], default="quantile")
    ap.add_argument("--log2", action="store_true", help="log2(x + pseudocount) 変換")
    ap.add_argument("--pseudocount", type=float, default=1.0)
    ap.add_argument("--max-vars", type=int, default=0,
                    help="分散上位 N 変数だけ残す (0 = 制限なし)")
    ap.add_argument("--keep-vars", default=None,
                    help="必ず残す変数名のファイル (1 行 1 名)")
    ap.add_argument("--drop-constant", action="store_true",
                    help="値が 1 種類しかない列を落とす")
    ap.add_argument("--transpose", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    names, X = read_matrix(args.input, args.transpose)
    n_samples, n_vars = X.shape
    if not args.quiet:
        log(f"入力: {args.input} ({n_samples} サンプル x {n_vars} 変数)")

    if args.log2:
        X = np.log2(np.clip(X, 0, None) + args.pseudocount)

    # --- 変数の選択 ---------------------------------------------------------
    keep = set()
    if args.keep_vars and os.path.isfile(args.keep_vars):
        with open(args.keep_vars, encoding="utf-8") as fp:
            keep = {l.strip().strip('"') for l in fp if l.strip()}
        present = keep & set(names)
        if not args.quiet:
            log(f"必ず残す変数: {len(present)} / 指定 {len(keep)}")
        keep = present

    var = np.nanvar(X, axis=0)
    sel = np.arange(n_vars)
    if args.max_vars and args.max_vars < n_vars:
        forced = np.array([i for i, nm in enumerate(names) if nm in keep], dtype=int)
        n_extra = max(0, args.max_vars - len(forced))
        rest = np.array([i for i in range(n_vars) if i not in set(forced.tolist())],
                        dtype=int)
        rest = rest[np.argsort(-var[rest])][:n_extra]
        sel = np.sort(np.concatenate([forced, rest])) if len(forced) else np.sort(rest)
        if not args.quiet:
            log(f"分散上位で {len(sel)} 変数に絞りました "
                f"(必ず残す {len(forced)} を含む)")

    names = [names[i] for i in sel]
    X = X[:, sel]

    # --- 離散化 -------------------------------------------------------------
    codes = np.zeros(X.shape, dtype=np.int64)
    used = np.zeros(X.shape[1], dtype=int)
    for j in range(X.shape[1]):
        codes[:, j], used[j] = discretize_column(X[:, j], args.bins, args.method)

    if args.drop_constant:
        ok = used > 1
        if not np.all(ok) and not args.quiet:
            log(f"値が 1 種類の列を {int((~ok).sum())} 個落としました")
        names = [nm for nm, k in zip(names, ok) if k]
        codes, used = codes[:, ok], used[ok]

    short = int((used < args.bins).sum())
    if short and not args.quiet:
        log(f"注意: {short} 変数は分位点の重複で {args.bins} 段階に分かれませんでした "
            f"(最小 {int(used.min())} 段階)")

    # --- 出力 ---------------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fp:
        fp.write("\t".join(names) + "\n")
        for row in codes:
            fp.write("\t".join(map(str, row)) + "\n")
    if not args.quiet:
        log(f"出力: {args.out} ({codes.shape[0]} サンプル x {codes.shape[1]} 変数, "
            f"{args.bins} 段階 {args.method})")

    if args.out_map:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_map)) or ".",
                    exist_ok=True)
        with open(args.out_map, "w", encoding="utf-8") as fp:
            fp.write("index\tcolumn_name\tgene_id\tgene_name\tvariance\t"
                     "detected_frac\tused_levels\twhitelisted\n")
            for j, nm in enumerate(names):
                col = X[:, j]
                fp.write(f"{j}\t{nm}\t{nm}\t{nm}\t{np.nanvar(col):.6g}\t"
                         f"{float(np.mean(codes[:, j] > 0)):.6g}\t{used[j]}\t"
                         f"{1 if nm in keep else 0}\n")
        if not args.quiet:
            log(f"対応表: {args.out_map}")


if __name__ == "__main__":
    main()

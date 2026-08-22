#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_data.py
===============
example_sc 専用のデータ準備スクリプト。

この例題のデータ (Tabula Muris Senis 由来の単一細胞発現量) は **すでに離散化済み**
なので、汎用パイプラインの前処理 (``script/preprocess_expr.py``: 正規化 -> log ->
フィルタ -> 離散化) は不要である。その代わり、離散化済み行列を
``script/*.sh`` がそのまま扱える形に整える:

  1. ``data/expr_disc.tsv``  : 選んだ行列 (all_disc*.tsv) へのシンボリックリンク
                               (作れない環境ではコピー)
  2. ``data/var_map.tsv``    : 列インデックス <-> 遺伝子名 の対応表
                               (index / column_name / gene_id / gene_name /
                                variance / detected_frac / used_levels / whitelisted)
  3. ``data/samples.tsv``    : 行番号 <-> サンプル ID / 群ラベル
                               (row_index / sample_id / group; **群 = 組織**)

``all_disc<N>.tsv`` の行は ``tissue/`` 以下の組織別ファイルをファイル名順に連結した
ものと一致する (本スクリプトが検証する)。この対応を使って「組織」を群ラベルとする
samples.tsv を生成し、``script/importance_groups.sh`` による組織別エッジ重要度を
可能にしている。

使い方 (通常は 01prepare.sh 経由で呼ばれる):
  python3 prepare_data.py --matrix data_bbknn_r_tissue_disc/all_disc100.tsv \
      --tissue-dir data_bbknn_r_tissue_disc/tissue \
      --outdir run_bbknn_bin100/data
"""

import argparse
import os
import re
import shutil
import sys

# 組織ファイル名から取り除くアッセイ名のサフィックス
# (例: abdominal_aorta.smartseq.tsv -> abdominal_aorta)
ASSAY_SUFFIX_RE = re.compile(r"\.(smartseq|smart_seq|tenx|10x|droplet|facs)$",
                             re.IGNORECASE)


def log(*args):
    print("[prepare]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[prepare] エラー: {msg}")


def read_header(path):
    with open(path, encoding="utf-8") as fp:
        line = fp.readline()
    if not line:
        die(f"{path} が空です")
    return line.rstrip("\n").split("\t")


def tissue_label(path):
    """組織ファイルのパスから群ラベルを作る。"""
    base = os.path.basename(path)
    base = os.path.splitext(base)[0]          # .tsv を除去
    return ASSAY_SUFFIX_RE.sub("", base)      # .smartseq などを除去


def scan_matrix(path):
    """行列を 1 度だけ走査して列ごとの統計量を集める。

    戻り値: (n_rows, sums, sumsq, nonzero, levels)
      levels[j] は列 j に出現した値の集合
    """
    header = read_header(path)
    d = len(header)
    sums = [0.0] * d
    sumsq = [0.0] * d
    nonzero = [0] * d
    levels = [set() for _ in range(d)]
    n_rows = 0

    with open(path, encoding="utf-8") as fp:
        fp.readline()  # ヘッダ
        for lineno, line in enumerate(fp, start=2):
            line = line.rstrip("\n")
            if line == "":
                continue
            cells = line.split("\t")
            if len(cells) != d:
                die(f"{path}:{lineno} の列数 ({len(cells)}) がヘッダ ({d}) と違います")
            n_rows += 1
            for j, c in enumerate(cells):
                try:
                    v = int(c)
                except ValueError:
                    die(f"{path}:{lineno} 列 {j + 1} が整数ではありません: '{c}'\n"
                        "  この例題の入力は離散化済み (整数コード) である必要があります。")
                sums[j] += v
                sumsq[j] += v * v
                if v != 0:
                    nonzero[j] += 1
                levels[j].add(v)

    if n_rows == 0:
        die(f"{path} にデータ行がありません")
    return header, n_rows, sums, sumsq, nonzero, levels


def build_samples(tissue_dir, matrix_header, n_rows):
    """組織別ファイルから (sample_id, group) のリストを作る。

    行列の行順 = 組織ファイルをファイル名順に連結した順、という前提を検証する。
    検証に失敗した場合は None を返す (群別解析だけがスキップされる)。
    """
    if not tissue_dir or not os.path.isdir(tissue_dir):
        log(f"警告: 組織ディレクトリがありません ({tissue_dir})。"
            " samples.tsv は作成しません (組織別解析はスキップされます)")
        return None

    files = sorted(os.path.join(tissue_dir, f)
                   for f in os.listdir(tissue_dir) if f.endswith(".tsv"))
    if not files:
        log(f"警告: {tissue_dir} に .tsv がありません。samples.tsv は作成しません")
        return None

    # 列順の検証: 組織ファイルの先頭 D 列が行列のヘッダと一致するか
    d = len(matrix_header)
    head0 = read_header(files[0])
    if len(head0) < d or head0[:d] != matrix_header:
        log(f"警告: {files[0]} の列順が行列と一致しません。"
            " samples.tsv は作成しません (組織別解析はスキップされます)")
        return None

    rows = []
    for path in files:
        label = tissue_label(path)
        n = 0
        with open(path, encoding="utf-8") as fp:
            fp.readline()
            for line in fp:
                if line.strip() != "":
                    n += 1
        for i in range(n):
            rows.append((f"{label}_{i + 1}", label))

    if len(rows) != n_rows:
        log(f"警告: 組織ファイルの合計行数 ({len(rows)}) が行列の行数 ({n_rows}) と"
            " 一致しません。samples.tsv は作成しません")
        return None

    n_groups = len(set(g for _, g in rows))
    log(f"組織 (群): {n_groups} 件 / 合計 {len(rows)} サンプル")
    return rows


def link_or_copy(src, dst):
    """src -> dst の相対シンボリックリンクを作る。失敗したらコピーする。"""
    if os.path.islink(dst) or os.path.exists(dst):
        os.remove(dst)
    rel = os.path.relpath(os.path.abspath(src), os.path.dirname(os.path.abspath(dst)))
    try:
        os.symlink(rel, dst)
        log(f"入力をリンク: {dst} -> {rel}")
    except OSError:
        shutil.copyfile(src, dst)
        log(f"入力をコピー: {dst} ({os.path.getsize(dst) / 1e6:.1f} MB)")


def main():
    ap = argparse.ArgumentParser(
        description="離散化済み行列を script/ パイプラインの入力形式に整える")
    ap.add_argument("--matrix", required=True,
                    help="離散化済み行列 (all_disc*.tsv; 行=サンプル, 列=遺伝子)")
    ap.add_argument("--tissue-dir", default=None,
                    help="組織別ファイルのディレクトリ (群ラベルの生成に使う)")
    ap.add_argument("--outdir", required=True, help="出力先 (通常 <RUNDIR>/data)")
    ap.add_argument("--targets", default="",
                    help="注目遺伝子 (カンマ区切り)。var_map の whitelisted 列と "
                         "--out-targets のファイルに反映される")
    ap.add_argument("--out-targets", default=None,
                    help="注目遺伝子リストの出力先 (target_genes.txt)")
    args = ap.parse_args()

    if not os.path.isfile(args.matrix):
        die(f"行列 {args.matrix} がありません。先に ./00download.sh を実行してください。")

    os.makedirs(args.outdir, exist_ok=True)

    log(f"行列を走査中: {args.matrix}")
    header, n_rows, sums, sumsq, nonzero, levels = scan_matrix(args.matrix)
    d = len(header)
    log(f"{n_rows} サンプル x {d} 遺伝子")

    # --- 1. 入力 (expr_disc.tsv) ------------------------------------------
    input_path = os.path.join(args.outdir, "expr_disc.tsv")
    link_or_copy(args.matrix, input_path)

    # --- 2. 注目遺伝子 ------------------------------------------------------
    wanted = [t.strip() for t in args.targets.split(",") if t.strip()]
    present = [t for t in wanted if t in header]
    missing = [t for t in wanted if t not in header]
    if missing:
        log(f"警告: 注目遺伝子のうち {len(missing)} 件はこのデータに含まれません: "
            f"{', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}")
    target_set = set(present)
    if args.out_targets:
        if present:
            with open(args.out_targets, "w", encoding="utf-8") as fp:
                fp.write("\n".join(present) + "\n")
            log(f"注目遺伝子 {len(present)} 件 -> {args.out_targets}")
        elif os.path.exists(args.out_targets):
            os.remove(args.out_targets)  # 古い設定の残骸を消す

    # --- 3. var_map.tsv -----------------------------------------------------
    var_map = os.path.join(args.outdir, "var_map.tsv")
    with open(var_map, "w", encoding="utf-8") as fp:
        fp.write("index\tcolumn_name\tgene_id\tgene_name\tvariance\t"
                 "detected_frac\tused_levels\twhitelisted\n")
        for j, name in enumerate(header):
            mean = sums[j] / n_rows
            var = max(0.0, sumsq[j] / n_rows - mean * mean)
            fp.write(f"{j}\t{name}\t{name}\t{name}\t{var:.6g}\t"
                     f"{nonzero[j] / n_rows:.6g}\t{len(levels[j])}\t"
                     f"{1 if name in target_set else 0}\n")
    log(f"対応表を出力: {var_map}")

    const_cols = sum(1 for s in levels if len(s) < 2)
    if const_cols:
        log(f"注意: 値が 1 種類しかない列が {const_cols} 件あります "
            "(これらのノードにはエッジが張られません)")

    # --- 4. samples.tsv -----------------------------------------------------
    samples = build_samples(args.tissue_dir, header, n_rows)
    samples_path = os.path.join(args.outdir, "samples.tsv")
    if samples is None:
        if os.path.exists(samples_path):
            os.remove(samples_path)
    else:
        with open(samples_path, "w", encoding="utf-8") as fp:
            fp.write("row_index\tsample_id\tgroup\n")
            for i, (sid, grp) in enumerate(samples):
                fp.write(f"{i}\t{sid}\t{grp}\n")
        log(f"サンプル表を出力: {samples_path}")


if __name__ == "__main__":
    main()

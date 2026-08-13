#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_groups.py
==============
離散化済み入力 TSV (preprocess_expr.py の出力; 行=サンプル, 列=遺伝子) を
実験群 (サブセット) ごとに分割する汎用スクリプト。各群ファイルはヘッダ
(遺伝子名) を保持し、その群に属するサンプル行のみを含む。fast_bn の
`--score-dataset` に渡して群別のエッジ重要度を計算するために使う。

群の指定方法は 3 通り (上から優先):
  1. --samples data/samples.tsv        : preprocess_expr.py の --out-samples 出力
                                         (row_index / sample_id / group 列) を使う
  2. --meta sample_meta.tsv            : サンプル ID と群ラベルの表を直接指定
                                         (入力の行順とメタデータのサンプル順が
                                          一致している前提)
  3. --labels / --sizes                : 「先頭から n 件ずつ」を手で指定

いずれの方法でも、**群はサンプル行が連続していなくても構わない**
(行番号を群ごとに集めて出力する)。

出力:
  <outdir>/<prefix>_g{N}_<label>.tsv   : 群ごとの score-dataset
  <outdir>/groups_manifest.tsv         : group_no / label / n_samples / file
"""

import argparse
import os
import sys


def log(*args):
    print("[make_groups]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[make_groups] エラー: {msg}")


def sanitize_label(label):
    """群ラベルをファイル名に使える形に正規化する。"""
    out = [ch if (ch.isalnum() or ch in "-.") else "_" for ch in str(label)]
    return "".join(out).strip("_") or "NA"


def read_delim(path):
    """TSV/CSV を [header, row, ...] のリストで返す。"""
    sep = "," if os.path.splitext(path)[1].lower() == ".csv" else "\t"
    rows = []
    with open(path, encoding="utf-8") as fp:
        for line in fp:
            line = line.rstrip("\n")
            if line.strip() == "":
                continue
            rows.append([c.strip() for c in line.split(sep)])
    if not rows:
        die(f"{path} が空です")
    return rows


def col_index(header, spec, default_idx, what):
    """列指定 (列名 or 0 始まり位置) を列番号に解決する。"""
    if spec is None or spec == "":
        return default_idx
    if spec in header:
        return header.index(spec)
    try:
        return int(spec)
    except ValueError:
        die(f"{what} に指定した列 '{spec}' が {header} に見つかりません")


def groups_from_samples_file(path, n_rows, group_col):
    """samples.tsv (row_index / sample_id / group) から (label, [行番号]) を作る。"""
    rows = read_delim(path)
    header, body = rows[0], rows[1:]
    gi = col_index(header, group_col, 2 if len(header) > 2 else len(header) - 1,
                   "--group-col")
    ri = header.index("row_index") if "row_index" in header else None
    order, members = [], {}
    for i, r in enumerate(body):
        if len(r) <= gi:
            die(f"{path} の {i + 2} 行目に群列がありません: {r}")
        label = sanitize_label(r[gi])
        row_no = int(r[ri]) if ri is not None and len(r) > ri else i
        if label not in members:
            members[label] = []
            order.append(label)
        members[label].append(row_no)
    if sum(len(v) for v in members.values()) != n_rows:
        die(f"{path} のサンプル数 ({sum(len(v) for v in members.values())}) が "
            f"入力の行数 ({n_rows}) と一致しません")
    return [(lab, members[lab]) for lab in order]


def groups_from_meta(path, n_rows, sample_col, group_col):
    """メタデータ (サンプル ID + 群) の行順が入力の行順と一致する前提で群を作る。"""
    rows = read_delim(path)
    header, body = rows[0], rows[1:]
    si = col_index(header, sample_col, 0, "--sample-col")
    gi = col_index(header, group_col, 1 if len(header) > 1 else 0, "--group-col")
    if len(body) != n_rows:
        die(f"{path} の行数 ({len(body)}) が入力のサンプル数 ({n_rows}) と一致しません。"
            " preprocess_expr.py の --out-samples 出力 (--samples) を使うか、"
            " メタデータの行順を入力に合わせてください。")
    order, members = [], {}
    for i, r in enumerate(body):
        label = sanitize_label(r[gi])
        if label not in members:
            members[label] = []
            order.append(label)
        members[label].append(i)
    log(f"メタデータ列: サンプル='{header[si]}', 群='{header[gi]}'")
    return [(lab, members[lab]) for lab in order]


def groups_from_sizes(labels, sizes, n_rows):
    """先頭から sizes 件ずつ切って群を作る (従来方式)。"""
    labs = [sanitize_label(s) for s in labels.split(",") if s.strip()]
    szs = [int(s) for s in sizes.split(",") if s.strip()]
    if len(labs) != len(szs):
        die(f"--labels ({len(labs)}) と --sizes ({len(szs)}) の数が一致しません")
    if sum(szs) != n_rows:
        die(f"--sizes の合計 ({sum(szs)}) が入力のサンプル数 ({n_rows}) と一致しません")
    out, start = [], 0
    for lab, n in zip(labs, szs):
        out.append((lab, list(range(start, start + n))))
        start += n
    return out


def main():
    ap = argparse.ArgumentParser(
        description="離散化入力を実験群 (サブセット) ごとに分割する")
    ap.add_argument("--input", required=True,
                    help="離散化済み入力 TSV (行=サンプル, 列=遺伝子, ヘッダ付き)")
    ap.add_argument("--outdir", required=True, help="群ファイルの出力先")
    ap.add_argument("--prefix", default="expr", help="群ファイル名の接頭辞")
    ap.add_argument("--samples", default=None,
                    help="preprocess_expr.py --out-samples の出力 "
                         "(row_index / sample_id / group)")
    ap.add_argument("--meta", default=None,
                    help="サンプル ID と群ラベルの表 (行順が入力と一致している前提)")
    ap.add_argument("--sample-col", default=None, help="メタデータのサンプル ID 列")
    ap.add_argument("--group-col", default=None, help="群ラベルの列 (既定 'group')")
    ap.add_argument("--labels", default=None, help="群ラベル (カンマ区切り; 手動指定)")
    ap.add_argument("--sizes", default=None, help="各群のサンプル数 (カンマ区切り)")
    ap.add_argument("--min-samples", type=int, default=2,
                    help="この件数未満の群は警告してスキップする")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        die(f"入力 {args.input} がありません")
    with open(args.input, encoding="utf-8") as fp:
        header = fp.readline()
        rows = [line for line in fp if line.strip() != ""]
    n_rows = len(rows)
    log(f"入力: {args.input} ({n_rows} サンプル, {len(header.split(chr(9)))} 列)")

    if args.samples and os.path.exists(args.samples):
        groups = groups_from_samples_file(args.samples, n_rows, args.group_col)
        log(f"群の定義元: {args.samples}")
    elif args.meta:
        groups = groups_from_meta(args.meta, n_rows, args.sample_col, args.group_col)
        log(f"群の定義元: {args.meta}")
    elif args.labels and args.sizes:
        groups = groups_from_sizes(args.labels, args.sizes, n_rows)
        log("群の定義元: --labels / --sizes")
    else:
        die("群の指定がありません。--samples / --meta / (--labels と --sizes) の"
            "いずれかを指定してください。")

    if len(groups) < 2:
        log(f"警告: 群が {len(groups)} 件しかありません (群間比較には 2 群以上が必要)")

    os.makedirs(args.outdir, exist_ok=True)
    manifest = []
    gno = 0
    for label, idxs in groups:
        if len(idxs) < args.min_samples:
            log(f"警告: 群 '{label}' は {len(idxs)} サンプルのみ "
                f"(< --min-samples {args.min_samples}) のためスキップ")
            continue
        gno += 1
        fname = f"{args.prefix}_g{gno}_{label}.tsv"
        path = os.path.join(args.outdir, fname)
        with open(path, "w", encoding="utf-8") as out:
            out.write(header)
            out.writelines(rows[i] for i in idxs)
        manifest.append((gno, label, len(idxs), path))
        log(f"群 {gno} '{label}': {len(idxs)} サンプル -> {path}")

    if not manifest:
        die("出力できる群がありませんでした (--min-samples を下げてください)")

    man_path = os.path.join(args.outdir, "groups_manifest.tsv")
    with open(man_path, "w", encoding="utf-8") as out:
        out.write("group_no\tlabel\tn_samples\tfile\n")
        for gno, label, n, path in manifest:
            out.write(f"{gno}\t{label}\t{n}\t{path}\n")
    log(f"マニフェスト: {man_path}")


if __name__ == "__main__":
    main()

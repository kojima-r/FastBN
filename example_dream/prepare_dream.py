#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_dream.py
================
DREAM チャレンジのネットワーク推定データを fast_bn の入力形式に整える。

対応するデータセット:

  dream4  : DREAM4 in silico network challenge (GeneNetWeaver 配布の zip)
            size10 x 5 と size100 x 5 の人工ネットワーク。
            発現量は複数の摂動実験 (multifactorial / knockouts / timeseries /
            wildtype) を縦に結合して 1 つの行列にする。
  dream5  : DREAM5 network inference challenge (Zenodo 17854236 の
            1_Challenge_Data_Supplement.zip)
            Network1 (in silico) / Network3 (E. coli) / Network4 (S. cerevisiae)。
            Network2 (S. aureus) は gold standard が採点に使われていないので除く。
  hpn     : HPN-DREAM breast cancer network inference challenge
            (Synapse syn1720047)。Synapse は認証が要るため自動取得できない。
            手元に展開したディレクトリがある場合のみ処理する。

各データセットについて次を書き出す:

  <out-dir>/<network>.tsv              fast_bn 入力 (整数コード, 行=サンプル)
  <out-dir>/<network>_varmap.tsv       列 -> 遺伝子名
  <truth-dir>/<network>_edges.tsv      正解エッジ (u -> v)
  <truth-dir>/<network>_pairs.tsv      評価対象ペア (gold standard が判定した組)

正解構造の注意: 遺伝子制御ネットワークはフィードバックループを含むので、
**DAG とは限りません**。評価スクリプトはその場合 SID を NA にします。
また DREAM5 の gold standard は TF x 遺伝子の一部ペアしか判定していないため、
`_pairs.tsv` を `evaluate_structure.py --eval-pairs` に渡して、判定対象の
ペアだけで Precision / Recall を計算します。
"""

import argparse
import os
import sys
import zipfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "script"))
from discretize_matrix import discretize_column  # noqa: E402


def log(*args):
    print("[dream]", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[dream] エラー: {msg}")


def dec(b):
    return b.decode("utf-8", errors="replace")


# ---------------------------------------------------------------------------
# 共通の書き出し
# ---------------------------------------------------------------------------

def write_dataset(name, gene_names, X, bins, method, out_dir, log2=False):
    """連続値行列を離散化して fast_bn 入力として書き出す。"""
    os.makedirs(out_dir, exist_ok=True)
    if log2:
        X = np.log2(np.clip(X, 0, None) + 1.0)
    codes = np.zeros(X.shape, dtype=np.int64)
    used = np.zeros(X.shape[1], dtype=int)
    for j in range(X.shape[1]):
        codes[:, j], used[j] = discretize_column(X[:, j], bins, method)

    keep = used > 1                      # 定数列は fast_bn で意味を持たない
    if not np.all(keep):
        log(f"  {name}: 値が 1 種類の列を {int((~keep).sum())} 個落とします")
    gene_names = [g for g, k in zip(gene_names, keep) if k]
    codes, used, X = codes[:, keep], used[keep], X[:, keep]

    path = os.path.join(out_dir, f"{name}.tsv")
    with open(path, "w", encoding="utf-8") as fp:
        fp.write("\t".join(gene_names) + "\n")
        for row in codes:
            fp.write("\t".join(map(str, row)) + "\n")
    with open(os.path.join(out_dir, f"{name}_varmap.tsv"), "w",
              encoding="utf-8") as fp:
        fp.write("index\tcolumn_name\tgene_id\tgene_name\tvariance\t"
                 "detected_frac\tused_levels\twhitelisted\n")
        for j, g in enumerate(gene_names):
            fp.write(f"{j}\t{g}\t{g}\t{g}\t{np.nanvar(X[:, j]):.6g}\t"
                     f"{float(np.mean(codes[:, j] > 0)):.6g}\t{used[j]}\t0\n")
    log(f"  {name}: {codes.shape[0]} サンプル x {codes.shape[1]} 遺伝子 -> {path}")
    return set(gene_names)


def write_truth(name, edges, pairs, present, truth_dir):
    """正解エッジと評価対象ペアを、解析対象に残った遺伝子だけに絞って書き出す。"""
    os.makedirs(truth_dir, exist_ok=True)
    e = [(u, v) for u, v in edges if u in present and v in present and u != v]
    epath = os.path.join(truth_dir, f"{name}_edges.tsv")
    with open(epath, "w", encoding="utf-8") as fp:
        for u, v in e:
            fp.write(f"{u}\t{v}\n")
    ppath = None
    if pairs is not None:
        p = [(u, v) for u, v in pairs if u in present and v in present and u != v]
        ppath = os.path.join(truth_dir, f"{name}_pairs.tsv")
        with open(ppath, "w", encoding="utf-8") as fp:
            for u, v in p:
                fp.write(f"{u}\t{v}\n")
        log(f"  {name}: 正解 {len(e)} 辺 / 評価対象 {len(p)} ペア")
    else:
        log(f"  {name}: 正解 {len(e)} 辺")
    st = set(e)
    cyc = {tuple(sorted(x)) for x in st if (x[1], x[0]) in st}
    if cyc:
        log(f"  {name}: 相互作用 {len(cyc)} 組を含むため DAG ではありません")
    return epath, ppath


def select_columns(gene_names, X, max_vars, priority=()):
    """分散上位で列を絞る。priority の遺伝子 (TF など) は優先して残す。"""
    if not max_vars or max_vars >= X.shape[1]:
        return gene_names, X
    pri = set(priority)
    idx_pri = [i for i, g in enumerate(gene_names) if g in pri]
    var = np.nanvar(X, axis=0)
    rest = [i for i in range(X.shape[1]) if i not in set(idx_pri)]
    rest.sort(key=lambda i: -var[i])
    # TF が多すぎる場合も分散で切る
    if len(idx_pri) > max_vars:
        idx_pri.sort(key=lambda i: -var[i])
        idx_pri = idx_pri[:max_vars]
        sel = sorted(idx_pri)
    else:
        sel = sorted(idx_pri + rest[:max_vars - len(idx_pri)])
    return [gene_names[i] for i in sel], X[:, sel]


# ---------------------------------------------------------------------------
# DREAM4
# ---------------------------------------------------------------------------

D4_PARTS = ["multifactorial", "knockouts", "knockdowns", "timeseries", "wildtype"]


def read_d4_table(text):
    """DREAM4 の TSV を (ヘッダ, 配列) にする。時系列は Time 列を落とし空行で区切る。"""
    lines = [l for l in text.splitlines() if l.strip()]
    header = [h.strip().strip('"') for h in lines[0].split("\t")]
    skip = 1 if header and header[0].lower() == "time" else 0
    rows = []
    for l in lines[1:]:
        cells = l.split("\t")[skip:]
        if len(cells) != len(header) - skip:
            continue
        try:
            rows.append([float(c) for c in cells])
        except ValueError:
            continue
    return header[skip:], np.asarray(rows, dtype=np.float64)


def prepare_dream4(zip_path, sizes, parts, out_dir, truth_dir, bins, method):
    if not os.path.isfile(zip_path):
        die(f"{zip_path} がありません (先に ./00download.sh)")
    z = zipfile.ZipFile(zip_path)
    names = z.namelist()
    done = []
    for size in sizes:
        for k in range(1, 6):
            net = f"insilico_size{size}_{k}"
            train = [n for n in names
                     if f"DREAM4 training data/{net}/" in n and n.endswith(".tsv")]
            if not train:
                log(f"  {net}: 訓練データが見つからないのでスキップ")
                continue
            gene_names, blocks = None, []
            for part in parts:
                hit = [n for n in train if n.endswith(f"{net}_{part}.tsv")]
                if not hit:
                    continue
                h, X = read_d4_table(dec(z.read(hit[0])))
                if X.size == 0:
                    continue
                if gene_names is None:
                    gene_names = h
                elif h != gene_names:
                    die(f"{net}: {part} の列順が他と違います")
                blocks.append(X)
            if not blocks:
                log(f"  {net}: 使えるデータがないのでスキップ")
                continue
            X = np.vstack(blocks)
            gs = [n for n in names
                  if n.endswith(f"DREAM4 gold standards/{net}_goldstandard.tsv")]
            if not gs:
                gs = [n for n in names if n.endswith(f"{net}_goldstandard.tsv")]
            edges, pairs = [], []
            for line in dec(z.read(gs[0])).splitlines():
                p = line.split("\t")
                if len(p) < 3:
                    continue
                pairs.append((p[0], p[1]))
                if p[2].strip() == "1":
                    edges.append((p[0], p[1]))
            present = write_dataset(net, gene_names, X, bins, method, out_dir)
            write_truth(net, edges, pairs, present, truth_dir)
            done.append(net)
    return done


# ---------------------------------------------------------------------------
# DREAM5
# ---------------------------------------------------------------------------

D5_LABELS = {"1": "insilico", "2": "saureus", "3": "ecoli", "4": "scerevisiae"}


def prepare_dream5(zip_path, networks, out_dir, truth_dir, bins, method, max_vars):
    if not os.path.isfile(zip_path):
        die(f"{zip_path} がありません (先に ./00download.sh)")
    z = zipfile.ZipFile(zip_path)
    names = z.namelist()
    done = []
    for num in networks:
        tag = f"dream5_net{num}_{D5_LABELS.get(str(num), 'net')}"
        expr = [n for n in names if n.endswith(f"net{num}_expression_data.tsv")]
        gs = [n for n in names
              if f"Network{num}/gold standard/" in n and n.endswith((".tsv", ".txt"))]
        if not expr or not gs:
            log(f"  Network{num}: データまたは gold standard が無いのでスキップ")
            continue
        lines = dec(z.read(expr[0])).splitlines()
        gene_names = [h.strip().strip('"') for h in lines[0].split("\t")]
        X = np.asarray([[float(c) for c in l.split("\t")] for l in lines[1:]
                        if l.strip()], dtype=np.float64)
        tf_file = [n for n in names if n.endswith(f"net{num}_transcription_factors.tsv")]
        tfs = dec(z.read(tf_file[0])).split() if tf_file else []
        log(f"  Network{num}: {X.shape[0]} サンプル x {X.shape[1]} 遺伝子 "
            f"(TF {len(tfs)} 個)")
        gene_names, X = select_columns(gene_names, X, max_vars, priority=tfs)

        edges, pairs = [], []
        for line in dec(z.read(sorted(gs)[0])).splitlines():
            p = line.split("\t")
            if len(p) < 2:
                continue
            u, v = p[0].strip(), p[1].strip()
            label = p[2].strip() if len(p) > 2 else "1"
            pairs.append((u, v))
            if label == "1":
                edges.append((u, v))
        present = write_dataset(tag, gene_names, X, bins, method, out_dir)
        write_truth(tag, edges, pairs, present, truth_dir)
        done.append(tag)
    return done


# ---------------------------------------------------------------------------
# HPN-DREAM (認証が要るため手動配置されている場合のみ)
# ---------------------------------------------------------------------------

def prepare_hpn(src_dir, out_dir, truth_dir, bins, method, max_vars):
    """HPN-DREAM の in silico サブチャレンジ用データが置かれていれば処理する。

    期待する配置 (Synapse からダウンロードして展開したもの):
        <src_dir>/expression.tsv   行 = サンプル, 列 = 変数, 1 行目 = 変数名
        <src_dir>/true_edges.tsv   正解エッジ (u <TAB> v)
    Synapse は認証が必要で自動取得できないため、無ければ何もしない。
    """
    expr = os.path.join(src_dir, "expression.tsv")
    truth = os.path.join(src_dir, "true_edges.tsv")
    if not (os.path.isfile(expr) and os.path.isfile(truth)):
        log("  HPN-DREAM: データが置かれていないのでスキップ "
            f"({expr} と {truth} を配置すると処理します)")
        return []
    with open(expr, encoding="utf-8") as fp:
        gene_names = [h.strip().strip('"') for h in
                      fp.readline().rstrip("\n").split("\t")]
        X = np.asarray([[float(c) for c in l.split("\t")] for l in fp if l.strip()],
                       dtype=np.float64)
    gene_names, X = select_columns(gene_names, X, max_vars)
    edges = []
    with open(truth, encoding="utf-8") as fp:
        for line in fp:
            p = line.split()
            if len(p) >= 2:
                edges.append((p[0], p[1]))
    present = write_dataset("hpn_dream", gene_names, X, bins, method, out_dir)
    write_truth("hpn_dream", edges, None, present, truth_dir)
    return ["hpn_dream"]


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="DREAM データを fast_bn 入力に整える")
    ap.add_argument("--dataset", required=True, choices=["dream4", "dream5", "hpn"])
    ap.add_argument("--zip", default=None, help="dream4 / dream5 の zip")
    ap.add_argument("--src-dir", default=None, help="hpn の展開済みディレクトリ")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--bins", type=int, default=3)
    ap.add_argument("--method", choices=["quantile", "uniform"], default="quantile")
    ap.add_argument("--sizes", default="10 100", help="dream4: 10 / 100")
    ap.add_argument("--parts", default=" ".join(D4_PARTS),
                    help="dream4: 結合する実験の種類")
    ap.add_argument("--networks", default="1 3 4", help="dream5: 対象ネットワーク番号")
    ap.add_argument("--max-vars", type=int, default=300,
                    help="dream5 / hpn: 使用する変数の上限 (0 = 制限なし)。"
                         "TF を優先し、残りを分散上位で選ぶ")
    ap.add_argument("--out-list", default=None,
                    help="作成したネットワーク名を 1 行ずつ書き出すファイル")
    args = ap.parse_args()

    if args.dataset == "dream4":
        done = prepare_dream4(args.zip, args.sizes.split(), args.parts.split(),
                              args.out_dir, args.truth_dir, args.bins, args.method)
    elif args.dataset == "dream5":
        done = prepare_dream5(args.zip, args.networks.split(), args.out_dir,
                              args.truth_dir, args.bins, args.method, args.max_vars)
    else:
        done = prepare_hpn(args.src_dir, args.out_dir, args.truth_dir,
                           args.bins, args.method, args.max_vars)

    if args.out_list:
        with open(args.out_list, "w", encoding="utf-8") as fp:
            for n in done:
                fp.write(n + "\n")
    log(f"{args.dataset}: {len(done)} ネットワークを準備しました")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_edges.py
================
学習したネットワークを**既知の正解構造**と比較して精度を評価する。
ダミーデータ (make_dummy_expr.py の true_edges.tsv) や、文献由来の
既知パスウェイと突き合わせる用途に使う。

比較は**遺伝子名**で行う (ノード番号は前処理のフィルタで変わるため)。
学習結果は `edges.tsv` (インデックス) と `edges_named.tsv` (名前) の行対応から
名前付きエッジとして読み込む。

評価指標:
  * 有向エッジ (向きも一致):  precision / recall / F1
  * 無向エッジ (向きは無視):  precision / recall / F1  (骨格の一致度)
  * SHD 相当 (対称差)         : 誤検出 + 未検出 の本数
  * 正解に含まれる遺伝子のうち、解析対象として残ったものの割合

使い方:
  python3 compare_edges.py --true data/true_edges.tsv \\
      --edges out/edges.tsv --edges-named out/edges_named.tsv \\
      --input data/expr_disc.tsv --out out/eval_vs_true.tsv
"""

import argparse
import os
import sys


def log(*args):
    print("[compare]", *args, file=sys.stderr)


def load_true_edges(path):
    """真のエッジ (parent child [weight]) を読み込む。ヘッダは自動判定。"""
    edges = set()
    with open(path, encoding="utf-8") as fp:
        for i, line in enumerate(fp):
            a = line.rstrip("\n").split("\t")
            if len(a) < 2:
                continue
            if i == 0 and a[0].lower() in ("parent", "u", "from", "source"):
                continue
            edges.add((a[0].strip(), a[1].strip()))
    return edges


def load_learned_edges(edges_path, named_path):
    """edges.tsv <-> edges_named.tsv の行対応から名前付きエッジ集合を作る。"""
    edges = set()
    with open(edges_path, encoding="utf-8") as fe, \
            open(named_path, encoding="utf-8") as fn:
        for le, ln in zip(fe, fn):
            n = ln.rstrip("\n").split("\t")
            if len(n) >= 2 and len(le.split()) >= 2:
                edges.add((n[0].strip(), n[1].strip()))
    return edges


def base_name(name):
    """前処理で付与された "__<gene_id>" 接尾辞を除いた素の遺伝子名。"""
    return name.split("__")[0]


def prf(tp, n_pred, n_true):
    prec = tp / n_pred if n_pred else 0.0
    rec = tp / n_true if n_true else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return prec, rec, f1


def main():
    ap = argparse.ArgumentParser(description="学習構造を既知の正解構造と比較する")
    ap.add_argument("--true", required=True, dest="true_path",
                    help="正解エッジ (parent<TAB>child[<TAB>weight])")
    ap.add_argument("--edges", default="out/edges.tsv",
                    help="学習エッジ (インデックス表記)")
    ap.add_argument("--edges-named", default="out/edges_named.tsv",
                    help="学習エッジ (遺伝子名表記; --edges と行対応)")
    ap.add_argument("--input", default=None,
                    help="離散化済み入力 TSV (解析対象に残った遺伝子の判定に使用)")
    ap.add_argument("--out", default=None, help="評価結果 TSV の出力先")
    ap.add_argument("--out-edges", default=None,
                    help="エッジ単位の判定 (TP/FP/FN) を書き出す TSV")
    ap.add_argument("--restrict-to-analyzed", action="store_true",
                    help="正解エッジのうち、両端が解析対象に残っているものだけを"
                         "評価対象にする (--input が必要)")
    args = ap.parse_args()

    for pth in (args.true_path, args.edges, args.edges_named):
        if not os.path.exists(pth):
            sys.exit(f"[compare] エラー: ファイルがありません: {pth}")

    true_edges = load_true_edges(args.true_path)
    learned = {(base_name(u), base_name(v))
               for u, v in load_learned_edges(args.edges, args.edges_named)}
    log(f"正解エッジ {len(true_edges)} 件, 学習エッジ {len(learned)} 件")

    analyzed = None
    if args.input and os.path.exists(args.input):
        with open(args.input, encoding="utf-8") as fp:
            analyzed = {base_name(c.strip())
                        for c in fp.readline().rstrip("\n").split("\t")}
        true_nodes = {n for e in true_edges for n in e}
        kept = true_nodes & analyzed
        log(f"正解の遺伝子 {len(true_nodes)} 件のうち {len(kept)} 件 "
            f"({len(kept) / max(1, len(true_nodes)):.1%}) が解析対象に残存 "
            f"(解析対象は全 {len(analyzed)} 変数)")

    eval_true = true_edges
    if args.restrict_to_analyzed:
        if analyzed is None:
            sys.exit("[compare] エラー: --restrict-to-analyzed には --input が必要です")
        eval_true = {(u, v) for u, v in true_edges if u in analyzed and v in analyzed}
        log(f"評価対象の正解エッジ: {len(eval_true)} 件 (両端が解析対象に残るもの)")

    # 有向での一致
    tp_d = learned & eval_true
    fp_d = learned - eval_true
    fn_d = eval_true - learned
    p_d, r_d, f_d = prf(len(tp_d), len(learned), len(eval_true))

    # 無向 (骨格) での一致
    und_true = {frozenset(e) for e in eval_true}
    und_learn = {frozenset(e) for e in learned}
    tp_u = und_learn & und_true
    p_u, r_u, f_u = prf(len(tp_u), len(und_learn), len(und_true))

    # 向きだけ誤り (骨格は一致しているが向きが逆)
    rev = {(u, v) for (u, v) in learned if (v, u) in eval_true and (u, v) not in eval_true}
    shd = len(fp_d) + len(fn_d)

    rows = [
        ("n_true_edges", len(eval_true)),
        ("n_learned_edges", len(learned)),
        ("directed_TP", len(tp_d)),
        ("directed_FP", len(fp_d)),
        ("directed_FN", len(fn_d)),
        ("directed_precision", f"{p_d:.4f}"),
        ("directed_recall", f"{r_d:.4f}"),
        ("directed_f1", f"{f_d:.4f}"),
        ("undirected_TP", len(tp_u)),
        ("undirected_precision", f"{p_u:.4f}"),
        ("undirected_recall", f"{r_u:.4f}"),
        ("undirected_f1", f"{f_u:.4f}"),
        ("reversed_edges", len(rev)),
        ("symmetric_difference (SHD 相当)", shd),
    ]
    width = max(len(k) for k, _ in rows)
    print("=" * (width + 14))
    print(" 学習構造 vs 正解構造")
    print("=" * (width + 14))
    for k, v in rows:
        print(f" {k:<{width}} : {v}")
    print("=" * (width + 14))
    print(" 注: 少サンプル・離散化・スコア基準の違いから、有向一致は一般に"
          " 骨格一致より低くなります。")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fp:
            fp.write("metric\tvalue\n")
            for k, v in rows:
                fp.write(f"{k}\t{v}\n")
        log(f"出力: {args.out}")

    if args.out_edges:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_edges)) or ".",
                    exist_ok=True)
        with open(args.out_edges, "w", encoding="utf-8") as fp:
            fp.write("u\tv\tstatus\n")
            for u, v in sorted(tp_d):
                fp.write(f"{u}\t{v}\tTP\n")
            for u, v in sorted(fp_d):
                st = "FP_reversed" if (u, v) in rev else "FP"
                fp.write(f"{u}\t{v}\t{st}\n")
            for u, v in sorted(fn_d):
                fp.write(f"{u}\t{v}\tFN\n")
        log(f"出力: {args.out_edges}")


if __name__ == "__main__":
    main()

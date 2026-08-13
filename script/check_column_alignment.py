#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_column_alignment.py
=========================
fast_bn のノード番号は**入力データの列位置**に対応する。したがって
`--score-dataset` に別のデータ (群別サブセットや検証データ) を渡して既存の
ネットワークを評価するときは、そのデータの列 (遺伝子) 順が、対象ネットワークを
学習したときの入力と完全に一致していなければならない。

このスクリプトは、
  * 評価に使う入力 TSV のヘッダ (列名の並び) と
  * 対象ネットワークの「インデックス表記エッジ」と「名前表記エッジ」の行対応
    (edges.tsv <-> edges_named.tsv) から復元した idx -> 名前
を突き合わせ、不一致があれば異常終了する (終了コード 2)。

使い方:
  python3 check_column_alignment.py <input.tsv> <edges.tsv> <edges_named.tsv>
"""

import sys


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: check_column_alignment.py <input.tsv> <edges.tsv> "
                 "<edges_named.tsv>")
    inp, edges_path, named_path = sys.argv[1:4]

    with open(inp, encoding="utf-8") as fp:
        header = fp.readline().rstrip("\n").split("\t")

    idx2name = {}
    with open(edges_path, encoding="utf-8") as fe, \
            open(named_path, encoding="utf-8") as fn:
        for le, ln in zip(fe, fn):
            e = le.split()
            n = ln.rstrip("\n").split("\t")
            if len(e) >= 2 and len(n) >= 2:
                try:
                    idx2name[int(e[0])] = n[0]
                    idx2name[int(e[1])] = n[1]
                except ValueError:
                    continue  # ヘッダ行等はスキップ

    if not idx2name:
        print("[align] 警告: 参照エッジが空です (検証をスキップ)", file=sys.stderr)
        return

    bad = [(i, header[i] if i < len(header) else "(範囲外)", nm)
           for i, nm in sorted(idx2name.items())
           if i >= len(header) or header[i] != nm]
    if bad:
        print(f"[align] 不一致 {len(bad)} 件 (例: {bad[:3]})", file=sys.stderr)
        print(f"[align] 入力の列数={len(header)}, 参照ノード数={len(idx2name)}",
              file=sys.stderr)
        print("[align] エラー: 入力の列順が対象ネットワークの学習入力と一致しません。",
              file=sys.stderr)
        print("        対象網を学習したときの入力 (列=遺伝子の順序が同じファイル) を",
              file=sys.stderr)
        print("        指定してください。前処理のパラメータを変えて学習結果を作り直した",
              file=sys.stderr)
        print("        場合は、その学習に使った入力をそのまま使います。", file=sys.stderr)
        sys.exit(2)
    print(f"[align] OK: {len(idx2name)} ノードが一致 (列数={len(header)})")


if __name__ == "__main__":
    main()

#!/usr/bin/env bash
# =============================================================================
# 01preprocess.sh — ステップ 1: 前処理
#   生カウント -> CPM 正規化 -> log2 変換 -> フィルタ -> 3 段階離散化 ->
#   fast_bn 入力 TSV。サンプルはメタデータの群順に並べ替えられる。
#
# 出力:
#   data/expr_disc.tsv : fast_bn 入力 (行=サンプル, 列=遺伝子, 値=0..2)
#   data/var_map.tsv   : 列インデックス <-> gene id / name / 統計量
#   data/samples.tsv   : 行番号 <-> サンプル ID / 群ラベル (群別解析に使う)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

"${BN_SCRIPTS}/preprocess.sh" "$@"

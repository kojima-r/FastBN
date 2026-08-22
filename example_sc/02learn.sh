#!/usr/bin/env bash
# =============================================================================
# 02learn.sh — ステップ 2: ネットワーク構造学習
#   fast_bn の Hill-Climb + Tabu サーチで DAG を推定する。
#   スコア・最大親数などは config.sh (SCORE / ESS / MAX_PARENTS / ITERS ...)。
#
# 出力 (${OUTDIR} = ${RUNDIR}/out):
#   edges.tsv       : 学習された DAG (ノード=列インデックス)
#   edges_named.tsv : 同じエッジを遺伝子名で表記 (edges.tsv と行対応)
#   all_counts.tsv  : カウント表 (CPT 推定・エッジ重要度評価に使う)
#   log_learn.txt   : fast_bn のログ
#
# 使い方:
#   ./02learn.sh
#   SCORE=bic MAX_PARENTS=2 ITERS=20000 ./02learn.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

"${BN_SCRIPTS}/learn_structure.sh" "$@"

echo "[02learn] 上位のハブ遺伝子 (出次数):"
cut -f1 "${OUTDIR}/edges_named.tsv" | sort | uniq -c | sort -rn | head -5 \
  | sed 's/^/   /' || true

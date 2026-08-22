#!/usr/bin/env bash
# =============================================================================
# run_all.sh — ステップ 0〜6 を順番に実行する
#   ./run_all.sh                              # 既定 (DREAM4 + DREAM5)
#   DATASETS=dream4 ./run_all.sh              # DREAM4 だけ
#   D5_MAX_VARS=0 ./run_all.sh                # DREAM5 を全遺伝子で (非常に重い)
#   SKIP_DOWNLOAD=1 ./run_all.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
start=$(date +%s)

echo "=================================================================="
echo " DREAM チャレンジ ベンチマーク"
echo "   データセット : ${DATASETS}"
echo "   DREAM4 サイズ: ${D4_SIZES}"
echo "   DREAM5 網    : ${D5_NETWORKS} (変数上限 ${D5_MAX_VARS})"
echo "   離散化       : ${BINS} 段階 (${DISC_METHOD})"
echo "   スコア       : ${SCORES}"
echo "   出力         : ${RUNDIR}/"
echo "=================================================================="

[[ "${SKIP_DOWNLOAD:-0}" == "0" ]] && ./00download.sh
./01prepare.sh
./02learn.sh
./03evaluate.sh
./04summarize.sh
./05compare.sh
./06report.sh

echo "=================================================================="
echo " 完了 ($(( $(date +%s) - start )) 秒)"
echo "   ${REPORT_HTML}  : HTML レポート"
echo "   ${BENCHMARK}    : 実行ごとの生の指標"
echo "=================================================================="

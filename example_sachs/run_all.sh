#!/usr/bin/env bash
# =============================================================================
# run_all.sh — ステップ 0〜6 を順番に実行する
#   ./run_all.sh                     # 既定 (18 実行, 数分)
#   PRESETS=obs BINS_LIST=3 ./run_all.sh
#   SKIP_DOWNLOAD=1 ./run_all.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
start=$(date +%s)

echo "=================================================================="
echo " Sachs タンパク質シグナル伝達データ ベンチマーク"
echo "   条件セット : ${PRESETS}"
echo "   離散化     : ${BINS_LIST} 段階"
echo "   スコア     : ${SCORES}"
echo "   出力       : ${RUNDIR}/"
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

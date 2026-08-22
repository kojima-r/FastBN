#!/usr/bin/env bash
# =============================================================================
# 04summarize.sh — ステップ 4: 集約
#   ${BENCHMARK} を (dataset, network, score) ごとに集約する。
#   データセット単位の要約も出す。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
[[ -s "${BENCHMARK}" ]] || { echo "[04summarize] エラー: ${BENCHMARK} がありません" >&2; exit 1; }

python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" --group-by "${SUMMARY_GROUP_BY}" \
  --metrics "${SUMMARY_METRICS}" --plot-metrics "${SUMMARY_PLOT_METRICS}" \
  --out "${SUMMARY_TSV}" --markdown "${SUMMARY_MD}" \
  --plot "${SUMMARY_PNG}" --plot-x n_vars --plot-facet dataset --plot-series score "$@"

python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" --group-by dataset --metrics "${SUMMARY_METRICS}" \
  --out "${SUMMARY_OVERALL_TSV}"

echo
echo "[04summarize] データセットごとの平均:"
python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" --group-by dataset --metrics "${SUMMARY_METRICS}"
echo
echo "[04summarize] 出力: ${SUMMARY_TSV}, ${SUMMARY_MD}, ${SUMMARY_PNG}"

#!/usr/bin/env bash
# =============================================================================
# 04summarize.sh — ステップ 4: 集約
#   ${BENCHMARK} を (preset, bins, score) ごとに集約し、平均 ± 標準偏差の表と
#   グラフを出力する。ここでは反復が無いので標準偏差は 0 になり、表は実質
#   各条件の値そのものになる。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
[[ -s "${BENCHMARK}" ]] || { echo "[04summarize] エラー: ${BENCHMARK} がありません。先に ./03evaluate.sh" >&2; exit 1; }

python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" --group-by "${SUMMARY_GROUP_BY}" \
  --metrics "${SUMMARY_METRICS}" --plot-metrics "${SUMMARY_PLOT_METRICS}" \
  --out "${SUMMARY_TSV}" --markdown "${SUMMARY_MD}" \
  --plot "${SUMMARY_PNG}" --plot-x bins --plot-facet preset --plot-series score "$@"

python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" --group-by preset --metrics "${SUMMARY_METRICS}" \
  --out "${SUMMARY_OVERALL_TSV}"

echo
echo "[04summarize] 条件セットごとの平均:"
python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" --group-by preset --metrics "${SUMMARY_METRICS}"
echo
echo "[04summarize] 出力: ${SUMMARY_TSV}, ${SUMMARY_MD}, ${SUMMARY_PNG}"

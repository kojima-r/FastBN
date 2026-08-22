#!/usr/bin/env bash
# =============================================================================
# 04summarize.sh — ステップ 4: ベンチマーク結果の集約
#   ${BENCHMARK} を (network, n, score) ごとに集約し、繰り返しにわたる
#   平均 ± 標準偏差の表とグラフを出力する。
#
# 出力:
#   ${SUMMARY_TSV}         : 集約結果 (機械処理用)
#   ${SUMMARY_MD}          : 同じ内容の Markdown 表
#   ${SUMMARY_PNG}         : 指標 x ネットワークの折れ線グラフ (横軸 = サンプル数)
#   ${SUMMARY_OVERALL_TSV} : サンプル数だけで集約したもの (全体傾向)
#
# 表には Precision / Recall / F1 を directed・skeleton の両方について出す
# (SUMMARY_METRICS)。グラフは見やすさのため指標を絞る (SUMMARY_PLOT_METRICS)。
#
# 集約キーや指標は config.sh で変更できる。
#   SUMMARY_GROUP_BY=network,score ./04summarize.sh   # サンプル数をまとめる
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

[[ -s "${BENCHMARK}" ]] || { echo "[04summarize] エラー: ${BENCHMARK} がありません。先に ./03evaluate.sh" >&2; exit 1; }

python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" \
  --group-by "${SUMMARY_GROUP_BY}" \
  --metrics "${SUMMARY_METRICS}" \
  --plot-metrics "${SUMMARY_PLOT_METRICS}" \
  --out "${SUMMARY_TSV}" \
  --markdown "${SUMMARY_MD}" \
  --plot "${SUMMARY_PNG}" \
  --plot-x n --plot-facet network --plot-series score \
  "$@"

# サンプル数だけで集約した全体傾向 (レポートの冒頭に出す)
python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" --group-by n --metrics "${SUMMARY_METRICS}" \
  --out "${SUMMARY_OVERALL_TSV}"

echo
echo "[04summarize] サンプル数ごとの平均 (全ネットワーク・全スコア):"
python3 "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input "${BENCHMARK}" --group-by n --metrics "${SUMMARY_METRICS}"
echo
echo "[04summarize] 出力: ${SUMMARY_TSV}, ${SUMMARY_MD}, ${SUMMARY_PNG}"

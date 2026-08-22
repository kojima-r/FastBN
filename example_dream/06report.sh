#!/usr/bin/env bash
# =============================================================================
# 06report.sh — ステップ 6: HTML レポート
#     ./06report.sh --embed   # 画像を埋め込んで 1 ファイルにする
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
[[ -s "${BENCHMARK}" ]] || { echo "[06report] エラー: ${BENCHMARK} がありません" >&2; exit 1; }

info="${RUNDIR}/network_info.tsv"
{
  printf 'dataset\tnetwork\tsamples\tvariables\ttrue_edges\tevaluable_pairs\n'
  while IFS=$'\t' read -r ds net; do
    d="${DATADIR}/${net}.tsv"; e="${TRUTHDIR}/${net}_edges.tsv"; p="${TRUTHDIR}/${net}_pairs.tsv"
    [[ -s "$d" ]] || continue
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${ds}" "${net}" \
      "$(( $(wc -l < "$d") - 1 ))" "$(head -1 "$d" | awk -F'\t' '{print NF}')" \
      "$(wc -l < "$e" 2>/dev/null || echo 0)" "$(wc -l < "$p" 2>/dev/null || echo '-')"
  done < "${NETLIST}"
} > "${info}"

python3 "${BN_SCRIPTS}/make_benchmark_report.py" \
  --benchmark "${BENCHMARK}" --summary "${SUMMARY_TSV}" \
  --summary-overall "${SUMMARY_OVERALL_TSV}" --networks "${info}" \
  --plot "${SUMMARY_PNG}" --compare-dir "${FIGDIR}" --out "${REPORT_HTML}" \
  --title "DREAM チャレンジ 構造学習ベンチマーク" \
  --subtitle "データセット: ${DATASETS} / 離散化: ${BINS} 段階 / スコア: ${SCORES}" "$@"

echo "[06report] ブラウザで開いて確認してください: $(pwd)/${REPORT_HTML#./}"

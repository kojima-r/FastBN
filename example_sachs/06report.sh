#!/usr/bin/env bash
# =============================================================================
# 06report.sh — ステップ 6: HTML レポート
#   集約表・グラフ・比較図を 1 つの HTML にまとめる。
#     ./06report.sh --embed   # 画像を埋め込んで 1 ファイルにする
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
[[ -s "${BENCHMARK}" ]] || { echo "[06report] エラー: ${BENCHMARK} がありません" >&2; exit 1; }

info="${RUNDIR}/datasets.tsv"
{
  printf 'dataset\tcells\tproteins\tbins\n'
  for preset in ${PRESETS}; do
    for bins in ${BINS_LIST}; do
      f="${DATADIR}/sachs_${preset}_b${bins}.tsv"
      [[ -s "$f" ]] || continue
      printf '%s\t%s\t%s\t%s\n' "${preset}" "$(( $(wc -l < "$f") - 1 ))" \
        "$(head -1 "$f" | awk -F'\t' '{print NF}')" "${bins}"
    done
  done
} > "${info}"

python3 "${BN_SCRIPTS}/make_benchmark_report.py" \
  --benchmark "${BENCHMARK}" --summary "${SUMMARY_TSV}" \
  --summary-overall "${SUMMARY_OVERALL_TSV}" --networks "${info}" \
  --plot "${SUMMARY_PNG}" --compare-dir "${FIGDIR}" --out "${REPORT_HTML}" \
  --title "Sachs タンパク質シグナル伝達データ 構造学習ベンチマーク" \
  --subtitle "条件セット: ${PRESETS} / 離散化: ${BINS_LIST} 段階 / スコア: ${SCORES}" "$@"

echo "[06report] ブラウザで開いて確認してください: $(pwd)/${REPORT_HTML#./}"

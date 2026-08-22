#!/usr/bin/env bash
# =============================================================================
# 05compare.sh — ステップ 5: 正解経路と学習ネットワークの比較図
#   ../script/plot_dag_comparison.py で左に正解・右に学習結果を同じ配置で描く。
#   緑 = 向きまで一致 / 橙 = 向きが逆 / 赤 = 余分 / 灰破線 = 見落とし
#
# 出力: ${FIGDIR}/<preset>/<preset>_b<bins>_<score>.png
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
mkdir -p "${FIGDIR}"

viz_presets="${VIZ_PRESETS:-${PRESETS}}"
viz_bins="${VIZ_BINS:-${BINS_LIST}}"
viz_scores="${VIZ_SCORES:-${SCORES}}"
echo "[05compare] 対象: preset [${viz_presets}] / bins [${viz_bins}] / score [${viz_scores}]"

count=0
for preset in ${viz_presets}; do
  mkdir -p "${FIGDIR}/${preset}"
  for bins in ${viz_bins}; do
    data="${DATADIR}/sachs_${preset}_b${bins}.tsv"
    for score in ${viz_scores}; do
      run="${preset}_b${bins}_${score}"
      edges="${OUTROOT}/${run}/edges.tsv"
      [[ -f "${edges}" ]] || { echo "[05compare] 警告: ${run} をスキップ" >&2; continue; }
      python3 "${BN_SCRIPTS}/plot_dag_comparison.py" \
        --true-edges "${TRUE_EDGES}" --pred-edges "${edges}" --input "${data}" \
        --title "Sachs / ${preset} / ${bins} bins / ${score}" \
        --out "${FIGDIR}/${preset}/${run}.png"
      count=$(( count + 1 ))
    done
  done
done
echo "[05compare] 完了: ${count} 枚 -> ${FIGDIR}/"

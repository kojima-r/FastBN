#!/usr/bin/env bash
# =============================================================================
# 05compare.sh — ステップ 5: 正解 vs 学習ネットワークの比較図
#   ../script/plot_dag_comparison.py で左に正解・右に学習結果を同じ配置で描く。
#   緑 = 向きまで一致 / 橙 = 向きが逆 / 赤 = 余分 / 灰破線 = 見落とし
#
#   DREAM5 のように数千ノードある網は全体を描いても判読できないので、
#   VIZ_MAX_NODES を超える場合は**次数上位のノードだけを取り出した部分グラフ**を
#   描く (図の副題に明記される)。副題の指標はネットワーク全体に対する値のまま。
#
# 出力: ${FIGDIR}/<dataset>/<network>_<score>.png
#
#   VIZ_MAX_NODES=200 ./05compare.sh   # もう少し大きく描く
#   VIZ_MAX_NODES=0   ./05compare.sh   # 部分グラフにせず全ノード描く (巨大網は非推奨)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
[[ -s "${NETLIST}" ]] || { echo "[05compare] エラー: ${NETLIST} がありません" >&2; exit 1; }
mkdir -p "${FIGDIR}"

viz_scores="${VIZ_SCORES:-${SCORES}}"
count=0; subgraph=0
# DATASETS に含まれないデータセットは飛ばす (NETLIST は 01prepare 時点の全件)
in_datasets() { local d; for d in ${DATASETS}; do [[ "$d" == "$1" ]] && return 0; done; return 1; }
while IFS=$'\t' read -r ds net; do
  in_datasets "${ds}" || continue
  data="${DATADIR}/${net}.tsv"
  [[ -s "${data}" ]] || continue
  nv=$(head -1 "${data}" | awk -F'\t' '{print NF}')
  [[ "${VIZ_MAX_NODES}" -gt 0 && "${nv}" -gt "${VIZ_MAX_NODES}" ]] && subgraph=$(( subgraph + 1 ))
  mkdir -p "${FIGDIR}/${ds}"
  for score in ${viz_scores}; do
    pred="${OUTROOT}/${net}_${score}/edges.tsv"
    [[ -f "${pred}" ]] || continue
    python3 "${BN_SCRIPTS}/plot_dag_comparison.py" \
      --true-edges "${TRUTHDIR}/${net}_edges.tsv" --pred-edges "${pred}" \
      --input "${data}" --max-nodes "${VIZ_MAX_NODES}" \
      --title "${net} / ${score} (${nv} vars)" \
      --out "${FIGDIR}/${ds}/${net}_${score}.png"
    count=$(( count + 1 ))
  done
done < "${NETLIST}"
echo "[05compare] 完了: ${count} 枚 -> ${FIGDIR}/"
[[ "${subgraph}" -gt 0 ]] && \
  echo "[05compare] (うち ${subgraph} ネットワークは ${VIZ_MAX_NODES} ノードのハブ部分グラフとして描画)"
exit 0

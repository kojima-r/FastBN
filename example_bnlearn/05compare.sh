#!/usr/bin/env bash
# =============================================================================
# 05compare.sh — ステップ 5: 正解ネットワークと学習ネットワークの比較図
#   ../script/plot_dag_comparison.py で、左に正解 DAG・右に学習 DAG を
#   **同じノード配置**で並べた図を作る。学習側のエッジは判定で色分けする:
#     緑 = 向きまで一致 / 橙 = 向きが逆 / 赤 = 余分 / 灰破線 = 見落とし
#
#   全 300 実行ぶん描くと多すぎるので、既定では
#     「各ネットワーク x 各スコア x VIZ_N のサンプル数 x 反復 VIZ_REP」
#   に絞る (既定 = 5 ネットワーク x 3 スコア x 1 = 15 枚)。
#
# 出力: ${FIGDIR}/<network>/<network>_n<N>_r<R>_<score>.png
#
# 対象の変え方:
#   VIZ_N="100 5000" ./05compare.sh    # 少数データと多数データを見比べる
#   VIZ_N=all ./05compare.sh           # SAMPLE_SIZES 全部 (枚数に注意)
#   VIZ_SCORES=bic VIZ_REP=2 ./05compare.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${FIGDIR}"

# VIZ_N=all なら SAMPLE_SIZES 全部、未指定なら最大のサンプル数
if [[ "${VIZ_N:-}" == "all" ]]; then
  viz_sizes="${SAMPLE_SIZES}"
elif [[ -n "${VIZ_N:-}" ]]; then
  viz_sizes="${VIZ_N}"
else
  viz_sizes=$(echo ${SAMPLE_SIZES} | tr ' ' '\n' | sort -n | tail -1)
fi
viz_scores="${VIZ_SCORES:-${SCORES}}"
viz_rep="${VIZ_REP:-1}"

echo "[05compare] 対象: サンプル数 [${viz_sizes}] / スコア [${viz_scores}] / 反復 ${viz_rep}"

count=0
for net in ${NETWORKS}; do
  bif="${NETDIR}/${net}.bif"
  [[ -s "${bif}" ]] || { echo "[05compare] エラー: ${bif} がありません" >&2; exit 1; }
  mkdir -p "${FIGDIR}/${net}"
  for n in ${viz_sizes}; do
    data="${DATADIR}/${net}_n${n}_r${viz_rep}.tsv"
    for score in ${viz_scores}; do
      run="${net}_n${n}_r${viz_rep}_${score}"
      edges="${OUTROOT}/${run}/edges.tsv"
      if [[ ! -f "${edges}" ]]; then
        echo "[05compare] 警告: ${edges} が無いのでスキップ (${run})" >&2
        continue
      fi
      python3 "${BN_SCRIPTS}/plot_dag_comparison.py" \
        --true-bif "${bif}" \
        --pred-edges "${edges}" \
        --input "${data}" \
        --alpha "${KL_ALPHA}" \
        --max-states "${MAX_STATES}" \
        --title "${net} / ${score} / n=${n} (rep ${viz_rep})" \
        --out "${FIGDIR}/${net}/${run}.png"
      count=$(( count + 1 ))
    done
  done
done

echo "[05compare] 完了: ${count} 枚 -> ${FIGDIR}/"

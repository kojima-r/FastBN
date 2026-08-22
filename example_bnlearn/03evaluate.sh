#!/usr/bin/env bash
# =============================================================================
# 03evaluate.sh — ステップ 3: 正解ネットワークとの比較 (5 指標)
#   全実行について ../script/evaluate_structure.py を呼び、結果を 1 行ずつ
#   ${BENCHMARK} に追記する。
#
# 計算される指標:
#   shd                              Structural Hamming Distance
#   precision/recall/f1_directed     向きまで含めたエッジの一致
#   precision/recall/f1_skeleton     向きを無視した骨格の一致
#   sid / sid_normalized             Structural Intervention Distance
#   kl_divergence                    KL(P_true || P_learned)  (nat, 厳密計算)
#
# 出力:
#   ${EVALDIR}/<run>.tsv  : 実行ごとの全指標
#   ${BENCHMARK}          : 全実行を 1 行ずつまとめたもの (04report.sh の入力)
#
# KL は「学習した構造 + 学習に使ったデータから推定した CPT」の同時分布と
# 正解ネットワークの同時分布を全状態列挙で比較する。CPT の推定には
# KL_ALPHA の Dirichlet 平滑化を使う (0 だと発散しうる)。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${EVALDIR}" "$(dirname "${BENCHMARK}")"
rm -f "${BENCHMARK}"        # 追記なので毎回作り直す

count=0
for net in ${NETWORKS}; do
  bif="${NETDIR}/${net}.bif"
  [[ -s "${bif}" ]] || { echo "[03evaluate] エラー: ${bif} がありません" >&2; exit 1; }
  for n in ${SAMPLE_SIZES}; do
    for rep in $(seq 1 "${REPLICATES}"); do
      data="${DATADIR}/${net}_n${n}_r${rep}.tsv"
      for score in ${SCORES}; do
        run="${net}_n${n}_r${rep}_${score}"
        edges="${OUTROOT}/${run}/edges.tsv"
        [[ -f "${edges}" ]] || { echo "[03evaluate] エラー: ${edges} がありません。先に ./02learn.sh" >&2; exit 1; }
        python3 "${BN_SCRIPTS}/evaluate_structure.py" \
          --true-bif "${bif}" \
          --pred-edges "${edges}" \
          --input "${data}" \
          --alpha "${KL_ALPHA}" \
          --max-states "${MAX_STATES}" \
          --out "${EVALDIR}/${run}.tsv" \
          --append "${BENCHMARK}" \
          --extra "network=${net}" --extra "n=${n}" \
          --extra "rep=${rep}" --extra "score=${score}" \
          --quiet
        count=$(( count + 1 ))
        printf "\r[03evaluate] %d 実行を評価しました          " "${count}"
      done
    done
  done
done

echo
echo "[03evaluate] 完了: ${count} 行 -> ${BENCHMARK}"

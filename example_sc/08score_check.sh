#!/usr/bin/env bash
# =============================================================================
# 08score_check.sh — ステップ 8: モデルスコア (対数尤度) の確認
#   学習した構造とカウント表を固定し、データセットに対する対数尤度を計算する
#   (fast_bn の Scoring Mode)。ブートストラップ統合の前後で当てはまりが
#   どう変わったかを比較できる。
#
#   example_bulk の 08evaluate.sh (真の DAG との比較) に相当するステップだが、
#   実データには正解構造が無いため、代わりに尤度で評価する。
#
#   既定では**学習に使ったのと同じデータ**を評価データにしている (当てはまりの
#   確認)。汎化性能を見たい場合は EVAL_INPUT に別のデータを指定する:
#     EVAL_INPUT=./run_bbknn_bin100/groups/expr_g1_Aorta.tsv ./08score_check.sh
#   ※ 列 (遺伝子) の順序が学習時の入力と一致している必要がある。
#
# 出力:
#   ${OUTDIR}/score_hc.tsv : 学習網の尤度
#   ${OUTDIR}/score_bs.tsv : コンセンサス網の尤度 (04 実行後)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

FASTBN_BIN="${FASTBN_BIN:-../fast_bn}"
eval_input="${EVAL_INPUT:-${INPUT}}"
[[ -f "${eval_input}" ]] || { echo "[08score_check] エラー: ${eval_input} がありません" >&2; exit 1; }

score_one() {  # <ラベル> <edges> <counts> <出力>
  local label="$1" edges="$2" counts="$3" out="$4"
  echo "=================================================================="
  echo " ${label}"
  echo "=================================================================="
  "${FASTBN_BIN}" --score "${SCORE_IMP}" \
    --score-dataset "${eval_input}" \
    --init   "${edges}" \
    --counts "${counts}" \
    --alpha  "${ALPHA}" \
    | tee "${out}"
}

echo "[08score_check] 評価データ: ${eval_input}"
score_one "(a) Hill-Climb 学習網" \
  "${OUTDIR}/edges.tsv" "${OUTDIR}/all_counts.tsv" "${OUTDIR}/score_hc.tsv"

if [[ -f "${OUTDIR}/integ_edges2.tsv" && -f "${OUTDIR}/integ_all_counts.tsv" ]]; then
  score_one "(b) ブートストラップ・コンセンサス網" \
    "${OUTDIR}/integ_edges2.tsv" "${OUTDIR}/integ_all_counts.tsv" \
    "${OUTDIR}/score_bs.tsv"

  echo "=================================================================="
  echo " まとめ"
  echo "=================================================================="
  echo "   学習網          : エッジ数 $(wc -l < "${OUTDIR}/edges.tsv"), logL/サンプル $(awk -F'\t' '$1=="log_likelihood_per_sample"{print $2}' "${OUTDIR}/score_hc.tsv")"
  echo "   コンセンサス網  : エッジ数 $(wc -l < "${OUTDIR}/integ_edges2.tsv"), logL/サンプル $(awk -F'\t' '$1=="log_likelihood_per_sample"{print $2}' "${OUTDIR}/score_bs.tsv")"
  echo
  echo "   コンセンサス網は各リサンプルで繰り返し現れたエッジだけを残したもの。"
  echo "   閾値 (THRESHOLD_PROB / THRESHOLD_COUNT) を上げるとエッジが減って"
  echo "   学習データへの当てはまり (logL) は下がり、下げると逆に学習網より"
  echo "   エッジが増えることもある。logL の差が小さいまま構造が単純になるなら、"
  echo "   落としたエッジは過学習だった可能性が高い。"
else
  echo "[08score_check] コンセンサス網が無いため (b) はスキップ"
  echo "                (先に ./04bootstrap.sh を実行すると比較できます)"
fi

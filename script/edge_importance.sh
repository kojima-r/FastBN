#!/usr/bin/env bash
# =============================================================================
# edge_importance.sh
#   学習済みネットワークの各エッジの重要度 (エッジを除いたときのスコア変化) を
#   fast_bn の Edge Importance Mode で計算する (汎用版)。
#
# 入力 (環境変数で上書き可能):
#   INPUT    : 評価データセット (既定 ./data/expr_disc.tsv)
#   INIT     : 対象ネットワークの構造     (既定 ${OUTDIR}/edges.tsv)
#   COUNTS   : 対象ネットワークのカウント表 (既定 ${OUTDIR}/all_counts.tsv)
#   OUT_IMP  : 出力パス                   (既定 ${OUTDIR}/edge_importance.tsv)
#
# 出力: OUT_IMP
#   列: u v ΔlogL ΔBIC ΔK2 ΔBDeu meanΔlogL_per_sample stdΔlogL_per_sample
#
# 使い方:
#   ../script/edge_importance.sh                       # 学習網
#   INIT=out/integ_edges.tsv COUNTS=out/integ_all_counts.tsv \
#     OUT_IMP=out/integ_edge_importance.tsv ../script/edge_importance.sh
#
# 注意: INPUT の列 (遺伝子) の順序は、対象ネットワークを学習したときの入力と
#       一致していなければならない (ノード番号 = 列位置)。
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="importance"; bn_tag="importance"

init="${INIT:-${OUTDIR}/edges.tsv}"
counts="${COUNTS:-${OUTDIR}/all_counts.tsv}"
out_imp="${OUT_IMP:-${OUTDIR}/edge_importance.tsv}"

score="${SCORE_IMP:-${SCORE:-bic}}"   # 重要度評価に用いるスコア
alpha="${ALPHA:-1.0}"                 # スムージング係数
ess="${ESS:-10.0}"                    # BDeu の等価サンプルサイズ

require_bin
require_file "${INPUT}" "${init}" "${counts}"
mkdir -p "$(dirname "${out_imp}")"

log "init=${init} counts=${counts} score-dataset=${INPUT} score=${score}"

"${FASTBN_BIN}" --score "${score}" \
  --edge-importance \
  --score-dataset "${INPUT}" \
  --init   "${init}" \
  --counts "${counts}" \
  --alpha "${alpha}" \
  --ess "${ess}" \
  --save-edge-importance "${out_imp}"

log "完了: ${out_imp} ($(( $(wc -l < "${out_imp}") - 1 )) エッジ)"

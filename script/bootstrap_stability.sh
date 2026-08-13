#!/usr/bin/env bash
# =============================================================================
# bootstrap_stability.sh
#   ブートストラップ・リサンプリングによるエッジ安定性解析 (汎用版)。
#   サンプルを復元抽出して構造学習を多数回繰り返し、各エッジの出現頻度
#   (ブートストラップ確率) からコンセンサスネットワークを構築し、その網の
#   エッジ重要度まで計算する。
#
#   流れ (FastBN/example の 02→03→05 に対応):
#     A) fast_bn --bootstrap  : 複数シードを並列にリサンプリング学習
#                               -> ${BSDIR}/edges_seed####.tsv (u v count prob)
#     B) compute_bs_prob.py   : 全シードを統合しコンセンサスエッジを抽出
#                               -> ${OUTDIR}/integ_edges.tsv, integ_edges_score.tsv
#     C) fast_bn --iters 0    : コンセンサス構造のカウント表を再計算
#                               -> ${OUTDIR}/integ_edges2.tsv, integ_edges_named.tsv,
#                                  integ_all_counts.tsv
#     D) fast_bn --edge-importance : コンセンサス網のエッジ重要度
#                               -> ${OUTDIR}/integ_edge_importance.tsv
#
# 前提: INPUT (既定 ./data/expr_disc.tsv) が作成済み。
#       ${OUTDIR}/edges.tsv があれば初期構造として利用する (warm start)。
#
# 使い方:
#   ../script/bootstrap_stability.sh
#   BOOTSTRAP=20 SEEDS=10 MAX_JOBS=4 ITERS=2000 ../script/bootstrap_stability.sh
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="bootstrap"; bn_tag="bootstrap"

require_bin
require_file "${INPUT}"
[[ -f "${BS_PROB_PY}" ]] || die "compute_bs_prob.py が見つかりません: ${BS_PROB_PY}"
mkdir -p "${OUTDIR}" "${BSDIR}"

# --- パラメータ (環境変数で上書き可能) --------------------------------------
score="${SCORE:-bdeu}"
ess="${ESS:-10}"
tabu="${TABU:-20}"
iters="${ITERS_BS:-${ITERS:-1000}}"   # リサンプルごとの Hill-Climb 反復数
topk="${TOPK:-20}"
maxpar="${MAX_PARENTS:-3}"
jcache="${JINDEX_CACHE:-1024}"
bootstrap="${BOOTSTRAP:-10}"          # 1 シードあたりのリサンプリング回数
seeds="${SEEDS:-5}"                   # シード数 (総リサンプル数 = bootstrap x seeds)
max_jobs="${MAX_JOBS:-0}"             # 同時実行プロセス数 (0 = 全シード同時)
thr_prob="${THRESHOLD_PROB:-0.2}"     # コンセンサス採用のブートストラップ確率閾値
thr_cnt="${THRESHOLD_COUNT:-2}"       # コンセンサス採用の出現回数閾値
warm_start="${WARM_START:-1}"         # 1 で out/edges.tsv を初期構造に使う
score_imp="${SCORE_IMP:-bic}"         # 重要度評価に用いるスコア

init_opt=()
if [[ "${warm_start}" != "0" && -f "${OUTDIR}/edges.tsv" ]]; then
  init_opt=(--init "${OUTDIR}/edges.tsv")
  log "初期構造 ${OUTDIR}/edges.tsv を使用 (warm start)"
else
  log "初期構造なし (各リサンプルを de novo で学習)"
fi

total=$(( bootstrap * seeds ))
log "総リサンプル数 = ${bootstrap} x ${seeds} = ${total} (iters=${iters}, score=${score})"
log "注意: 1 リサンプルあたり概ね iters x (1 ノード走査時間) の計算量。"
log "      まず BOOTSTRAP/SEEDS/ITERS を小さくして所要時間を確認してください。"

# --- A) 並列ブートストラップ学習 --------------------------------------------
if [[ "${max_jobs}" -gt 0 && "${max_jobs}" -lt "${seeds}" ]]; then
  hr "A) ブートストラップ・リサンプリング (${seeds} シード, 同時実行最大 ${max_jobs})"
else
  hr "A) ブートストラップ・リサンプリング (${seeds} シード並列)"
fi
rm -f "${BSDIR}"/edges_seed*.tsv
fail=0
running=0
for seed in $(seq 1 "${seeds}"); do
  if [[ "${max_jobs}" -gt 0 && "${running}" -ge "${max_jobs}" ]]; then
    wait -n || fail=1
    running=$(( running - 1 ))
  fi
  "${FASTBN_BIN}" "${init_opt[@]}" \
    --input "${INPUT}" \
    --score "${score}" --ess "${ess}" \
    --bootstrap "${bootstrap}" \
    --save-bootstrap-counts "${BSDIR}/edges.tsv" \
    --max-parents "${maxpar}" \
    --topk "${topk}" --jindex-cache "${jcache}" --tabu "${tabu}" --iters "${iters}" \
    --seed "${seed}" > "${BSDIR}/log_seed${seed}.txt" 2>&1 &
  running=$(( running + 1 ))
done
while [[ "${running}" -gt 0 ]]; do
  wait -n || fail=1
  running=$(( running - 1 ))
done
[[ "${fail}" -eq 0 ]] || die "ブートストラップ実行に失敗しました。${BSDIR}/log_seed*.txt を確認してください。"
log "出力: ${BSDIR}/edges_seed*.tsv ($(ls -1 "${BSDIR}"/edges_seed*.tsv | wc -l) ファイル)"

# --- B) 統合 -> コンセンサス構造 --------------------------------------------
hr "B) コンセンサス構造の抽出 (compute_bs_prob.py; prob>=${thr_prob}, count>=${thr_cnt})"
"${PYTHON_BIN}" "${BS_PROB_PY}" \
  --input "${BSDIR}"/edges_seed*.tsv \
  --out-edge "${OUTDIR}/integ_edges.tsv" \
  --out      "${OUTDIR}/integ_edges_score.tsv" \
  --threshold-prob "${thr_prob}" --threshold-count "${thr_cnt}" \
  --remove-cycle --sort-by-prob
n_integ=$(wc -l < "${OUTDIR}/integ_edges.tsv")
[[ "${n_integ}" -gt 0 ]] || die "コンセンサスエッジが 0 件でした。THRESHOLD_PROB を下げるか、BOOTSTRAP/SEEDS を増やしてください。"
log "コンセンサスエッジ: ${n_integ} 件 -> ${OUTDIR}/integ_edges.tsv (頻度つき: integ_edges_score.tsv)"

# --- C) コンセンサス構造のカウント表を再計算 (探索なし: iters=0) ------------
hr "C) コンセンサス構造のカウント表を再計算 (iters=0)"
"${FASTBN_BIN}" --input "${INPUT}" --score "${score}" --ess "${ess}" \
  --init "${OUTDIR}/integ_edges.tsv" --iters 0 --topk "${topk}" \
  --jindex-cache "${jcache}" \
  --save        "${OUTDIR}/integ_edges2.tsv" \
  --save-names  "${OUTDIR}/integ_edges_named.tsv" \
  --save-counts "${OUTDIR}/integ_all_counts.tsv" \
  > "${OUTDIR}/log_integ.txt" 2>&1
log "出力: ${OUTDIR}/integ_edges2.tsv, integ_edges_named.tsv, integ_all_counts.tsv"

# --- D) コンセンサス網のエッジ重要度 ----------------------------------------
hr "D) コンセンサス網のエッジ重要度 (Edge Importance Mode)"
INIT="${OUTDIR}/integ_edges.tsv" \
COUNTS="${OUTDIR}/integ_all_counts.tsv" \
OUT_IMP="${OUTDIR}/integ_edge_importance.tsv" \
SCORE_IMP="${score_imp}" \
  "${BN_SCRIPT_DIR}/edge_importance.sh"

hr "完了。主な出力:"
echo "   ${BSDIR}/edges_seed*.tsv          : シードごとのエッジ出現頻度"
echo "   ${OUTDIR}/integ_edges_score.tsv   : 統合エッジ (u v count prob)"
echo "   ${OUTDIR}/integ_edges_named.tsv   : コンセンサス網 (遺伝子名)"
echo "   ${OUTDIR}/integ_edge_importance.tsv : コンセンサス網のエッジ重要度"
echo "=================================================================="

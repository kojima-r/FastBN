#!/usr/bin/env bash
# =============================================================================
# 02learn.sh — ステップ 2: 構造学習
#   ${NETLIST} の全ネットワーク x SCORES で fast_bn を実行する。
#
#   反復数は変数数に応じて自動で決める:
#       ITERS = max(ITERS_MIN, ITERS_PER_VAR x 変数数)
#   Hill-Climb は成長段階で 1 反復 1 辺しか足せないので、固定反復数だと大きい網が
#   「反復上限で打ち切られた構造」になってしまうため。各実行が収束したか
#   打ち切られたかをログに出す。ITERS を明示すると固定値で上書きできる。
#
# 出力: ${OUTROOT}/<network>_<score>/edges.tsv ほか
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
[[ -s "${NETLIST}" ]] || { echo "[02learn] エラー: ${NETLIST} がありません。先に ./01prepare.sh" >&2; exit 1; }

mkdir -p "${OUTROOT}"
start=$(date +%s); count=0
# DATASETS に含まれないデータセットは飛ばす (NETLIST は 01prepare 時点の全件)
in_datasets() { local d; for d in ${DATASETS}; do [[ "$d" == "$1" ]] && return 0; done; return 1; }
while IFS=$'\t' read -r ds net; do
  in_datasets "${ds}" || continue
  data="${DATADIR}/${net}.tsv"
  [[ -s "${data}" ]] || { echo "[02learn] 警告: ${data} が無いのでスキップ" >&2; continue; }
  nv=$(head -1 "${data}" | awk -F'\t' '{print NF}')
  # 反復数を変数数に応じて決める (ITERS が明示されていればそれを使う)
  if [[ -n "${ITERS:-}" ]]; then
    iters="${ITERS}"
  else
    iters=$(( ITERS_PER_VAR * nv ))
    [[ "${iters}" -lt "${ITERS_MIN}" ]] && iters="${ITERS_MIN}"
  fi
  for score in ${SCORES}; do
    run="${net}_${score}"
    t0=$(date +%s)
    INPUT="${data}" OUTDIR="${OUTROOT}/${run}" SCORE="${score}" ITERS="${iters}" \
      "${BN_SCRIPTS}/learn_structure.sh" > /dev/null 2>&1 \
      || { echo "[02learn] エラー: ${run} の学習に失敗しました" >&2; exit 1; }
    # 反復予算が足りているかの判定。tabu 探索は局所最適でも動き続けるので
    # [stop] はまず出ない。代わりに「最後の 10% の反復でまだ最良スコアが
    # 更新されたか」を見る。更新されていれば反復を増やす余地がある。
    state=$("${BN_SCRIPTS}/../example_dream/iter_state.sh" "${OUTROOT}/${run}/log_learn.txt")
    count=$(( count + 1 ))
    echo "[02learn] ${run} (${nv} 変数, iters=${iters}, $(( $(date +%s) - t0 )) 秒, ${state})"
  done
done < "${NETLIST}"
echo "[02learn] 完了: ${count} 実行 / $(( $(date +%s) - start )) 秒 -> ${OUTROOT}/"

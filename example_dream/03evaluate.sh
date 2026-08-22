#!/usr/bin/env bash
# =============================================================================
# 03evaluate.sh — ステップ 3: 正解ネットワークとの比較
#   ../script/evaluate_structure.py で SHD / Directed・Skeleton の P/R/F1 を計算。
#
#   評価対象ペア: gold standard が判定したペア (<network>_pairs.tsv) だけで
#   Precision / Recall を計算する。DREAM5 の gold standard は TF x 遺伝子の
#   一部しか判定していないため、これをしないと Precision が過小評価される。
#   判定対象外だった学習エッジの本数は n_pred_not_evaluable 列に出る。
#
#   iters_used / budget_binding 列: 実際に回した反復数と、反復予算が律速して
#        いるか (1 = 最後の 10% の反復でまだ最良スコアが更新された)。
#        budget_binding=1 の指標は「探索性能」ではなく「反復予算」を反映している。
#
#   SID: 正解が DAG でない (フィードバックあり) 場合と、ノード数が
#        MAX_SID_NODES を超える場合は NA。
#   KL : 真の CPT が無いので計算しない。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
[[ -s "${NETLIST}" ]] || { echo "[03evaluate] エラー: ${NETLIST} がありません" >&2; exit 1; }

mkdir -p "${EVALDIR}"
rm -f "${BENCHMARK}"
count=0
# DATASETS に含まれないデータセットは飛ばす (NETLIST は 01prepare 時点の全件)
in_datasets() { local d; for d in ${DATASETS}; do [[ "$d" == "$1" ]] && return 0; done; return 1; }
while IFS=$'\t' read -r ds net; do
  in_datasets "${ds}" || continue
  data="${DATADIR}/${net}.tsv"
  edges_true="${TRUTHDIR}/${net}_edges.tsv"
  pairs="${TRUTHDIR}/${net}_pairs.tsv"
  [[ -s "${edges_true}" ]] || { echo "[03evaluate] 警告: ${edges_true} が無いのでスキップ" >&2; continue; }
  n=$(( $(wc -l < "${data}") - 1 ))
  nv=$(head -1 "${data}" | awk -F'\t' '{print NF}')
  pair_opt=()
  [[ -s "${pairs}" ]] && pair_opt=(--eval-pairs "${pairs}")
  for score in ${SCORES}; do
    run="${net}_${score}"
    pred="${OUTROOT}/${run}/edges.tsv"
    [[ -f "${pred}" ]] || { echo "[03evaluate] 警告: ${pred} が無いのでスキップ" >&2; continue; }
    # 反復上限で打ち切られたか (収束していれば log に [stop] が出る) も記録する。
    # 大きい網では「反復予算が足りていないだけ」の可能性を指標と一緒に見たいため。
    lg="${OUTROOT}/${run}/log_learn.txt"
    used_iters=$(grep -c "^\[it " "${lg}" 2>/dev/null); used_iters="${used_iters:-0}"
    # budget_binding=1 なら「最後の 10% の反復でまだ最良スコアが更新された」=
    # 反復を増やせばまだ良くなる = 指標は反復予算に律速されている
    binding=$("${BN_SCRIPTS}/../example_dream/iter_state.sh" --flag "${lg}")
    python3 "${BN_SCRIPTS}/evaluate_structure.py" \
      --true-edges "${edges_true}" --pred-edges "${pred}" --input "${data}" \
      "${pair_opt[@]}" --max-sid-nodes "${MAX_SID_NODES}" \
      --out "${EVALDIR}/${run}.tsv" --append "${BENCHMARK}" \
      --extra "dataset=${ds}" --extra "network=${net}" --extra "score=${score}" \
      --extra "n=${n}" --extra "n_vars=${nv}" \
      --extra "iters_used=${used_iters}" --extra "budget_binding=${binding}" --quiet
    count=$(( count + 1 ))
    printf "\r[03evaluate] %d 実行を評価しました          " "${count}"
  done
done < "${NETLIST}"
echo
echo "[03evaluate] 完了: ${count} 行 -> ${BENCHMARK}"

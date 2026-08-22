#!/usr/bin/env bash
# =============================================================================
# 03evaluate.sh — ステップ 3: 正解シグナル伝達経路との比較
#   ../script/evaluate_structure.py で SHD / Directed・Skeleton の P/R/F1 を計算し、
#   ${BENCHMARK} に 1 行ずつ追記する。
#
#   SID: 正解構造が PKA <-> PIP3 の相互作用を含み DAG ではないため NA になる。
#   KL : 真の CPT が与えられていないため計算しない (--true-bif が無い)。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${EVALDIR}" "$(dirname "${BENCHMARK}")"
rm -f "${BENCHMARK}"
[[ -s "${TRUE_EDGES}" ]] || { echo "[03evaluate] エラー: ${TRUE_EDGES} がありません。先に ./01prepare.sh" >&2; exit 1; }

count=0
for preset in ${PRESETS}; do
  for bins in ${BINS_LIST}; do
    data="${DATADIR}/sachs_${preset}_b${bins}.tsv"
    n=$(( $(wc -l < "${data}") - 1 ))
    for score in ${SCORES}; do
      run="${preset}_b${bins}_${score}"
      edges="${OUTROOT}/${run}/edges.tsv"
      [[ -f "${edges}" ]] || { echo "[03evaluate] エラー: ${edges} がありません。先に ./02learn.sh" >&2; exit 1; }
      python3 "${BN_SCRIPTS}/evaluate_structure.py" \
        --true-edges "${TRUE_EDGES}" \
        --pred-edges "${edges}" \
        --input "${data}" \
        --out "${EVALDIR}/${run}.tsv" \
        --append "${BENCHMARK}" \
        --extra "preset=${preset}" --extra "bins=${bins}" \
        --extra "score=${score}" --extra "n=${n}" \
        --quiet
      count=$(( count + 1 ))
      printf "\r[03evaluate] %d 実行を評価しました          " "${count}"
    done
  done
done
echo
echo "[03evaluate] 完了: ${count} 行 -> ${BENCHMARK}"

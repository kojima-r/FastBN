#!/usr/bin/env bash
# =============================================================================
# 01prepare.sh — ステップ 1: データセットの構築
#   条件の組み合わせ (PRESETS) と離散化の段階数 (BINS_LIST) ごとに
#   fast_bn 入力を作る。正解エッジ表は 1 回だけ書き出す。
#
# 出力 (${DATADIR}):
#   sachs_<preset>_b<bins>.tsv          fast_bn 入力 (整数コード)
#   sachs_<preset>_b<bins>_samples.tsv  行 -> 実験条件 の対応
#   sachs_<preset>_b<bins>_varmap.tsv   列 -> タンパク質名
#   true_edges.tsv                      正解エッジ (from -> to)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${DATADIR}"
log2_opt=()
[[ "${USE_LOG2}" != "0" ]] && log2_opt=(--log2)

first=1
for preset in ${PRESETS}; do
  for bins in ${BINS_LIST}; do
    tag="sachs_${preset}_b${bins}"
    edges_opt=()
    if [[ "${first}" == "1" ]]; then
      edges_opt=(--out-edges "${TRUE_EDGES}")
      first=0
    fi
    python3 ./prepare_sachs.py \
      --data-dir "${SRCDIR}/Data Files" \
      --ground-truth "${SRCDIR}/GroundTruth.csv" \
      --preset "${preset}" --bins "${bins}" --method "${DISC_METHOD}" \
      "${log2_opt[@]}" \
      --out "${DATADIR}/${tag}.tsv" \
      --out-samples "${DATADIR}/${tag}_samples.tsv" \
      --out-map "${DATADIR}/${tag}_varmap.tsv" \
      "${edges_opt[@]}"
  done
done

echo "[01prepare] 完了: $(ls -1 "${DATADIR}"/sachs_*_b*.tsv | grep -vc _samples\\\|_varmap) データセット -> ${DATADIR}/"

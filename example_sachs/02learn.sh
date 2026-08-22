#!/usr/bin/env bash
# =============================================================================
# 02learn.sh — ステップ 2: 構造学習
#   PRESETS x BINS_LIST x SCORES の全組み合わせで fast_bn を実行する。
#   1 実行につき ../script/learn_structure.sh を呼ぶだけ。
#
# 出力: ${OUTROOT}/<preset>_b<bins>_<score>/edges.tsv ほか
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${OUTROOT}"
start=$(date +%s); count=0
for preset in ${PRESETS}; do
  for bins in ${BINS_LIST}; do
    data="${DATADIR}/sachs_${preset}_b${bins}.tsv"
    [[ -s "${data}" ]] || { echo "[02learn] エラー: ${data} がありません。先に ./01prepare.sh" >&2; exit 1; }
    for score in ${SCORES}; do
      run="${preset}_b${bins}_${score}"
      INPUT="${data}" OUTDIR="${OUTROOT}/${run}" SCORE="${score}" \
        "${BN_SCRIPTS}/learn_structure.sh" > /dev/null 2>&1 \
        || { echo "[02learn] エラー: ${run} の学習に失敗しました" >&2; exit 1; }
      count=$(( count + 1 ))
      printf "\r[02learn] %d 実行完了 (最新: %s)          " "${count}" "${run}"
    done
  done
done
echo
echo "[02learn] 完了: ${count} 実行 / $(( $(date +%s) - start )) 秒 -> ${OUTROOT}/"

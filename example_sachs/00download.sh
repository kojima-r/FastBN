#!/usr/bin/env bash
# =============================================================================
# 00download.sh — ステップ 0: Sachs データの取得
#   Zenodo (https://zenodo.org/records/7681811, CC-BY-4.0) から
#   sachs.zip を取得して ${SRCDIR} に展開する。既にあればスキップ。
#
# 展開されるもの:
#   ${SRCDIR}/Data Files/*.csv   14 実験条件ごとの単一細胞測定値 (11 タンパク質)
#   ${SRCDIR}/GroundTruth.csv    既知のシグナル伝達経路 (from, to)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

if [[ -s "${SRCDIR}/GroundTruth.csv" ]]; then
  echo "[00download] ${SRCDIR} は取得済み (スキップ)"
else
  mkdir -p "${SRCDIR}"
  echo "[00download] 取得: ${SACHS_URL}"
  tmp="${SRCDIR}/sachs.zip"
  if command -v curl > /dev/null; then
    curl -fsSL --retry 3 -C - -o "${tmp}" "${SACHS_URL}"
  elif command -v wget > /dev/null; then
    wget -qc -O "${tmp}" "${SACHS_URL}"
  else
    echo "[00download] エラー: curl も wget もありません" >&2; exit 1
  fi
  unzip -oq "${tmp}" -d "${SRCDIR}"
  rm -f "${tmp}"
fi

[[ -s "${SRCDIR}/GroundTruth.csv" ]] || { echo "[00download] エラー: GroundTruth.csv がありません" >&2; exit 1; }
n_cond=$(ls -1 "${SRCDIR}/Data Files/"*.csv 2>/dev/null | wc -l)
n_edge=$(( $(wc -l < "${SRCDIR}/GroundTruth.csv") - 1 ))
echo "[00download] 完了: ${n_cond} 条件のデータ / 正解エッジ ${n_edge} 本"

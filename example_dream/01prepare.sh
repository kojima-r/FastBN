#!/usr/bin/env bash
# =============================================================================
# 01prepare.sh — ステップ 1: データセットの構築
#   各ネットワークについて、発現量行列を離散化した fast_bn 入力と、
#   正解エッジ・評価対象ペアを書き出す。
#
# 出力:
#   ${DATADIR}/<network>.tsv          fast_bn 入力 (整数コード)
#   ${DATADIR}/<network>_varmap.tsv   列 -> 遺伝子名
#   ${TRUTHDIR}/<network>_edges.tsv   正解エッジ
#   ${TRUTHDIR}/<network>_pairs.tsv   gold standard が判定したペア (評価の分母)
#   ${NETLIST}                        作成したネットワーク名の一覧
#
# DREAM4 は複数の摂動実験を縦に結合して 1 つの行列にします (D4_PARTS)。
# DREAM5 は数千遺伝子あるので既定で D5_MAX_VARS 個に絞ります (TF を優先)。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${DATADIR}" "${TRUTHDIR}" "$(dirname "${NETLIST}")"
: > "${NETLIST}"

for ds in ${DATASETS}; do
  echo "[01prepare] ${ds} を準備します"
  tmp="${RUNDIR}/.list_${ds}.txt"
  case "${ds}" in
    dream4)
      python3 ./prepare_dream.py --dataset dream4 --zip "${SRCDIR}/dream4.zip" \
        --out-dir "${DATADIR}" --truth-dir "${TRUTHDIR}" \
        --bins "${BINS}" --method "${DISC_METHOD}" \
        --sizes "${D4_SIZES}" --parts "${D4_PARTS}" --out-list "${tmp}" ;;
    dream5)
      python3 ./prepare_dream.py --dataset dream5 --zip "${SRCDIR}/dream5.zip" \
        --out-dir "${DATADIR}" --truth-dir "${TRUTHDIR}" \
        --bins "${BINS}" --method "${DISC_METHOD}" \
        --networks "${D5_NETWORKS}" --max-vars "${D5_MAX_VARS}" --out-list "${tmp}" ;;
    hpn)
      python3 ./prepare_dream.py --dataset hpn --src-dir "${HPN_SRCDIR}" \
        --out-dir "${DATADIR}" --truth-dir "${TRUTHDIR}" \
        --bins "${BINS}" --method "${DISC_METHOD}" \
        --max-vars "${D5_MAX_VARS}" --out-list "${tmp}" ;;
  esac
  while read -r net; do
    [[ -n "${net}" ]] && printf '%s\t%s\n' "${ds}" "${net}" >> "${NETLIST}"
  done < "${tmp}"
  rm -f "${tmp}"
done

echo "[01prepare] 完了: $(wc -l < "${NETLIST}") ネットワーク -> ${NETLIST}"

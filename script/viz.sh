#!/usr/bin/env bash
# =============================================================================
# viz.sh
#   Hill-Climb 学習網 (${OUTDIR}/edges.tsv, edges_named.tsv, edge_importance.tsv)
#   を可視化する簡易呼び出しスクリプト (汎用版)。
#   図は ${FIGDIR} (既定 ./figures) に出力する。
#   追加オプションはそのまま visualize.py へ渡される
#   (例: ../script/viz.sh --metrics dBIC --top-n 80 --layout kamada)。
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="viz"; bn_tag="viz"

FIGDIR="${FIGDIR:-./figures}"
require_file "${OUTDIR}/edges.tsv" "${OUTDIR}/edges_named.tsv"

"${PYTHON_BIN}" "$(py_tool visualize.py)" \
  --out-dir "${OUTDIR}" \
  --fig-dir "${FIGDIR}" \
  --var-map "${VARMAP}" \
  --target-file "${TARGET_FILE}" \
  --edges        edges.tsv \
  --edges-named  edges_named.tsv \
  --importance   edge_importance.tsv \
  "$@"

log "学習網の図を ${FIGDIR}/ に出力しました"

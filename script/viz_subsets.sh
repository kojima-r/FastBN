#!/usr/bin/env bash
# =============================================================================
# viz_subsets.sh
#   Hill-Climb 学習網について、群 (サブセット) ごとのエッジ重要度を
#   「全体グラフを背景に薄く + その群で重要なエッジのみ強調」で描画する。
#   ノード配置は全群で共通。図は ${FIGDIR}/subsets (既定 ./figures/subsets) へ。
#   追加オプションは viz_subsets.py へ素通し (例: --top-n 60 --metrics dBIC)。
#
#   前提: importance_groups.sh で ${OUTDIR}/edge_importance_g*_*.tsv を作成済み。
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="viz_subsets"; bn_tag="viz_subsets"

FIGDIR="${FIGDIR:-./figures}"
require_file "${OUTDIR}/edges.tsv" "${OUTDIR}/edges_named.tsv"

"${PYTHON_BIN}" "$(py_tool viz_subsets.py)" \
  --out-dir "${OUTDIR}" \
  --fig-dir "${FIGDIR}/subsets" \
  --var-map "${VARMAP}" \
  --target-file "${TARGET_FILE}" \
  --edges        edges.tsv \
  --edges-named  edges_named.tsv \
  --prefix       "${IMP_PREFIX:-edge_importance}" \
  --include-all \
  "$@"

log "学習網の群別重要度図を ${FIGDIR}/subsets/ に出力しました"

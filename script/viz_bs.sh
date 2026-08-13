#!/usr/bin/env bash
# =============================================================================
# viz_bs.sh
#   ブートストラップ・コンセンサス網を可視化する簡易呼び出しスクリプト (汎用版)。
#   bootstrap_stability.sh が出力した ${OUTDIR}/integ_* を読み込み、
#   図は ${FIGDIR_BS} (既定 ./figures_bs) に出力する。
#     構造/名前 : integ_edges2.tsv <-> integ_edges_named.tsv (行対応)
#     重要度     : integ_edge_importance.tsv
#     確率       : integ_edges_score.tsv (u v count prob) -> 06_bootstrap_prob.png
#   追加オプションはそのまま visualize.py へ渡される。
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="viz_bs"; bn_tag="viz_bs"

FIGDIR="${FIGDIR_BS:-./figures_bs}"
for f in integ_edges2.tsv integ_edges_named.tsv; do
  [[ -f "${OUTDIR}/${f}" ]] || die "${OUTDIR}/${f} がありません。先に bootstrap_stability.sh を実行してください。"
done

"${PYTHON_BIN}" "$(py_tool visualize.py)" \
  --out-dir "${OUTDIR}" \
  --fig-dir "${FIGDIR}" \
  --var-map "${VARMAP}" \
  --target-file "${TARGET_FILE}" \
  --edges        integ_edges2.tsv \
  --edges-named  integ_edges_named.tsv \
  --importance   integ_edge_importance.tsv \
  --edge-prob    integ_edges_score.tsv \
  "$@"

log "コンセンサス網の図を ${FIGDIR}/ に出力しました"

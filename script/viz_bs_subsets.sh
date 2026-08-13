#!/usr/bin/env bash
# =============================================================================
# viz_bs_subsets.sh
#   ブートストラップ・コンセンサス網について、群 (サブセット) ごとの
#   エッジ重要度を「全体グラフを背景に薄く + 重要部分のみ強調」で描画する。
#   ノード配置は全群で共通。図は ${FIGDIR_BS}/subsets (既定 ./figures_bs/subsets) へ。
#   線の太さにはブートストラップ確率 (integ_edges_score.tsv) を用いる。
#
#   前提: bootstrap_stability.sh と、コンセンサス網に対する importance_groups.sh
#         (OUT_PREFIX=integ_edge_importance) を実行済み。
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="viz_bs_subsets"; bn_tag="viz_bs_subsets"

FIGDIR="${FIGDIR_BS:-./figures_bs}"
for f in integ_edges2.tsv integ_edges_named.tsv; do
  [[ -f "${OUTDIR}/${f}" ]] || die "${OUTDIR}/${f} がありません。先に bootstrap_stability.sh を実行してください。"
done

"${PYTHON_BIN}" "$(py_tool viz_subsets.py)" \
  --out-dir "${OUTDIR}" \
  --fig-dir "${FIGDIR}/subsets" \
  --var-map "${VARMAP}" \
  --target-file "${TARGET_FILE}" \
  --edges        integ_edges2.tsv \
  --edges-named  integ_edges_named.tsv \
  --prefix       "${IMP_PREFIX_BS:-integ_edge_importance}" \
  --edge-prob    integ_edges_score.tsv \
  --include-all \
  "$@"

log "コンセンサス網の群別重要度図を ${FIGDIR}/subsets/ に出力しました"

#!/usr/bin/env bash
# =============================================================================
# 06visualize.sh — ステップ 6: ネットワークの可視化
#   networkx + matplotlib で図を出力する。存在する成果物に応じて自動的に
#   対象を増やす (コンセンサス網・組織別図は該当ファイルがある場合のみ)。
#
# 出力:
#   ${FIGDIR}/            : 学習網の図 (01_structure_full, 02_importance_full ...)
#   ${FIGDIR}/subsets/    : 学習網の組織別比較図
#   ${FIGDIR_BS}/         : コンセンサス網の図 (06_bootstrap_prob を含む)
#   ${FIGDIR_BS}/subsets/ : コンセンサス網の組織別比較図
#   ${FIGDIR}*/edge_importance_named_<metric>.tsv : 重要度降順の名前付きエッジ表
#
# メトリクスや本数は config.sh (VIZ_METRICS / VIZ_TOP_N) で変更できる。
#   VIZ_METRICS=dlogL,dBIC,dK2,dBDeu ./06visualize.sh
# 追加オプションはそのまま visualize.py / viz_subsets.py に渡る。
#   ./06visualize.sh --layout kamada
#
# 注: 既定の NVARS=all では図のノードが数千個になり、全体図はほぼ塗りつぶしになる
#     (描画自体は通る: bbknn/all で約 4 分)。--top-n を絞るか、注目遺伝子
#     (config.sh の TARGET_GENES) を指定した図を見ること。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

common_opt=(--metrics "${VIZ_METRICS}" --top-n "${VIZ_TOP_N}")

echo "=================================================================="
echo " (a) Hill-Climb 学習網 -> ${FIGDIR}/"
echo "=================================================================="
"${BN_SCRIPTS}/viz.sh" "${common_opt[@]}" "$@"

if [[ -f "${OUTDIR}/integ_edges2.tsv" ]]; then
  echo "=================================================================="
  echo " (b) ブートストラップ・コンセンサス網 -> ${FIGDIR_BS}/"
  echo "=================================================================="
  "${BN_SCRIPTS}/viz_bs.sh" "${common_opt[@]}" "$@"
fi

if compgen -G "${OUTDIR}/edge_importance_g*_*.tsv" > /dev/null; then
  echo "=================================================================="
  echo " (c) 組織別重要度 (学習網) -> ${FIGDIR}/subsets/"
  echo "=================================================================="
  "${BN_SCRIPTS}/viz_subsets.sh" "${common_opt[@]}" "$@"
fi

if compgen -G "${OUTDIR}/integ_edge_importance_g*_*.tsv" > /dev/null; then
  echo "=================================================================="
  echo " (d) 組織別重要度 (コンセンサス網) -> ${FIGDIR_BS}/subsets/"
  echo "=================================================================="
  "${BN_SCRIPTS}/viz_bs_subsets.sh" "${common_opt[@]}" "$@"
fi

echo "[06visualize] 完了。生成された図:"
find "${FIGDIR}" "${FIGDIR_BS}" -name '*.png' 2>/dev/null | sort | sed 's/^/   /'

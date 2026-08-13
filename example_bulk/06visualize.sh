#!/usr/bin/env bash
# =============================================================================
# 06visualize.sh — ステップ 6: ネットワークの可視化
#   networkx + matplotlib で図を出力する。存在する成果物に応じて自動的に
#   対象を増やす (コンセンサス網・群別図は該当ファイルがある場合のみ)。
#
# 出力:
#   figures/            : 学習網の図 (01_structure_full, 02_importance_full ...)
#   figures/subsets/    : 学習網の群別比較図
#   figures_bs/         : コンセンサス網の図 (06_bootstrap_prob を含む)
#   figures_bs/subsets/ : コンセンサス網の群別比較図
#   figures*/edge_importance_named_<metric>.tsv : 重要度降順の名前付きエッジ表
#
# メトリクスや本数は config.sh (VIZ_METRICS / VIZ_TOP_N) で変更できる。
#   VIZ_METRICS=dlogL,dBIC,dK2,dBDeu ./06visualize.sh
# 追加オプションはそのまま visualize.py / viz_subsets.py に渡る。
#   ./06visualize.sh --layout kamada
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

common_opt=(--metrics "${VIZ_METRICS}" --top-n "${VIZ_TOP_N}")

echo "=================================================================="
echo " (a) Hill-Climb 学習網 -> figures/"
echo "=================================================================="
"${BN_SCRIPTS}/viz.sh" "${common_opt[@]}" "$@"

if [[ -f "${OUTDIR}/integ_edges2.tsv" ]]; then
  echo "=================================================================="
  echo " (b) ブートストラップ・コンセンサス網 -> figures_bs/"
  echo "=================================================================="
  "${BN_SCRIPTS}/viz_bs.sh" "${common_opt[@]}" "$@"
fi

if compgen -G "${OUTDIR}/edge_importance_g*_*.tsv" > /dev/null; then
  echo "=================================================================="
  echo " (c) 群別重要度 (学習網) -> figures/subsets/"
  echo "=================================================================="
  "${BN_SCRIPTS}/viz_subsets.sh" "${common_opt[@]}" "$@"
fi

if compgen -G "${OUTDIR}/integ_edge_importance_g*_*.tsv" > /dev/null; then
  echo "=================================================================="
  echo " (d) 群別重要度 (コンセンサス網) -> figures_bs/subsets/"
  echo "=================================================================="
  "${BN_SCRIPTS}/viz_bs_subsets.sh" "${common_opt[@]}" "$@"
fi

echo "[06visualize] 完了。生成された図:"
find figures figures_bs -name '*.png' 2>/dev/null | sort | sed 's/^/   /'

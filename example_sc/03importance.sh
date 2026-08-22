#!/usr/bin/env bash
# =============================================================================
# 03importance.sh — ステップ 3: エッジ重要度 (全サンプル)
#   各エッジを除いたときのスコア変化を計算する (Edge Importance Mode)。
#
# 出力:
#   ${OUTDIR}/edge_importance.tsv
#     列: u v ΔlogL ΔBIC ΔK2 ΔBDeu meanΔlogL_per_sample stdΔlogL_per_sample
#   ノード番号 u, v は ${DATADIR}/var_map.tsv の index 列で遺伝子名に対応づく
#   (図・レポートでは自動的に名前に変換される)。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

"${BN_SCRIPTS}/edge_importance.sh" "$@"

echo "[03importance] 上位 5 エッジ (ΔlogL 昇順 = 除去の影響が大きい順):"
tail -n +2 "${OUTDIR}/edge_importance.tsv" \
  | sort -t$'\t' -k3,3g | head -5 | sed 's/^/   /' || true

#!/usr/bin/env bash
# =============================================================================
# 05importance_groups.sh — ステップ 5: 群 (条件) 別のエッジ重要度
#   ネットワーク構造を固定したまま、評価データを群ごとに差し替えて
#   (--score-dataset) エッジ重要度を計算する。
#   -> 「どのエッジがどの条件で効いているか」を比較できる。
#
#   群の定義は data/samples.tsv (前処理が出力した サンプル<->群 の表) から
#   自動的に読み込まれる。
#
# 出力:
#   groups/expr_g{N}_<label>.tsv               : 群ごとの評価データ
#   out/edge_importance_g{N}_<label>.tsv       : 学習網ベースの群別重要度
#   out/integ_edge_importance_g{N}_<label>.tsv : コンセンサス網ベース (04 実行後)
#
# 注意: 群あたりのサンプル数が少ないと ΔlogL のノイズが大きい (logL はサンプル数に
#       比例するため絶対値も小さくなる)。群間の相対比較として解釈すること。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

echo "=================================================================="
echo " (a) Hill-Climb 学習網ベースの群別重要度"
echo "=================================================================="
"${BN_SCRIPTS}/importance_groups.sh" "$@"

if [[ -f "${OUTDIR}/integ_edges.tsv" && -f "${OUTDIR}/integ_all_counts.tsv" ]]; then
  echo "=================================================================="
  echo " (b) ブートストラップ・コンセンサス網ベースの群別重要度"
  echo "=================================================================="
  INIT="${OUTDIR}/integ_edges.tsv" \
  COUNTS="${OUTDIR}/integ_all_counts.tsv" \
  REF_EDGES="${OUTDIR}/integ_edges2.tsv" \
  REF_NAMED="${OUTDIR}/integ_edges_named.tsv" \
  OUT_PREFIX=integ_edge_importance \
    "${BN_SCRIPTS}/importance_groups.sh" "$@"
else
  echo "[05importance_groups] コンセンサス網が無いため (b) はスキップ"
  echo "                      (先に ./04bootstrap.sh を実行すると計算されます)"
fi

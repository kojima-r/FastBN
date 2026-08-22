#!/usr/bin/env bash
# =============================================================================
# 05importance_groups.sh — ステップ 5: 組織 (群) 別のエッジ重要度
#   ネットワーク構造を固定したまま、評価データを組織ごとに差し替えて
#   (--score-dataset) エッジ重要度を計算する。
#   -> 「どのエッジがどの組織で効いているか」を比較できる。
#
#   群の定義は 01prepare.sh が作った ${SAMPLES} の group 列 (= 組織名) から
#   自動的に読み込まれる。
#
# 出力:
#   ${GROUPDIR}/expr_g{N}_<組織>.tsv               : 組織ごとの評価データ
#   ${OUTDIR}/edge_importance_g{N}_<組織>.tsv       : 学習網ベース
#   ${OUTDIR}/integ_edge_importance_g{N}_<組織>.tsv : コンセンサス網ベース (04 実行後)
#
# 注意: 1 組織あたりのサンプル数は 10 件程度と少ないため ΔlogL のノイズが大きく、
#       絶対値も全サンプル版より小さくなる (logL はサンプル数に比例)。
#       組織間の相対比較として解釈すること。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

if [[ ! -f "${SAMPLES}" ]]; then
  echo "[05importance_groups] ${SAMPLES} が無いためスキップします。" >&2
  echo "                      (01prepare.sh が組織ファイルを見つけられませんでした)" >&2
  exit 0
fi

echo "=================================================================="
echo " (a) Hill-Climb 学習網ベースの組織別重要度"
echo "=================================================================="
"${BN_SCRIPTS}/importance_groups.sh" "$@"

if [[ -f "${OUTDIR}/integ_edges.tsv" && -f "${OUTDIR}/integ_all_counts.tsv" ]]; then
  echo "=================================================================="
  echo " (b) ブートストラップ・コンセンサス網ベースの組織別重要度"
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

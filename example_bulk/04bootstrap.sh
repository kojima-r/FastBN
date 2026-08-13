#!/usr/bin/env bash
# =============================================================================
# 04bootstrap.sh — ステップ 4: ブートストラップによるエッジ安定性解析
#   サンプルを復元抽出して構造学習を BOOTSTRAP x SEEDS 回繰り返し、エッジの
#   出現頻度 (ブートストラップ確率) から安定なエッジだけを残した
#   コンセンサスネットワークを作り、その重要度まで計算する。
#
# 出力:
#   bs/edges_seed*.tsv            : シードごとのエッジ出現頻度 (u v count prob)
#   out/integ_edges.tsv           : コンセンサスエッジ (u v)
#   out/integ_edges_score.tsv     : コンセンサスエッジ + 出現頻度
#   out/integ_edges_named.tsv     : コンセンサス網 (遺伝子名)
#   out/integ_all_counts.tsv      : コンセンサス構造のカウント表
#   out/integ_edge_importance.tsv : コンセンサス網のエッジ重要度
#
# 計算時間が気になる場合は BOOTSTRAP / SEEDS / ITERS_BS を小さくして試す。
#   BOOTSTRAP=3 SEEDS=2 ./04bootstrap.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

"${BN_SCRIPTS}/bootstrap_stability.sh" "$@"

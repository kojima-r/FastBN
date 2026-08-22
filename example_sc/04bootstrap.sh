#!/usr/bin/env bash
# =============================================================================
# 04bootstrap.sh — ステップ 4: ブートストラップによるエッジ安定性解析
#   サンプル (細胞集団) を復元抽出して構造学習を BOOTSTRAP x SEEDS 回繰り返し、
#   エッジの出現頻度 (ブートストラップ確率) から安定なエッジだけを残した
#   コンセンサスネットワークを作り、その重要度まで計算する。
#
#   内部の流れ (../script/bootstrap_stability.sh):
#     A) fast_bn --bootstrap        : シード並列でリサンプリング学習
#     B) compute_bs_prob.py         : 全シードを統合しコンセンサスエッジを抽出
#     C) fast_bn --iters 0          : コンセンサス構造のカウント表を再計算
#     D) fast_bn --edge-importance  : コンセンサス網のエッジ重要度
#
# 出力:
#   ${BSDIR}/edges_seed*.tsv            : シードごとのエッジ出現頻度 (u v count prob)
#   ${OUTDIR}/integ_edges.tsv           : コンセンサスエッジ (u v)
#   ${OUTDIR}/integ_edges_score.tsv     : コンセンサスエッジ + 出現頻度
#   ${OUTDIR}/integ_edges2.tsv          : 同じ構造を fast_bn が書き直したもの
#   ${OUTDIR}/integ_edges_named.tsv     : コンセンサス網 (遺伝子名; 上と行対応)
#   ${OUTDIR}/integ_all_counts.tsv      : コンセンサス構造のカウント表
#   ${OUTDIR}/integ_edge_importance.tsv : コンセンサス網のエッジ重要度
#
# 計算時間はここが支配的 (BOOTSTRAP x SEEDS 回の構造学習)。まず小さく試すこと:
#   BOOTSTRAP=3 SEEDS=2 ITERS_BS=500 ./04bootstrap.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

"${BN_SCRIPTS}/bootstrap_stability.sh" "$@"

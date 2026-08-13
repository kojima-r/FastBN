#!/usr/bin/env bash
# =============================================================================
# 02learn.sh — ステップ 2: ネットワーク構造学習
#   fast_bn の Hill-Climb + Tabu サーチで DAG を推定する。
#   スコア・最大親数などは config.sh (SCORE / MAX_PARENTS / ITERS ...) で設定。
#
# 出力:
#   out/edges.tsv       : 学習された DAG (ノード=列インデックス)
#   out/edges_named.tsv : 同じエッジを遺伝子名で表記 (edges.tsv と行対応)
#   out/all_counts.tsv  : カウント表 (CPT 推定・エッジ重要度評価に使う)
#   out/log_learn.txt   : fast_bn のログ
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

"${BN_SCRIPTS}/learn_structure.sh" "$@"

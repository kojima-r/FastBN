#!/usr/bin/env bash
# =============================================================================
# bootstrap.sh — example_sachs の全実験 (離散化方法 2 通り x 最大親数 3 通り)
#   1 掃引あたり 3 条件セット x 3 離散化段階 x 3 スコア = 27 実行、6 掃引で 162 実行。
#   出力は results_<離散化方法>_p<最大親数>/ に分かれる。
#
# -----------------------------------------------------------------------------
# ITERS の決め方 — 変数数 D とサンプル数 N から決める
# -----------------------------------------------------------------------------
# (1) Hill-Climb は成長段階では 1 反復に 1 辺しか足せない。したがって必要な
#     反復数の下限は「学習される辺の本数 |E|」そのものである。
# (2) |E| <= D x P_eff。P_eff は 1 ノードが実際に持てる親の数で
#         P_eff = min(MAX_PARENTS, floor(log_r(N / 10)))
#     r は 1 変数の状態数 (= 離散化の段階数)。親を P 個持つノードの CPT は
#     r^P 通りの親設定を持ち、1 設定あたり 10 サンプル程度は無いと罰則付き
#     スコアはその親を保持しない。
# (3) 成長後の修正 (REMOVE / REVERSE) と Tabu が局所最適から抜ける分に、
#     成長段階の 2 倍を見込む:   ITERS = 3 x D x P_eff
#
#   Sachs は 11 タンパク質なので D = 11 で固定、サンプルは最少の obs でも 853 細胞
#   ある。一番厳しい組み合わせ (4 段階離散化 x obs) でも
#       P_eff = min(4, floor(log4(853/10))) = min(4, 3) = 3
#   なので必要反復は ITERS = 3 x 11 x 3 = 99。取りうる有向辺が 11 x 10 = 110 本、
#   正解の辺が 20 本・最大入次数 3 であることと整合する。
#   既定の ITERS=3000 はその 30 倍あり、下げても上げても結果は変わらないため
#   既定のまま使う (サンプルの多い all = 11671 細胞・2 段階なら P_eff = 4 まで
#   上がるが、それでも必要反復は 132)。
#
#   ブートストラップを行わない例題なので ITERS_BS は使わない。
# =============================================================================

# 全掃引に共通の設定 (離散化方法と最大親数だけを下の行で振る)
export PRESETS="obs int all"     # 刺激のみ 853 / 介入 10818 / 全 11671 細胞
export BINS_LIST="2 3 4"         # 原論文は 3 段階
export SCORES="bic bdeu k2"
export ITERS=3000                # = 既定値。上の規則からの必要量 99〜132 の 20 倍以上
export TABU=10

RUNDIR=./results_quantile_p2 DISC_METHOD=quantile MAX_PARENTS=2 ./run_all.sh
RUNDIR=./results_quantile_p3 DISC_METHOD=quantile MAX_PARENTS=3 ./run_all.sh   # 正解の最大入次数と同じ
RUNDIR=./results_quantile_p4 DISC_METHOD=quantile MAX_PARENTS=4 ./run_all.sh
RUNDIR=./results_uniform_p2  DISC_METHOD=uniform  MAX_PARENTS=2 ./run_all.sh
RUNDIR=./results_uniform_p3  DISC_METHOD=uniform  MAX_PARENTS=3 ./run_all.sh
RUNDIR=./results_uniform_p4  DISC_METHOD=uniform  MAX_PARENTS=4 ./run_all.sh

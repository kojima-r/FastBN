#!/usr/bin/env bash
# =============================================================================
# bootstrap.sh — example_dream の全実験 (DREAM4 / DREAM5 300 変数 / DREAM5 全遺伝子)
#   nohup ./bootstrap.sh > bootstrap.log 2>&1 &
#   出力は results_<パス名>/ に分かれる。
#
# -----------------------------------------------------------------------------
# 反復数の決め方 — 変数数 D とサンプル数 N から決める
# -----------------------------------------------------------------------------
# (1) Hill-Climb は成長段階では 1 反復に 1 辺しか足せない。したがって必要な
#     反復数の下限は「学習される辺の本数 |E|」そのものである。
# (2) |E| <= D x P_eff。P_eff は 1 ノードが実際に持てる親の数で
#         P_eff = min(MAX_PARENTS, floor(log_r(N / 10)))
#     r は 1 変数の状態数 (= BINS)。親を P 個持つノードの CPT は r^P 通りの
#     親設定を持ち、1 設定あたり 10 サンプル程度は無いと罰則付きスコアは
#     その親を保持しない。
# (3) 成長後の修正 (REMOVE / REVERSE) と Tabu が局所最適から抜ける分に、
#     成長段階の 2 倍を見込む:   ITERS = 3 x D x P_eff
#
#   この例題は D が 10 から 5950 まで 3 桁ちがうので、ITERS を固定せず
#       ITERS = max(ITERS_MIN, ITERS_PER_VAR x D)
#   と自動スケールさせる。上の (3) は D に比例する形なので、規則はそのまま
#       ITERS_PER_VAR = 3 x P_eff
#   と読み替えられる。BINS=3 で各データのサンプル数は
#       dream4 size10   N=136  -> P_eff = min(3, floor(log3 13.6)) = 2
#       dream4 size100  N=411  -> P_eff = min(3, floor(log3 41.1)) = 3
#       dream5 net1/3   N=805  -> P_eff = min(3, floor(log3 80.5)) = 3
#       dream5 net4     N=536  -> P_eff = min(3, floor(log3 53.6)) = 3
#   なので 3 x P_eff = 6〜9、既定の ITERS_PER_VAR=10 がそのまま妥当な値になる
#   (ITERS_MIN=2000 は D=10 の size10 で必要量 60 を大きく上回る)。よって
#   下の行では反復数を上書きしない (ITERS を渡すと自動スケールが固定値で
#   潰れてしまうので渡してはいけない)。
#
#   結果として各パスの反復数は
#       dream4 size10 / size100 : 2000 / 2000        (必要量 60 / 900)
#       dream5 300 変数         : 3000               (必要量 2700)
#       dream5 全遺伝子         : 16430/45110/59500  (必要量 14787/40599/53550)
#   となり、いずれも必要量を満たす。実際に満たせたかは benchmark.tsv の
#   budget_binding 列 (1 = 反復を増やせばまだ改善する) で確認できる。
#   なお dream4 で budget_binding=1 が出ることがあるが、必要量に対して反復は
#   2 倍以上あるので、これは Tabu の彷徨いによる終盤の微小更新である。
#
#   ブートストラップを行わない例題なので ITERS_BS は使わない。
# =============================================================================

# 全パスに共通の設定 (データセットだけを下の行で振る)
export SCORES="bic bdeu k2"
export BINS=3                # r = 3。上の P_eff の計算はこの値に依存する
export DISC_METHOD="quantile"

RUNDIR=./results_dream4      DATASETS=dream4 D4_SIZES="10 100"                   ./run_all.sh
RUNDIR=./results_dream5_300  DATASETS=dream5 D5_NETWORKS="1 3 4" D5_MAX_VARS=300 ./run_all.sh
RUNDIR=./results_dream5_full DATASETS=dream5 D5_NETWORKS="1 3 4" D5_MAX_VARS=0   ./run_all.sh

# HPN-DREAM は Synapse (syn1720047) の認証が要るため自動取得できない。
# source/hpn/ に手動配置した場合だけ次の行のコメントを外す。
#RUNDIR=./results_hpn        DATASETS=hpn    D5_MAX_VARS=0                       ./run_all.sh

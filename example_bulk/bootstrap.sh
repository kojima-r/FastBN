#!/usr/bin/env bash
# =============================================================================
# bootstrap.sh — example_bulk の全実験
#   データの乱数シード 5 通り x スコア 3 通り x 離散化段階 2 通り x 最大親数 2 通り
#   = 60 ケースを、それぞれフルパイプライン (00 前処理 〜 08 真の構造との比較) で回す。
#     nohup ./bootstrap.sh > bootstrap.log 2>&1 &
#   出力は run_s<seed>_<score>_b<段階>_p<最大親数>/ に分かれる。
#   ダミーデータは真の DAG が既知なので、各ケースの out/eval_hc.tsv (学習網) と
#   out/eval_bs.tsv (コンセンサス網) を並べれば設定の良し悪しを直接比較できる。
#
# -----------------------------------------------------------------------------
# ITERS / ITERS_BS の決め方 — 変数数 D とサンプル数 N から決める
# -----------------------------------------------------------------------------
# (1) Hill-Climb は成長段階では 1 反復に 1 辺しか足せない。したがって必要な
#     反復数の下限は「学習される辺の本数 |E|」そのものである。
# (2) |E| <= D x P_eff。P_eff は 1 ノードが実際に持てる親の数で
#         P_eff = min(MAX_PARENTS, floor(log_r(N / 10)))
#     r は 1 変数の状態数 (= N_BINS)。親を P 個持つノードの CPT は r^P 通りの
#     親設定を持ち、1 設定あたり 10 サンプル程度は無いと罰則付きスコアは
#     その親を保持しない。
# (3) 成長後の修正 (REMOVE / REVERSE) と Tabu が局所最適から抜ける分に、
#     成長段階の 2 倍を見込む:   ITERS = 3 x D x P_eff
# (4) ブートストラップは out/edges.tsv からの warm start なので成長段階が要らず、
#     リサンプル間の差分を直せればよい:   ITERS_BS = D x P_eff  (= ITERS の 1/3)
# (5) リサンプル総数 B = BOOTSTRAP x SEEDS は、エッジ出現確率の標準誤差
#     sqrt(p(1-p)/B) <= 0.5/sqrt(B) で決める。
#
#   このデータは D = TOP_VAR_GENES = 60 変数、N = 4 群 x 30 反復 = 120 サンプル。
#     N_BINS=2 (r=2): P_eff = min(3, floor(log2 12)) = 3 -> ITERS 540 / ITERS_BS 180
#     N_BINS=3 (r=3): P_eff = min(3, floor(log3 12)) = 2 -> ITERS 360 / ITERS_BS 120
#   真の DAG の辺が 60 本前後 (DUMMY_EDGE_PROB=0.06 x DUMMY_MAX_PARENTS=2) なので
#   妥当な桁である。既定の ITERS=5000 / ITERS_BS=1500 は必要量の 9 倍 / 8 倍あり、
#   下げても上げても結果は変わらないため既定のまま使う。
#   3 段階の方が必要反復が少ないのはサンプル数のため: N=120 で 3 値だと親 3 個は
#   3^3 = 27 通りの親設定に対し 4.4 サンプル/設定しかなく、スコアが 3 個目の親を
#   保持しない (2 値なら 2^3 = 8 通りで 15 サンプル/設定あり、親 3 個を支えられる)。
#
#   リサンプル総数は B = 10 x 20 = 200 とした。SE <= 0.035 で、THRESHOLD_PROB=0.3 の
#   採否が +-0.07 (2SE) の精度で決まる。N=120 と小さくリサンプル 1 回が安いので、
#   sc のような大きなデータ (B=100) より分解能を上げてある。
#   (MAX_JOBS は結果に影響しない実行上の都合。CPU コア数で決める)
#
# 軽く試すときは 1 ケースだけ回す:
#   RUNDIR=./run_test DUMMY_SEED=1 SCORE=bdeu N_BINS=3 MAX_PARENTS=2 ./run_all.sh
# =============================================================================

# 全ケースに共通の設定 (下の行では乱数シード・スコア・段階・最大親数だけを振る)
export ITERS=5000        # = 既定値。上の規則からの必要量 360〜540 の 9 倍以上
export ITERS_BS=1500     # = 既定値。同 120〜180 の 8 倍以上
export TABU=30
export BOOTSTRAP=10                          # 1 シードあたりのリサンプル回数
export SEEDS=20                              # 総リサンプル数 B = 10 x 20 = 200
export MAX_JOBS=20                           # 同時実行プロセス数 (= SEEDS で頭打ち)
export VIZ_METRICS="dlogL,dBIC,dK2,dBDeu"    # 全指標の図を作る

# --- データ乱数シード 1 (真の DAG を引き直した反復) ---
RUNDIR=./run_s1_bic_b2_p2   DUMMY_SEED=1 SCORE=bic  N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s1_bic_b2_p3   DUMMY_SEED=1 SCORE=bic  N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s1_bic_b3_p2   DUMMY_SEED=1 SCORE=bic  N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s1_bic_b3_p3   DUMMY_SEED=1 SCORE=bic  N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s1_bdeu_b2_p2  DUMMY_SEED=1 SCORE=bdeu N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s1_bdeu_b2_p3  DUMMY_SEED=1 SCORE=bdeu N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s1_bdeu_b3_p2  DUMMY_SEED=1 SCORE=bdeu N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s1_bdeu_b3_p3  DUMMY_SEED=1 SCORE=bdeu N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s1_k2_b2_p2    DUMMY_SEED=1 SCORE=k2   N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s1_k2_b2_p3    DUMMY_SEED=1 SCORE=k2   N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s1_k2_b3_p2    DUMMY_SEED=1 SCORE=k2   N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s1_k2_b3_p3    DUMMY_SEED=1 SCORE=k2   N_BINS=3 MAX_PARENTS=3 ./run_all.sh

# --- データ乱数シード 2 (真の DAG を引き直した反復) ---
RUNDIR=./run_s2_bic_b2_p2   DUMMY_SEED=2 SCORE=bic  N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s2_bic_b2_p3   DUMMY_SEED=2 SCORE=bic  N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s2_bic_b3_p2   DUMMY_SEED=2 SCORE=bic  N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s2_bic_b3_p3   DUMMY_SEED=2 SCORE=bic  N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s2_bdeu_b2_p2  DUMMY_SEED=2 SCORE=bdeu N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s2_bdeu_b2_p3  DUMMY_SEED=2 SCORE=bdeu N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s2_bdeu_b3_p2  DUMMY_SEED=2 SCORE=bdeu N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s2_bdeu_b3_p3  DUMMY_SEED=2 SCORE=bdeu N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s2_k2_b2_p2    DUMMY_SEED=2 SCORE=k2   N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s2_k2_b2_p3    DUMMY_SEED=2 SCORE=k2   N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s2_k2_b3_p2    DUMMY_SEED=2 SCORE=k2   N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s2_k2_b3_p3    DUMMY_SEED=2 SCORE=k2   N_BINS=3 MAX_PARENTS=3 ./run_all.sh

# --- データ乱数シード 3 (真の DAG を引き直した反復) ---
RUNDIR=./run_s3_bic_b2_p2   DUMMY_SEED=3 SCORE=bic  N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s3_bic_b2_p3   DUMMY_SEED=3 SCORE=bic  N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s3_bic_b3_p2   DUMMY_SEED=3 SCORE=bic  N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s3_bic_b3_p3   DUMMY_SEED=3 SCORE=bic  N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s3_bdeu_b2_p2  DUMMY_SEED=3 SCORE=bdeu N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s3_bdeu_b2_p3  DUMMY_SEED=3 SCORE=bdeu N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s3_bdeu_b3_p2  DUMMY_SEED=3 SCORE=bdeu N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s3_bdeu_b3_p3  DUMMY_SEED=3 SCORE=bdeu N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s3_k2_b2_p2    DUMMY_SEED=3 SCORE=k2   N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s3_k2_b2_p3    DUMMY_SEED=3 SCORE=k2   N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s3_k2_b3_p2    DUMMY_SEED=3 SCORE=k2   N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s3_k2_b3_p3    DUMMY_SEED=3 SCORE=k2   N_BINS=3 MAX_PARENTS=3 ./run_all.sh

# --- データ乱数シード 4 (真の DAG を引き直した反復) ---
RUNDIR=./run_s4_bic_b2_p2   DUMMY_SEED=4 SCORE=bic  N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s4_bic_b2_p3   DUMMY_SEED=4 SCORE=bic  N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s4_bic_b3_p2   DUMMY_SEED=4 SCORE=bic  N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s4_bic_b3_p3   DUMMY_SEED=4 SCORE=bic  N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s4_bdeu_b2_p2  DUMMY_SEED=4 SCORE=bdeu N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s4_bdeu_b2_p3  DUMMY_SEED=4 SCORE=bdeu N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s4_bdeu_b3_p2  DUMMY_SEED=4 SCORE=bdeu N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s4_bdeu_b3_p3  DUMMY_SEED=4 SCORE=bdeu N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s4_k2_b2_p2    DUMMY_SEED=4 SCORE=k2   N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s4_k2_b2_p3    DUMMY_SEED=4 SCORE=k2   N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s4_k2_b3_p2    DUMMY_SEED=4 SCORE=k2   N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s4_k2_b3_p3    DUMMY_SEED=4 SCORE=k2   N_BINS=3 MAX_PARENTS=3 ./run_all.sh

# --- データ乱数シード 5 (真の DAG を引き直した反復) ---
RUNDIR=./run_s5_bic_b2_p2   DUMMY_SEED=5 SCORE=bic  N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s5_bic_b2_p3   DUMMY_SEED=5 SCORE=bic  N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s5_bic_b3_p2   DUMMY_SEED=5 SCORE=bic  N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s5_bic_b3_p3   DUMMY_SEED=5 SCORE=bic  N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s5_bdeu_b2_p2  DUMMY_SEED=5 SCORE=bdeu N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s5_bdeu_b2_p3  DUMMY_SEED=5 SCORE=bdeu N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s5_bdeu_b3_p2  DUMMY_SEED=5 SCORE=bdeu N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s5_bdeu_b3_p3  DUMMY_SEED=5 SCORE=bdeu N_BINS=3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s5_k2_b2_p2    DUMMY_SEED=5 SCORE=k2   N_BINS=2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s5_k2_b2_p3    DUMMY_SEED=5 SCORE=k2   N_BINS=2 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./run_s5_k2_b3_p2    DUMMY_SEED=5 SCORE=k2   N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s5_k2_b3_p3    DUMMY_SEED=5 SCORE=k2   N_BINS=3 MAX_PARENTS=3 ./run_all.sh


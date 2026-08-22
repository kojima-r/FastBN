#!/usr/bin/env bash
# =============================================================================
# bootstrap.sh — example_sc の全実験 (DATASET x DISC の 4 通り, 全遺伝子)
#   nohup ./bootstrap.sh > bootstrap.log 2>&1 &
#   出力は run_<DATASET>_<DISC>_all/ に分かれる。
#
# -----------------------------------------------------------------------------
# ITERS / ITERS_BS の決め方 — 変数数 D とサンプル数 N から決める
# -----------------------------------------------------------------------------
# (1) Hill-Climb は成長段階では 1 反復に 1 辺しか足せない。したがって必要な
#     反復数の下限は「学習される辺の本数 |E|」そのものである。
# (2) |E| <= D x P_eff。P_eff は 1 ノードが実際に持てる親の数で
#         P_eff = min(MAX_PARENTS, floor(log_r(N / 10)))
#     r は 1 変数の状態数 (= 離散化の段階数)。親を P 個持つノードの CPT は
#     r^P 通りの親設定を持ち、1 設定あたり 10 サンプル程度は無いと罰則付き
#     スコア (BDeu / BIC / K2) はその親を保持しない。つまり P_eff は N が決める。
# (3) 成長後の修正 (REMOVE / REVERSE) と Tabu が局所最適から抜ける分に、
#     成長段階の 2 倍を見込む:   ITERS = 3 x D x P_eff
# (4) ブートストラップは out/edges.tsv からの warm start なので成長段階が要らず、
#     リサンプル間の差分を直せればよい:   ITERS_BS = D x P_eff  (= ITERS の 1/3)
# (5) リサンプル総数 B = BOOTSTRAP x SEEDS は、エッジ出現確率の標準誤差
#     sqrt(p(1-p)/B) <= 0.5/sqrt(B) で決める。B = 5 x 20 = 100 で SE <= 0.05 なので、
#     THRESHOLD_PROB=0.3 の採否は十分な余裕をもって判定できる。
#     (MAX_JOBS は結果に影響しない実行上の都合。CPU コア数とメモリで決める)
#
#   これを 4 ケースに当てはめた値:
#
#     ケース      D     N    r  P_eff   ITERS=3D*P_eff   ITERS_BS=D*P_eff
#     bbknn/bin  2488   240  2    3         22392  ->22000     7464 -> 7500
#     bbknn/tri  2488   240  3    2         14928  ->15000     4976 -> 5000
#     ss/bin     6862  4500  2    3         61758  ->62000    20586 ->20000
#     ss/tri     6862  4500  3    3         61758  ->62000    20586 ->20000
#
#   bbknn/tri だけ反復が少ないのはサンプル数のため: N=240 で 3 値だと親 3 個は
#   3^3 = 27 通りの親設定に対し 1 設定あたり 8.9 サンプルしかなく、スコアが 3 個目の
#   親を保持しない (2 値なら 2^3 = 8 通りで 30 サンプル/設定あり、親 3 個を支えられる)。
#   逆に ss は N=4500 と多いので、離散化段階によらず MAX_PARENTS=3 を支えられる。
#
#   注意: ss の 2 行は上の規則どおりの反復数で、遺伝子数が多いぶん非常に長時間
#   (数週間規模) になる。1 反復のコストは辺が増えるほど上がるためで、途中で
#   打ち切って構わない場合や先に骨格だけ見たい場合は、代わりに末尾のコメント行
#   (ITERS=6862 = D、平均入次数 1 まで成長させる縮小版) を使う。
# =============================================================================

DATASET=bbknn DISC=bin ITERS=22000 ITERS_BS=7500  BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh
DATASET=bbknn DISC=tri ITERS=15000 ITERS_BS=5000  BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh
DATASET=ss    DISC=bin ITERS=62000 ITERS_BS=20000 BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh
DATASET=ss    DISC=tri ITERS=62000 ITERS_BS=20000 BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh

# ss の縮小版 (ITERS を D = 6862 = 平均入次数 1 にとどめる。ITERS_BS はその 1/3。
#  規則どおりの値では回しきれない場合に、上の ss 2 行と入れ替えて使う)
#DATASET=ss   DISC=bin ITERS=6862  ITERS_BS=2300  BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh
#DATASET=ss   DISC=tri ITERS=6862  ITERS_BS=2300  BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh

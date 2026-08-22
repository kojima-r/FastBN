#!/usr/bin/env bash
# =============================================================================
# config.sh — example_bnlearn の全設定
#   各ステップスクリプト (00〜04) は先頭でこのファイルを source する。
#   すべて環境変数で上書きできる (このファイルを編集しなくてよい)。
#     NETWORKS=asia SAMPLE_SIZES=1000 ./run_all.sh
#     SCORES=bic REPLICATES=1 ./02learn.sh
# =============================================================================

# --- ベンチマークの対象 -------------------------------------------------------
# bnlearn Bayesian Network Repository の discrete-small (ノード数 20 未満) 全 5 種
#   asia(8ノード/8辺) cancer(5/4) earthquake(5/4) sachs(11/17) survey(6/6)
export NETWORKS="${NETWORKS:-asia cancer earthquake sachs survey}"

# 各ネットワークから生成するサンプル数
export SAMPLE_SIZES="${SAMPLE_SIZES:-100 500 1000 5000}"

# サンプル数ごとの繰り返し回数 (乱数シード 1..REPLICATES)
export REPLICATES="${REPLICATES:-5}"

# 比較するスコア関数
export SCORES="${SCORES:-bic bdeu k2}"

# 総実行回数 = |NETWORKS| x |SAMPLE_SIZES| x REPLICATES x |SCORES|
# 既定は 5 x 4 x 5 x 3 = 300 回 (数分)

# --- 出力先 -------------------------------------------------------------------
export NETDIR="${NETDIR:-./networks}"        # ダウンロードした .bif
export RUNDIR="${RUNDIR:-./results}"
export DATADIR="${DATADIR:-${RUNDIR}/data}"  # 生成したサンプルと正解エッジ
export OUTROOT="${OUTROOT:-${RUNDIR}/out}"   # 学習結果 (実行ごとにサブディレクトリ)
export EVALDIR="${EVALDIR:-${RUNDIR}/eval}"  # 実行ごとの評価
export BENCHMARK="${BENCHMARK:-${RUNDIR}/benchmark.tsv}"   # 全実行を 1 行ずつ
export SUMMARY_TSV="${SUMMARY_TSV:-${RUNDIR}/summary.tsv}"
export SUMMARY_MD="${SUMMARY_MD:-${RUNDIR}/summary.md}"
export SUMMARY_PNG="${SUMMARY_PNG:-${RUNDIR}/summary.png}"
export SUMMARY_OVERALL_TSV="${SUMMARY_OVERALL_TSV:-${RUNDIR}/summary_overall.tsv}"
export FIGDIR="${FIGDIR:-${RUNDIR}/figures}"               # 正解 vs 学習の比較図
export REPORT_HTML="${REPORT_HTML:-${RUNDIR}/report.html}"

# --- 構造学習のパラメータ (../script/learn_structure.sh に渡る) --------------
export MAX_PARENTS="${MAX_PARENTS:-3}"   # 真の最大入次数は 3 (sachs)
export ITERS="${ITERS:-2000}"            # 小規模なので十分収束する
export TABU="${TABU:-10}"
export TOPK="${TOPK:-20}"                # ノード数より大きいので実質無制限
export ESS="${ESS:-1}"                   # BDeu の等価サンプルサイズ (bnlearn 既定と同じ)
export ALPHA="${ALPHA:-1.0}"
export CAND_METRIC="${CAND_METRIC:-mi}"
export JINDEX_CACHE="${JINDEX_CACHE:-1024}"
export REACH="${REACH:-lazy}"

# --- 評価のパラメータ ---------------------------------------------------------
# KL 用に学習構造の CPT を推定するときの Dirichlet 平滑化。
# 0 にすると未観測の親設定で確率 0 が出て KL が発散しうる。
export KL_ALPHA="${KL_ALPHA:-1.0}"
# KL を厳密計算する状態空間の上限 (sachs で 3^11 = 177147)
export MAX_STATES="${MAX_STATES:-2000000}"

# 表に載せる指標 (Precision / Recall / F1 を directed・skeleton の両方について出す)
export SUMMARY_METRICS="${SUMMARY_METRICS:-shd,precision_directed,recall_directed,f1_directed,precision_skeleton,recall_skeleton,f1_skeleton,sid_normalized,kl_divergence}"
# グラフに描く指標 (多すぎると読めないので絞る)
export SUMMARY_PLOT_METRICS="${SUMMARY_PLOT_METRICS:-shd,f1_directed,f1_skeleton,sid_normalized,kl_divergence}"
export SUMMARY_GROUP_BY="${SUMMARY_GROUP_BY:-network,n,score}"

# --- 比較図 (05compare.sh) ----------------------------------------------------
# 既定では最大サンプル数・反復 1・全スコアぶんだけ描く。
#   VIZ_N="100 5000"  : 描くサンプル数 ("all" で SAMPLE_SIZES 全部)
#   VIZ_SCORES=bic    : 描くスコア
#   VIZ_REP=1         : 描く反復
export VIZ_N="${VIZ_N:-}"
export VIZ_SCORES="${VIZ_SCORES:-}"
export VIZ_REP="${VIZ_REP:-1}"

# --- データ取得 ---------------------------------------------------------------
export BNLEARN_BASE_URL="${BNLEARN_BASE_URL:-https://www.bnlearn.com/bnrepository}"

# 汎用スクリプトの場所 (このディレクトリからの相対)
export BN_SCRIPTS="${BN_SCRIPTS:-../script}"

#!/usr/bin/env bash
# =============================================================================
# config.sh — example_dream の全設定
#   各ステップスクリプト (00〜06) は先頭でこのファイルを source する。
#   すべて環境変数で上書きできる。
#     DATASETS=dream4 D4_SIZES=10 ./run_all.sh
# =============================================================================

# --- 対象データセット ---------------------------------------------------------
#   dream4 : DREAM4 in silico network challenge (人工ネットワーク 10 個)
#   dream5 : DREAM5 network inference challenge (in silico / E.coli / S.cerevisiae)
#   hpn    : HPN-DREAM (Synapse 認証が要るため、手動配置した場合のみ)
export DATASETS="${DATASETS:-dream4 dream5}"

# DREAM4: 対象サイズと結合する実験の種類
export D4_SIZES="${D4_SIZES:-10 100}"
export D4_PARTS="${D4_PARTS:-multifactorial knockouts knockdowns timeseries wildtype}"

# DREAM5: 対象ネットワーク (2 = S.aureus は gold standard が採点に未使用のため除く)
export D5_NETWORKS="${D5_NETWORKS:-1 3 4}"
# DREAM5 / HPN は数千遺伝子ある。既定は 0 = 全遺伝子を使う (計算時間は下の表を参照)。
# 手早く試したいときは D5_MAX_VARS=300 のように絞る (TF を優先し残りは分散上位)。
export D5_MAX_VARS="${D5_MAX_VARS:-0}"

# 離散化
export BINS="${BINS:-3}"
export DISC_METHOD="${DISC_METHOD:-quantile}"

# 比較するスコア関数
export SCORES="${SCORES:-bic bdeu k2}"

# --- 出力先 -------------------------------------------------------------------
export SRCDIR="${SRCDIR:-./source}"           # ダウンロードした zip
export HPN_SRCDIR="${HPN_SRCDIR:-./source/hpn}"  # HPN-DREAM を手動配置する場所
export RUNDIR="${RUNDIR:-./results}"
export DATADIR="${DATADIR:-${RUNDIR}/data}"
export TRUTHDIR="${TRUTHDIR:-${RUNDIR}/truth}"
export OUTROOT="${OUTROOT:-${RUNDIR}/out}"
export EVALDIR="${EVALDIR:-${RUNDIR}/eval}"
export NETLIST="${NETLIST:-${RUNDIR}/networks.txt}"
export BENCHMARK="${BENCHMARK:-${RUNDIR}/benchmark.tsv}"
export SUMMARY_TSV="${SUMMARY_TSV:-${RUNDIR}/summary.tsv}"
export SUMMARY_MD="${SUMMARY_MD:-${RUNDIR}/summary.md}"
export SUMMARY_PNG="${SUMMARY_PNG:-${RUNDIR}/summary.png}"
export SUMMARY_OVERALL_TSV="${SUMMARY_OVERALL_TSV:-${RUNDIR}/summary_overall.tsv}"
export FIGDIR="${FIGDIR:-${RUNDIR}/figures}"
export REPORT_HTML="${REPORT_HTML:-${RUNDIR}/report.html}"

# --- 構造学習のパラメータ -----------------------------------------------------
export MAX_PARENTS="${MAX_PARENTS:-3}"

# 反復数は**変数数に応じて自動スケール**する:
#     ITERS = max(ITERS_MIN, ITERS_PER_VAR x 変数数)
# Hill-Climb は成長段階で 1 反復あたり 1 辺しか足せないため、固定の反復数だと
# 大きい網では「反復上限で打ち切られた構造」を評価してしまう (DREAM5 の全遺伝子で
# 実際にそうなる: 正解 4012 辺に対し ITERS=2000 では最大 2000 辺しか置けない)。
# ITERS を明示すると自動スケールを上書きして固定値になる。
export ITERS_PER_VAR="${ITERS_PER_VAR:-10}"
export ITERS_MIN="${ITERS_MIN:-2000}"
export TABU="${TABU:-10}"
export TOPK="${TOPK:-20}"
export ESS="${ESS:-1}"
export ALPHA="${ALPHA:-1.0}"
export CAND_METRIC="${CAND_METRIC:-mi}"
export JINDEX_CACHE="${JINDEX_CACHE:-1024}"
export REACH="${REACH:-lazy}"

# --- 評価 ---------------------------------------------------------------------
# 正解は遺伝子制御ネットワークでフィードバックを含むため SID は多くの場合 NA。
# KL は真の CPT が無いので計算しない。
export MAX_SID_NODES="${MAX_SID_NODES:-150}"   # これより大きいと SID を省略
export SUMMARY_METRICS="${SUMMARY_METRICS:-shd,precision_directed,recall_directed,f1_directed,precision_skeleton,recall_skeleton,f1_skeleton}"
export SUMMARY_PLOT_METRICS="${SUMMARY_PLOT_METRICS:-shd,f1_directed,f1_skeleton}"
export SUMMARY_GROUP_BY="${SUMMARY_GROUP_BY:-dataset,network,score}"

# --- 比較図 -------------------------------------------------------------------
# ノードが多い図は読めないので、これを超える網は次数上位のハブ部分グラフを描く
export VIZ_MAX_NODES="${VIZ_MAX_NODES:-100}"
export VIZ_SCORES="${VIZ_SCORES:-}"

# --- データ取得 ---------------------------------------------------------------
export DREAM4_URL="${DREAM4_URL:-https://gnw.sourceforge.net/resources/DREAM4%20in%20silico%20challenge.zip}"
export DREAM5_URL="${DREAM5_URL:-https://zenodo.org/records/17854236/files/1_Challenge_Data_Supplement.zip?download=1}"

export BN_SCRIPTS="${BN_SCRIPTS:-../script}"

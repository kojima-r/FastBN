#!/usr/bin/env bash
# =============================================================================
# config.sh — example_sachs の全設定
#   各ステップスクリプト (00〜06) は先頭でこのファイルを source する。
#   すべて環境変数で上書きできる。
#     PRESETS=obs BINS=3 ./run_all.sh
#     SCORES=bic ./02learn.sh
# =============================================================================

# --- ベンチマークの対象 -------------------------------------------------------
# Sachs のフローサイトメトリー・データ (11 タンパク質 x 14 実験条件)。
#   obs : 一般刺激のみ (cd3cd28; 853 細胞)  = 介入なしの古典的な設定
#   int : 阻害剤などの介入条件のみ          (10818 細胞)
#   all : 14 条件すべて                     (11671 細胞)
export PRESETS="${PRESETS:-obs int all}"

# 離散化の段階数 (Sachs らの原論文は 3 段階)
export BINS_LIST="${BINS_LIST:-2 3}"

# 比較するスコア関数
export SCORES="${SCORES:-bic bdeu k2}"

# 離散化の方法と前処理
export DISC_METHOD="${DISC_METHOD:-quantile}"   # quantile | uniform
export USE_LOG2="${USE_LOG2:-1}"                # 1 で log2(x+1) 変換

# 総実行回数 = |PRESETS| x |BINS_LIST| x |SCORES| (既定 3 x 2 x 3 = 18)

# --- 出力先 -------------------------------------------------------------------
export SRCDIR="${SRCDIR:-./source}"          # 展開した sachs.zip
export RUNDIR="${RUNDIR:-./results}"
export DATADIR="${DATADIR:-${RUNDIR}/data}"
export OUTROOT="${OUTROOT:-${RUNDIR}/out}"
export EVALDIR="${EVALDIR:-${RUNDIR}/eval}"
export BENCHMARK="${BENCHMARK:-${RUNDIR}/benchmark.tsv}"
export SUMMARY_TSV="${SUMMARY_TSV:-${RUNDIR}/summary.tsv}"
export SUMMARY_MD="${SUMMARY_MD:-${RUNDIR}/summary.md}"
export SUMMARY_PNG="${SUMMARY_PNG:-${RUNDIR}/summary.png}"
export SUMMARY_OVERALL_TSV="${SUMMARY_OVERALL_TSV:-${RUNDIR}/summary_overall.tsv}"
export FIGDIR="${FIGDIR:-${RUNDIR}/figures}"
export REPORT_HTML="${REPORT_HTML:-${RUNDIR}/report.html}"
export TRUE_EDGES="${TRUE_EDGES:-${DATADIR}/true_edges.tsv}"

# --- 構造学習のパラメータ -----------------------------------------------------
export MAX_PARENTS="${MAX_PARENTS:-3}"   # 正解の最大入次数は 3
export ITERS="${ITERS:-3000}"
export TABU="${TABU:-10}"
export TOPK="${TOPK:-20}"                # ノード数 11 より大きいので実質無制限
export ESS="${ESS:-1}"
export ALPHA="${ALPHA:-1.0}"
export CAND_METRIC="${CAND_METRIC:-mi}"
export JINDEX_CACHE="${JINDEX_CACHE:-1024}"
export REACH="${REACH:-lazy}"

# --- 評価 ---------------------------------------------------------------------
# 正解構造 (GroundTruth.csv) は PKA <-> PIP3 の相互作用を含み DAG ではないため、
# SID は自動的に NA になる。KL は真の CPT が無いので計算しない。
export SUMMARY_METRICS="${SUMMARY_METRICS:-shd,precision_directed,recall_directed,f1_directed,precision_skeleton,recall_skeleton,f1_skeleton}"
export SUMMARY_PLOT_METRICS="${SUMMARY_PLOT_METRICS:-shd,f1_directed,f1_skeleton}"
export SUMMARY_GROUP_BY="${SUMMARY_GROUP_BY:-preset,bins,score}"

# --- 比較図 -------------------------------------------------------------------
export VIZ_PRESETS="${VIZ_PRESETS:-}"    # 既定は PRESETS 全部
export VIZ_BINS="${VIZ_BINS:-3}"
export VIZ_SCORES="${VIZ_SCORES:-}"      # 既定は SCORES 全部

# --- データ取得 ---------------------------------------------------------------
export SACHS_URL="${SACHS_URL:-https://zenodo.org/records/7681811/files/sachs.zip?download=1}"

# 汎用スクリプトの場所
export BN_SCRIPTS="${BN_SCRIPTS:-../script}"

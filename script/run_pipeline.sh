#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh
#   バルク RNA 発現量データに対する解析を一括実行する汎用ドライバ。
#   各ステップは独立スクリプトなので、個別に実行しても構わない。
#
#   1) preprocess.sh           前処理 (正規化 -> log -> フィルタ -> 離散化)
#   2) learn_structure.sh      構造学習 (Hill-Climb + Tabu)
#   3) edge_importance.sh      エッジ重要度 (全サンプル)
#   4) bootstrap_stability.sh  ブートストラップ -> コンセンサス網 -> 重要度
#   5) importance_groups.sh    群別エッジ重要度 (学習網 / コンセンサス網)
#   6) viz*.sh                 可視化
#   7) make_report.sh          HTML レポート
#
# 使い方 (解析ディレクトリで実行):
#   source ./config.sh && ../script/run_pipeline.sh
#   DO_BOOTSTRAP=0 DO_GROUPS=0 ../script/run_pipeline.sh   # 学習+重要度のみ
#
# ステップの ON/OFF (既定はすべて 1 = 実行):
#   DO_PREPROCESS / DO_LEARN / DO_IMPORTANCE / DO_BOOTSTRAP /
#   DO_GROUPS / DO_VIZ / DO_REPORT
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="pipeline"; bn_tag="pipeline"

S="${BN_SCRIPT_DIR}"

do_preprocess="${DO_PREPROCESS:-1}"
do_learn="${DO_LEARN:-1}"
do_importance="${DO_IMPORTANCE:-1}"
do_bootstrap="${DO_BOOTSTRAP:-1}"
do_groups="${DO_GROUPS:-1}"
do_viz="${DO_VIZ:-1}"
do_report="${DO_REPORT:-1}"

step=0
next() { step=$(( step + 1 )); hr "STEP ${step}: $*"; }

if [[ "${do_preprocess}" != "0" ]]; then
  next "前処理 (preprocess.sh)"
  "${S}/preprocess.sh"
fi

if [[ "${do_learn}" != "0" ]]; then
  next "構造学習 (learn_structure.sh)"
  "${S}/learn_structure.sh"
fi

if [[ "${do_importance}" != "0" ]]; then
  next "エッジ重要度 (edge_importance.sh)"
  "${S}/edge_importance.sh"
fi

if [[ "${do_bootstrap}" != "0" ]]; then
  next "ブートストラップ安定性解析 (bootstrap_stability.sh)"
  "${S}/bootstrap_stability.sh"
fi

if [[ "${do_groups}" != "0" ]]; then
  next "群別エッジ重要度 (importance_groups.sh)"
  # (a) 学習網ベース
  "${S}/importance_groups.sh"
  # (b) コンセンサス網ベース (ブートストラップ済みの場合のみ)
  if [[ -f "${OUTDIR}/integ_edges.tsv" && -f "${OUTDIR}/integ_all_counts.tsv" ]]; then
    log "コンセンサス網ベースの群別重要度も計算します"
    INIT="${OUTDIR}/integ_edges.tsv" \
    COUNTS="${OUTDIR}/integ_all_counts.tsv" \
    REF_EDGES="${OUTDIR}/integ_edges2.tsv" \
    REF_NAMED="${OUTDIR}/integ_edges_named.tsv" \
    OUT_PREFIX=integ_edge_importance \
      "${S}/importance_groups.sh"
  fi
fi

if [[ "${do_viz}" != "0" ]]; then
  next "可視化 (viz*.sh)"
  "${S}/viz.sh"
  if [[ -f "${OUTDIR}/integ_edges2.tsv" ]]; then
    "${S}/viz_bs.sh"
  fi
  # 群別図は群別重要度がある場合のみ
  if compgen -G "${OUTDIR}/edge_importance_g*_*.tsv" > /dev/null; then
    "${S}/viz_subsets.sh"
  fi
  if compgen -G "${OUTDIR}/integ_edge_importance_g*_*.tsv" > /dev/null; then
    "${S}/viz_bs_subsets.sh"
  fi
fi

if [[ "${do_report}" != "0" ]]; then
  next "HTML レポート (make_report.sh)"
  "${S}/make_report.sh"
fi

hr "パイプライン完了"
show() { printf "   %-36s : %s\n" "$1" "$2"; }
show "${INPUT}" "離散化済み入力"
show "${OUTDIR}/edges_named.tsv" "学習ネットワーク (遺伝子名)"
show "${OUTDIR}/edge_importance.tsv" "エッジ重要度"
[[ -f "${OUTDIR}/integ_edges_named.tsv" ]] && \
  show "${OUTDIR}/integ_edges_named.tsv" "コンセンサス網 (遺伝子名)"
[[ -f "./report.html" ]] && show "./report.html" "HTML レポート"
echo "=================================================================="

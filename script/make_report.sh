#!/usr/bin/env bash
# =============================================================================
# make_report.sh
#   解析成果物 (図・重要度テーブル・データ要約) を 1 つの HTML レポートに
#   まとめる簡易呼び出しスクリプト (汎用版)。既定では図を相対リンクで参照する
#   (軽量・即開ける)。単一ファイルで共有したい場合は --embed。
#
#   前提: 先に図を生成しておくこと。
#     ../script/viz.sh              # 学習網        -> ./figures/
#     ../script/viz_bs.sh           # コンセンサス網 -> ./figures_bs/
#     ../script/viz_subsets.sh      # 群別 (任意)
#     ../script/viz_bs_subsets.sh   # 群別 (任意)
#
#   追加オプションはそのまま make_report.py へ渡される。
#     ../script/make_report.sh --embed                # base64 埋め込み
#     ../script/make_report.sh --metrics dlogL,dBIC   # メトリクスを限定
#     ../script/make_report.sh --top-edges 60 --title "私の解析"
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="report"; bn_tag="report"

out_html="${REPORT_HTML:-./report.html}"

"${PYTHON_BIN}" "$(py_tool make_report.py)" \
  --base-dir . \
  --out "${out_html}" \
  --figures "${FIGDIR:-./figures}" \
  --figures-bs "${FIGDIR_BS:-./figures_bs}" \
  --input-tsv "${INPUT}" \
  --var-map "${VARMAP}" \
  --samples "${SAMPLES}" \
  --groups-manifest "${GROUPDIR}/groups_manifest.tsv" \
  --target-file "${TARGET_FILE}" \
  --edges "${OUTDIR}/edges.tsv" \
  --integ-edges "${OUTDIR}/integ_edges2.tsv" \
  ${REPORT_TITLE:+--title "${REPORT_TITLE}"} \
  "$@"

log "HTML レポートを出力しました: ${out_html}"

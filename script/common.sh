#!/usr/bin/env bash
# =============================================================================
# common.sh
#   script/ 以下の汎用パイプラインスクリプトが共通で読み込む設定・ユーティリティ。
#   各スクリプトの先頭で
#       source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
#   のように読み込む。
#
# 設計方針:
#   * スクリプトは「カレントディレクトリ = 解析ディレクトリ」として動作する。
#     (gssg_analysis のようにスクリプト自身の場所へ cd はしない)
#     -> 任意のデータセット用ディレクトリを作り、そこから ../script/*.sh を
#        呼ぶだけで解析できる。
#   * 設定は全て環境変数で上書きする。解析ディレクトリに config.sh を置き、
#     `source ./config.sh` してから呼ぶ運用を推奨 (example_bulk 参照)。
# =============================================================================

# --- 場所の解決 --------------------------------------------------------------
# BN_SCRIPT_DIR : この script/ ディレクトリ (python ツールの置き場)
# BN_REPO_DIR   : FastBN リポジトリのルート (fast_bn バイナリの既定の置き場)
BN_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BN_REPO_DIR="$(cd "${BN_SCRIPT_DIR}/.." && pwd)"
export BN_SCRIPT_DIR BN_REPO_DIR

# fast_bn バイナリ (FASTBN_BIN で明示指定可能)
FASTBN_BIN="${FASTBN_BIN:-${BN_REPO_DIR}/fast_bn}"
# ブートストラップ統合スクリプト (リポジトリ直下)
BS_PROB_PY="${BS_PROB_PY:-${BN_REPO_DIR}/compute_bs_prob.py}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

# --- 既定のディレクトリ / ファイル配置 --------------------------------------
# 解析ディレクトリ (カレント) 以下のレイアウト。全て環境変数で上書き可能。
DATADIR="${DATADIR:-./data}"          # 前処理の入出力
OUTDIR="${OUTDIR:-./out}"            # 構造学習・重要度の出力
BSDIR="${BSDIR:-./bs}"               # ブートストラップの生出力
GROUPDIR="${GROUPDIR:-./groups}"     # 群別 score-dataset
INPUT="${INPUT:-${DATADIR}/expr_disc.tsv}"      # fast_bn 入力 (離散化済み)
VARMAP="${VARMAP:-${DATADIR}/var_map.tsv}"      # 列インデックス <-> 遺伝子対応表
SAMPLES="${SAMPLES:-${DATADIR}/samples.tsv}"    # サンプル順と群ラベル
TARGET_FILE="${TARGET_FILE:-./target_genes.txt}" # 注目遺伝子 (任意)

# --- ログ --------------------------------------------------------------------
bn_tag="${BN_TAG:-fastbn}"

log()  { echo "[${bn_tag}] $*"; }
warn() { echo "[${bn_tag}] 警告: $*" >&2; }
die()  { echo "[${bn_tag}] エラー: $*" >&2; exit 1; }

hr() {
  echo "=================================================================="
  [[ $# -gt 0 ]] && echo " $*"
  [[ $# -gt 0 ]] && echo "=================================================================="
  return 0
}

# --- 前提チェック ------------------------------------------------------------
require_bin() {
  [[ -x "${FASTBN_BIN}" ]] || die "fast_bn バイナリが見つかりません: ${FASTBN_BIN}
  リポジトリ直下で ./compile.sh を実行してビルドするか、FASTBN_BIN で場所を指定してください。"
}

require_file() {
  local f
  for f in "$@"; do
    [[ -f "${f}" ]] || die "必須ファイルがありません: ${f}"
  done
}

# python ツールのパスを返す (例: py_tool preprocess_expr.py)
py_tool() { echo "${BN_SCRIPT_DIR}/$1"; }

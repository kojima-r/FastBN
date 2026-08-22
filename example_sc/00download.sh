#!/usr/bin/env bash
# =============================================================================
# 00download.sh — ステップ 0: データの取得
#   離散化済みの単一細胞発現量データ (Tabula Muris Senis 由来) を取得する。
#   既に展開済みのディレクトリがある場合は何もしない。
#
# 取得先:
#   ${DATA_BASE_URL}/<ディレクトリ名>.tar.gz
#   (config.sh の DATA_BASE_URL で変更できる。手元にアーカイブがある場合は
#    このディレクトリで tar xf するだけでよい)
#
# 展開されるディレクトリ (config.sh の DATASET で選ぶ):
#   data_bbknn_r_tissue_disc/  : 複数プロトコル統合 (BBKNN 補正後)
#   data_ss_r_tissue_disc/     : Smart-seq2 のみ
#
# それぞれの中身:
#   all_disc.tsv / all_disc10.tsv / all_disc100.tsv / all_disc1000.tsv
#       2 値離散化した行列 (行=サンプル, 列=遺伝子)。数字は先頭 N 遺伝子。
#   all_disc_tri*.tsv
#       同じものを 3 値離散化したもの。
#   tissue/ , tissue_tri/
#       組織ごとに分けた同じデータ (群別解析の群ラベルの元になる)。
#
# 使い方:
#   ./00download.sh          # config.sh の DATASET のぶんだけ取得
#   ./00download.sh --all    # 両方のデータセットを取得
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

want_all=0
[[ "${1:-}" == "--all" ]] && want_all=1

if [[ "${want_all}" -eq 1 ]]; then
  dirs=(data_bbknn_r_tissue_disc data_ss_r_tissue_disc)
else
  dirs=("$(basename "${DATA_ROOT}")")
fi

fetch() {
  local name="$1"
  if [[ -d "./${name}" ]]; then
    echo "[00download] ./${name} は既にあります (スキップ)"
    return 0
  fi
  local url="${DATA_BASE_URL}/${name}.tar.gz"
  echo "[00download] 取得: ${url}"
  if command -v wget > /dev/null; then
    wget -q --show-progress "${url}" -O "${name}.tar.gz"
  elif command -v curl > /dev/null; then
    curl -fL --progress-bar "${url}" -o "${name}.tar.gz"
  else
    echo "[00download] エラー: wget も curl もありません。" >&2
    echo "             ${url} を手動で取得し、このディレクトリで展開してください。" >&2
    exit 1
  fi
  echo "[00download] 展開: ${name}.tar.gz"
  tar xf "${name}.tar.gz"
  rm -f "${name}.tar.gz"
  [[ -d "./${name}" ]] || {
    echo "[00download] エラー: 展開しても ./${name} ができませんでした。" >&2
    echo "             アーカイブのディレクトリ構成を確認してください。" >&2
    exit 1
  }
}

for d in "${dirs[@]}"; do
  fetch "${d}"
done

echo "[00download] 完了。現在のデータ:"
for d in data_bbknn_r_tissue_disc data_ss_r_tissue_disc; do
  [[ -d "./${d}" ]] && echo "   ./${d}/ ($(ls -1 "./${d}"/*.tsv 2>/dev/null | wc -l) 行列, 組織別 $(ls -1 "./${d}"/tissue/*.tsv 2>/dev/null | wc -l) ファイル)"
done
echo "[00download] 次に選択されている入力: ${SRC_MATRIX}"

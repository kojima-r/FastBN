#!/usr/bin/env bash
# =============================================================================
# 00download.sh — ステップ 0: 正解ネットワークの取得
#   bnlearn Bayesian Network Repository の discrete-small カテゴリ
#   (https://www.bnlearn.com/bnrepository/discrete-small.html) から
#   BIF 形式のネットワーク定義をダウンロードして展開する。
#   既に取得済みのものはスキップする。
#
# 取得先: ${BNLEARN_BASE_URL}/<net>/<net>.bif.gz
# 出力  : ${NETDIR}/<net>.bif
#
# ダウンロードされるのは「正解のベイジアンネットワーク (構造 + CPT)」であって
# データではない。データは 01sample.sh がこの CPT からサンプリングして作る。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${NETDIR}"

fetch() {
  local net="$1" out="${NETDIR}/$1.bif"
  if [[ -s "${out}" ]]; then
    echo "[00download] ${out} は取得済み (スキップ)"
    return 0
  fi
  local url="${BNLEARN_BASE_URL}/${net}/${net}.bif.gz"
  echo "[00download] 取得: ${url}"
  if command -v curl > /dev/null; then
    curl -fsSL "${url}" | gunzip > "${out}"
  elif command -v wget > /dev/null; then
    wget -qO- "${url}" | gunzip > "${out}"
  else
    echo "[00download] エラー: curl も wget もありません。" >&2
    echo "             ${url} を手動で取得し ${out} に展開してください。" >&2
    exit 1
  fi
  [[ -s "${out}" ]] || { echo "[00download] エラー: ${out} が空です" >&2; exit 1; }
}

for net in ${NETWORKS}; do
  fetch "${net}"
done

echo "[00download] 取得したネットワーク:"
for net in ${NETWORKS}; do
  python3 "${BN_SCRIPTS}/bif_io.py" info --bif "${NETDIR}/${net}.bif" \
    | awk -F'\t' -v n="${net}" '
        $1=="nodes"{nodes=$2} $1=="arcs"{arcs=$2}
        $1=="parameters"{par=$2} $1=="state_space"{ss=$2}
        END{printf "   %-12s ノード %3s / 辺 %3s / パラメータ %4s / 状態空間 %s\n",
                   n, nodes, arcs, par, ss}'
done

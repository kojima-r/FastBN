#!/usr/bin/env bash
# =============================================================================
# 01sample.sh — ステップ 1: 正解ネットワークからのデータ生成
#   各ネットワークの CPT から祖先サンプリングでデータを作る。
#   サンプル数 x 繰り返し (乱数シード) の組み合わせぶん生成する。
#
# 出力 (${DATADIR}):
#   <net>_n<N>_r<R>.tsv   : fast_bn 入力 (行=サンプル, 列=変数, 値=状態コード)
#   <net>_true_edges.tsv  : 正解エッジ (変数名, u -> v)
#
# 列の順序は BIF の変数宣言順で全ファイル共通。fast_bn のノード番号 = 列位置
# なので、評価時にこのヘッダで番号 -> 変数名を復元できる。
#
# 出現しない状態がある場合 (asia の稀な状態など) は警告が出る。fast_bn 側の
# 基数はデータの最大値+1 になるが、評価スクリプトは正解ネットワークの基数を
# 使うので KL は正しく計算される。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${DATADIR}"

n_files=0
for net in ${NETWORKS}; do
  bif="${NETDIR}/${net}.bif"
  [[ -s "${bif}" ]] || { echo "[01sample] エラー: ${bif} がありません。先に ./00download.sh" >&2; exit 1; }
  for n in ${SAMPLE_SIZES}; do
    for rep in $(seq 1 "${REPLICATES}"); do
      out="${DATADIR}/${net}_n${n}_r${rep}.tsv"
      edges_opt=()
      # 正解エッジはネットワークごとに 1 回だけ書けばよい
      [[ "${n}" == "$(echo ${SAMPLE_SIZES} | awk '{print $1}')" && "${rep}" == "1" ]] \
        && edges_opt=(--out-edges "${DATADIR}/${net}_true_edges.tsv")
      python3 "${BN_SCRIPTS}/bif_io.py" sample \
        --bif "${bif}" --n "${n}" --seed "${rep}" \
        --out "${out}" "${edges_opt[@]}"
      n_files=$(( n_files + 1 ))
    done
  done
done

echo "[01sample] 完了: ${n_files} 個のデータセットを ${DATADIR}/ に生成しました"

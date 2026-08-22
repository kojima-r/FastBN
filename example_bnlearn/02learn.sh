#!/usr/bin/env bash
# =============================================================================
# 02learn.sh — ステップ 2: 構造学習 (ベンチマーク本体)
#   生成した全データセット x 全スコアについて fast_bn の Hill-Climb + Tabu を
#   走らせる。1 実行につき ../script/learn_structure.sh を呼ぶだけ
#   (INPUT と OUTDIR を実行ごとに切り替える)。
#
# 出力 (${OUTROOT}/<net>_n<N>_r<R>_<score>/):
#   edges.tsv        : 学習された DAG (ノード = 列インデックス)
#   edges_named.tsv  : 同じエッジを変数名で (edges.tsv と行対応)
#   all_counts.tsv   : カウント表
#   log_learn.txt    : fast_bn のログ (実行時間もここに出る)
#
# 総実行時間もここに出す。1 実行あたりの時間は log_learn.txt の
# "total time" 行を参照。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

mkdir -p "${OUTROOT}"
start=$(date +%s)
count=0

for net in ${NETWORKS}; do
  for n in ${SAMPLE_SIZES}; do
    for rep in $(seq 1 "${REPLICATES}"); do
      data="${DATADIR}/${net}_n${n}_r${rep}.tsv"
      [[ -s "${data}" ]] || { echo "[02learn] エラー: ${data} がありません。先に ./01sample.sh" >&2; exit 1; }
      for score in ${SCORES}; do
        run="${net}_n${n}_r${rep}_${score}"
        INPUT="${data}" OUTDIR="${OUTROOT}/${run}" SCORE="${score}" \
          "${BN_SCRIPTS}/learn_structure.sh" > /dev/null 2>&1 \
          || { echo "[02learn] エラー: ${run} の学習に失敗しました" >&2
               INPUT="${data}" OUTDIR="${OUTROOT}/${run}" SCORE="${score}" \
                 "${BN_SCRIPTS}/learn_structure.sh" >&2 || true
               exit 1; }
        count=$(( count + 1 ))
        printf "\r[02learn] %d 実行完了 (最新: %s)          " "${count}" "${run}"
      done
    done
  done
done

echo
echo "[02learn] 完了: ${count} 実行 / $(( $(date +%s) - start )) 秒 -> ${OUTROOT}/"

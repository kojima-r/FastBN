#!/usr/bin/env bash
# =============================================================================
# run_all.sh — ステップ 0〜4 を順番に実行する
#   いずれかで失敗したら止まる (set -e)。
#
#   ./run_all.sh                                   # 既定 (300 実行, 数分)
#   NETWORKS=asia SAMPLE_SIZES=1000 ./run_all.sh   # 1 ネットワークだけ
#   SCORES=bic REPLICATES=2 ./run_all.sh           # 軽く試す
#   SKIP_DOWNLOAD=1 ./run_all.sh                   # 取得済みの .bif を使う
#
#   各ステップは単独でも実行できるので、学習パラメータを変えて試すときは
#   config.sh を編集 (または環境変数を指定) して 02 以降だけ再実行すればよい。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

start=$(date +%s)

n_runs=0
for net in ${NETWORKS}; do
  for n in ${SAMPLE_SIZES}; do
    for score in ${SCORES}; do
      n_runs=$(( n_runs + REPLICATES ))
    done
  done
done

echo "=================================================================="
echo " bnlearn discrete-small ベンチマーク"
echo "   ネットワーク : ${NETWORKS}"
echo "   サンプル数   : ${SAMPLE_SIZES}  (各 ${REPLICATES} 反復)"
echo "   スコア       : ${SCORES}"
echo "   総実行回数   : ${n_runs}"
echo "   出力         : ${RUNDIR}/"
echo "=================================================================="

if [[ "${SKIP_DOWNLOAD:-0}" == "0" ]]; then
  ./00download.sh
fi
./01sample.sh
./02learn.sh
./03evaluate.sh
./04summarize.sh
./05compare.sh
./06report.sh

echo "=================================================================="
echo " 完了 ($(( $(date +%s) - start )) 秒)"
echo "   ${REPORT_HTML}  : HTML レポート (表 + グラフ + 正解との比較図)"
echo "   ${SUMMARY_MD}   : 集約表 (Markdown)"
echo "   ${FIGDIR}/      : 正解 vs 学習ネットワークの比較図"
echo "   ${BENCHMARK}    : 実行ごとの生の指標"
echo "=================================================================="

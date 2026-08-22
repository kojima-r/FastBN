#!/usr/bin/env bash
# =============================================================================
# run_all.sh — ステップ 0〜8 を順番に実行する
#   いずれかで失敗したら止まる (set -e)。
#
#   ./run_all.sh                 # 全ステップ
#   SKIP_MAKE_DATA=1 ./run_all.sh  # 既存の data/counts.tsv を使う
#
#   各ステップは単独でも実行できるので、パラメータを変えて試すときは
#   config.sh を編集して該当ステップだけ再実行すればよい。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

start=$(date +%s)

if [[ "${SKIP_MAKE_DATA:-0}" == "0" ]]; then
  ./00make_data.sh
fi
./01preprocess.sh
./02learn.sh
./03importance.sh
./04bootstrap.sh
./05importance_groups.sh
./06visualize.sh
./07report.sh
./08evaluate.sh

echo "=================================================================="
echo " 全ステップ完了 ($(( $(date +%s) - start )) 秒)"
echo "   ${REPORT_HTML} をブラウザで開くと結果を一覧できます"
echo "=================================================================="

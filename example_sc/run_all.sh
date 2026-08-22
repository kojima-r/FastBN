#!/usr/bin/env bash
# =============================================================================
# run_all.sh — ステップ 0〜8 を順番に実行する
#   いずれかで失敗したら止まる (set -e)。
#
#   ./run_all.sh                       # 既定 (bbknn / bin / 全遺伝子; 2 時間程度)
#   NVARS=100 ./run_all.sh             # まず流れを確認するならこちら (数分)
#   DATASET=ss ./run_all.sh            # Smart-seq2 データに切り替え
#   DISC=tri NVARS=1000 ./run_all.sh   # 3 値離散化・1000 遺伝子
#   SKIP_DOWNLOAD=1 ./run_all.sh       # 既に展開済みのデータを使う
#
#   結果は組み合わせごとに別ディレクトリ (${RUNDIR}) に出るので、設定を変えて
#   何度実行しても互いに上書きしない。
#
#   各ステップは単独でも実行できるので、パラメータを変えて試すときは
#   config.sh を編集 (または環境変数を指定) して該当ステップだけ再実行すればよい。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

start=$(date +%s)

echo "=================================================================="
echo " DATASET=${DATASET}  DISC=${DISC}  NVARS=${NVARS}"
echo " 入力: ${SRC_MATRIX}"
echo " 出力: ${RUNDIR}/"
echo "=================================================================="

if [[ "${SKIP_DOWNLOAD:-0}" == "0" ]]; then
  ./00download.sh
fi
./01prepare.sh
./02learn.sh
./03importance.sh
./04bootstrap.sh
./05importance_groups.sh
./06visualize.sh
./07report.sh
./08score_check.sh

echo "=================================================================="
echo " 全ステップ完了 ($(( $(date +%s) - start )) 秒)"
echo "   ${RUNDIR}/report.html をブラウザで開くと結果を一覧できます"
echo "=================================================================="

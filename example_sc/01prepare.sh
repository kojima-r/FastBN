#!/usr/bin/env bash
# =============================================================================
# 01prepare.sh — ステップ 1: 入力の準備
#   この例題のデータは**すでに離散化済み**なので、example_bulk のような前処理
#   (正規化 -> log -> フィルタ -> 離散化) は行わない。代わりに、選択された行列を
#   ../script/ の汎用スクリプトがそのまま扱える形に整える。
#
#   選択は config.sh の 3 つのスイッチで決まる:
#     DATASET (bbknn | ss) / DISC (bin | tri) / NVARS (10 | 100 | 1000 | all)
#
# 出力 (${DATADIR} = ${RUNDIR}/data):
#   expr_disc.tsv : fast_bn 入力 (選択した行列へのリンク; 行=サンプル, 列=遺伝子)
#   var_map.tsv   : 列インデックス <-> 遺伝子名 / 分散 / 出現水準数
#   samples.tsv   : 行番号 <-> サンプル ID / 群ラベル (**群 = 組織**)
#   ${TARGET_FILE}: 注目遺伝子リスト (config.sh の TARGET_GENES を指定した場合)
#
# samples.tsv は tissue/ 以下の組織別ファイルから作る。all_disc*.tsv の行順が
# 「組織別ファイルをファイル名順に連結したもの」と一致することを検証したうえで
# 組織名を群ラベルにするので、05importance_groups.sh で「どのエッジがどの組織で
# 効いているか」を比較できる。
#
# 使い方:
#   ./01prepare.sh
#   DATASET=ss DISC=tri NVARS=1000 ./01prepare.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

if [[ ! -f "${SRC_MATRIX}" ]]; then
  echo "[01prepare] エラー: ${SRC_MATRIX} がありません。" >&2
  echo "            先に ./00download.sh を実行するか、DATASET / DISC / NVARS を" >&2
  echo "            見直してください。" >&2
  exit 1
fi

echo "[01prepare] DATASET=${DATASET} DISC=${DISC} NVARS=${NVARS}"
echo "[01prepare] 入力行列  : ${SRC_MATRIX}"
echo "[01prepare] 組織ファイル: ${SRC_TISSUE_DIR}"
echo "[01prepare] 出力先    : ${RUNDIR}/"

mkdir -p "${RUNDIR}" "${DATADIR}"

python3 ./prepare_data.py \
  --matrix "${SRC_MATRIX}" \
  --tissue-dir "${SRC_TISSUE_DIR}" \
  --outdir "${DATADIR}" \
  --targets "${TARGET_GENES}" \
  --out-targets "${TARGET_FILE}" \
  "$@"

n_samples=$(( $(wc -l < "${INPUT}") - 1 ))
n_vars=$(head -1 "${INPUT}" | awk -F'\t' '{print NF}')
echo "[01prepare] 完了: ${INPUT} (${n_samples} サンプル x ${n_vars} 遺伝子)"
if [[ -f "${SAMPLES}" ]]; then
  echo "[01prepare]       群 (組織) 数: $(tail -n +2 "${SAMPLES}" | cut -f3 | sort -u | wc -l)"
else
  echo "[01prepare]       samples.tsv なし -> 05importance_groups.sh はスキップされます"
fi

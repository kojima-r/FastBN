#!/usr/bin/env bash
# =============================================================================
# 00make_data.sh — ステップ 0: ダミーデータの生成
#   バルク RNA-seq 風の生カウント行列を、既知の DAG から生成する。
#   実データで解析する場合はこのステップを飛ばし、config.sh の EXPR_INPUT /
#   SAMPLE_META を自分のファイルに向けるだけでよい。
#
# 出力:
#   data/counts.tsv       : 行=遺伝子, 列=サンプル の生カウント (gene_id/gene_name/gene_length 付き)
#   data/sample_meta.tsv  : sample_id / group / replicate / library_size
#   data/true_edges.tsv   : 真の DAG (08evaluate.sh で学習結果と比較する)
#   target_genes.txt      : 注目遺伝子リスト (ハブ + 群応答遺伝子)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

python3 "${BN_SCRIPTS}/make_dummy_expr.py" \
  --outdir "${DATADIR}" \
  --out-targets "${TARGET_FILE}" \
  --n-genes "${DUMMY_GENES}" \
  --n-noise "${DUMMY_NOISE}" \
  --groups "${DUMMY_GROUPS}" \
  --n-replicates "${DUMMY_REPLICATES}" \
  --max-parents "${DUMMY_MAX_PARENTS}" \
  --edge-prob "${DUMMY_EDGE_PROB}" \
  --signal-frac "${DUMMY_SIGNAL_FRAC}" \
  --seed "${DUMMY_SEED}" \
  "$@"

echo "[00make_data] 完了: ${DATADIR}/counts.tsv, ${DATADIR}/sample_meta.tsv,"
echo "              ${DATADIR}/true_edges.tsv, ${TARGET_FILE}"

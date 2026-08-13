#!/usr/bin/env bash
# =============================================================================
# 08evaluate.sh — ステップ 8: 正解構造との比較 (ダミーデータ限定)
#   ダミーデータは真の DAG が既知 (data/true_edges.tsv) なので、学習結果の
#   精度を測れる。実データにこのステップは無いが、既知パスウェイを正解として
#   与えれば同じように使える。
#
# 出力:
#   out/eval_hc.tsv         : 学習網の評価指標
#   out/eval_hc_edges.tsv   : エッジ単位の判定 (TP / FP / FP_reversed / FN)
#   out/eval_bs.tsv         : コンセンサス網の評価指標 (04 実行後)
#   out/eval_bs_edges.tsv
#
# --restrict-to-analyzed: 正解エッジのうち「両端の遺伝子が前処理後も残っている」
#   ものだけを分母にする。分散フィルタで落ちた遺伝子のエッジは原理的に
#   復元できないため、こちらが学習器自体の精度に近い。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

true_edges="${DATADIR}/true_edges.tsv"
if [[ ! -f "${true_edges}" ]]; then
  echo "[08evaluate] ${true_edges} がありません (ダミーデータ以外ではスキップ)" >&2
  exit 0
fi

echo "=================================================================="
echo " (a) Hill-Climb 学習網 vs 真の DAG"
echo "=================================================================="
python3 "${BN_SCRIPTS}/compare_edges.py" \
  --true "${true_edges}" \
  --edges "${OUTDIR}/edges.tsv" --edges-named "${OUTDIR}/edges_named.tsv" \
  --input "${INPUT}" --restrict-to-analyzed \
  --out "${OUTDIR}/eval_hc.tsv" --out-edges "${OUTDIR}/eval_hc_edges.tsv"

if [[ -f "${OUTDIR}/integ_edges2.tsv" ]]; then
  echo "=================================================================="
  echo " (b) ブートストラップ・コンセンサス網 vs 真の DAG"
  echo "=================================================================="
  python3 "${BN_SCRIPTS}/compare_edges.py" \
    --true "${true_edges}" \
    --edges "${OUTDIR}/integ_edges2.tsv" \
    --edges-named "${OUTDIR}/integ_edges_named.tsv" \
    --input "${INPUT}" --restrict-to-analyzed \
    --out "${OUTDIR}/eval_bs.tsv" --out-edges "${OUTDIR}/eval_bs_edges.tsv"
fi

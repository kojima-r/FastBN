#!/usr/bin/env bash
# =============================================================================
# importance_groups.sh
#   学習済みネットワーク (構造 + カウント表) を**固定したまま**、サンプルを
#   実験群 (サブセット) に分けて `--score-dataset` に各群を与え、群ごとの
#   エッジ重要度を計算する (汎用版)。
#   「どのエッジがどの条件で効いているか」を見るための解析。
#
# 群の定義:
#   make_groups.py に渡す情報を環境変数で指定する (いずれか)。
#     SAMPLES      : preprocess_expr.py --out-samples の出力 (既定 ./data/samples.tsv)
#     GROUP_META   : サンプル ID + 群ラベルの表 (行順が INPUT と一致している前提)
#     GROUP_LABELS + GROUP_SIZES : 「先頭から n 件ずつ」を手動指定
#
# 対象ネットワークの切り替え (環境変数):
#   ・学習網 (既定):      INIT=${OUTDIR}/edges.tsv       COUNTS=${OUTDIR}/all_counts.tsv
#                        REF_EDGES=${OUTDIR}/edges.tsv   REF_NAMED=${OUTDIR}/edges_named.tsv
#   ・コンセンサス網:      INIT=${OUTDIR}/integ_edges.tsv COUNTS=${OUTDIR}/integ_all_counts.tsv
#                        REF_EDGES=${OUTDIR}/integ_edges2.tsv
#                        REF_NAMED=${OUTDIR}/integ_edges_named.tsv
#                        OUT_PREFIX=integ_edge_importance
#
# 出力:
#   ${GROUPDIR}/expr_g{N}_<label>.tsv          : 群ごとの score-dataset
#   ${OUTDIR}/<OUT_PREFIX>_g{N}_<label>.tsv    : 群ごとのエッジ重要度
#
# 注意: ノード番号 = 列位置 なので、INPUT の列順は対象網の学習入力と一致して
#       いなければならない。実行前に check_column_alignment.py で検証する。
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="groups"; bn_tag="groups"

init="${INIT:-${OUTDIR}/edges.tsv}"
counts="${COUNTS:-${OUTDIR}/all_counts.tsv}"
ref_edges="${REF_EDGES:-${OUTDIR}/edges.tsv}"
ref_named="${REF_NAMED:-${OUTDIR}/edges_named.tsv}"
out_prefix="${OUT_PREFIX:-edge_importance}"
group_prefix="${GROUP_PREFIX:-expr}"
score="${SCORE_IMP:-${SCORE:-bic}}"
alpha="${ALPHA:-1.0}"
ess="${ESS:-10.0}"
min_group_samples="${MIN_GROUP_SAMPLES:-2}"

require_bin
require_file "${INPUT}" "${init}" "${counts}" "${ref_edges}" "${ref_named}"
mkdir -p "${OUTDIR}" "${GROUPDIR}"

log "対象網: init=${init}, counts=${counts}, 出力prefix=${out_prefix}"

# --- 列順アライメント検証 ---------------------------------------------------
log "アライメント検証中 (${INPUT} <-> ${ref_edges}) ..."
"${PYTHON_BIN}" "$(py_tool check_column_alignment.py)" \
  "${INPUT}" "${ref_edges}" "${ref_named}"

# --- 群分割 -----------------------------------------------------------------
group_opt=()
if [[ -n "${GROUP_LABELS:-}" && -n "${GROUP_SIZES:-}" ]]; then
  group_opt=(--labels "${GROUP_LABELS}" --sizes "${GROUP_SIZES}")
  log "群の定義: GROUP_LABELS / GROUP_SIZES"
elif [[ -n "${GROUP_META:-}" ]]; then
  require_file "${GROUP_META}"
  group_opt=(--meta "${GROUP_META}")
  [[ -n "${GROUP_META_SAMPLE_COL:-}" ]] && group_opt+=(--sample-col "${GROUP_META_SAMPLE_COL}")
  [[ -n "${GROUP_META_GROUP_COL:-}" ]] && group_opt+=(--group-col "${GROUP_META_GROUP_COL}")
  log "群の定義: ${GROUP_META}"
elif [[ -f "${SAMPLES}" ]]; then
  group_opt=(--samples "${SAMPLES}")
  log "群の定義: ${SAMPLES}"
else
  die "群の定義が見つかりません。preprocess_expr.py に --out-samples ${SAMPLES} を
  指定して実行するか、GROUP_META / (GROUP_LABELS と GROUP_SIZES) を設定してください。"
fi

"${PYTHON_BIN}" "$(py_tool make_groups.py)" \
  --input "${INPUT}" --outdir "${GROUPDIR}" --prefix "${group_prefix}" \
  --min-samples "${min_group_samples}" \
  "${group_opt[@]}"

manifest="${GROUPDIR}/groups_manifest.tsv"
require_file "${manifest}"

# --- 群ごとのエッジ重要度 ---------------------------------------------------
while IFS=$'\t' read -r gno label nsamp gfile; do
  outimp="${OUTDIR}/${out_prefix}_g${gno}_${label}.tsv"
  log "群 ${gno} '${label}' (${nsamp} サンプル) -> ${outimp}"
  "${FASTBN_BIN}" --score "${score}" \
    --edge-importance \
    --score-dataset "${gfile}" \
    --init   "${init}" \
    --counts "${counts}" \
    --alpha "${alpha}" \
    --ess "${ess}" \
    --save-edge-importance "${outimp}"
done < <(tail -n +2 "${manifest}")

log "完了: ${OUTDIR}/${out_prefix}_g*_*.tsv"
log "注意: 群あたりのサンプル数が少ないと ΔlogL のノイズが大きく、絶対値は"
log "      全サンプル版より小さくなります (logL はサンプル数に比例)。"
log "      群間の相対比較として解釈してください。"

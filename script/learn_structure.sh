#!/usr/bin/env bash
# =============================================================================
# learn_structure.sh
#   離散化済み発現データからベイジアンネットワーク構造を学習する (汎用版)。
#   fast_bn の Hill-Climb + Tabu サーチを呼ぶだけの薄いラッパで、設定は
#   すべて環境変数で上書きできる。
#
# 前提: カレントディレクトリが解析ディレクトリであること。
#       INPUT (既定 ./data/expr_disc.tsv) が preprocess_expr.py で作成済み。
#
# 出力 (OUTDIR, 既定 ./out):
#   edges.tsv        : 学習された DAG のエッジ (ノード=列インデックス: u v)
#   edges_named.tsv  : 同じエッジを遺伝子名で表記 (edges.tsv と行対応)
#   all_counts.tsv   : CPT 推定 / エッジ重要度評価に用いるカウント表
#   log_learn.txt    : fast_bn のログ (--verbose 時に詳しい)
#
# 使い方:
#   ../script/learn_structure.sh
#   SCORE=bic MAX_PARENTS=2 ITERS=20000 ../script/learn_structure.sh
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="learn"; bn_tag="learn"

require_bin
require_file "${INPUT}"
mkdir -p "${OUTDIR}"

# --- 学習パラメータ (環境変数で上書き可能) ----------------------------------
score="${SCORE:-bdeu}"          # bic | k2 | bdeu
ess="${ESS:-10}"                # BDeu の等価サンプルサイズ
alpha="${ALPHA:-1.0}"           # スムージング係数
tabu="${TABU:-30}"              # Tabu サーチの禁制期間
iters="${ITERS:-10000}"         # Hill-Climb の最大反復数
topk="${TOPK:-20}"              # 候補親の上位 K 制限
maxpar="${MAX_PARENTS:-3}"      # 各ノードの最大親数
maxchild="${MAX_CHILDREN:-0}"   # 各ノードの最大子数 (0 = 無制限)
candmetric="${CAND_METRIC:-mi}" # 候補親の関連度指標 (mi | chi2)
jcache="${JINDEX_CACHE:-1024}"  # 親配置キャッシュ
reach="${REACH:-lazy}"          # 到達可能性チェック (dense | lazy)
verbose="${VERBOSE:-0}"         # 1 で --verbose

opt=()
[[ "${maxchild}" -gt 0 ]] && opt+=(--max-children "${maxchild}")
[[ "${verbose}" != "0" ]] && opt+=(--verbose)
# 初期構造 (任意): INIT を指定すると温かくスタートする
[[ -n "${INIT:-}" ]] && { require_file "${INIT}"; opt+=(--init "${INIT}"); }

n_samples=$(( $(wc -l < "${INPUT}") - 1 ))
n_vars=$(head -1 "${INPUT}" | awk -F'\t' '{print NF}')
log "入力: ${INPUT} (${n_samples} サンプル x ${n_vars} 変数)"
log "score=${score} ess=${ess} max-parents=${maxpar} iters=${iters} topk=${topk}"

"${FASTBN_BIN}" \
  --input "${INPUT}" \
  --score "${score}" \
  --ess "${ess}" --alpha "${alpha}" \
  --tabu "${tabu}" --iters "${iters}" \
  --topk "${topk}" --cand-metric "${candmetric}" \
  --max-parents "${maxpar}" \
  --reach "${reach}" --jindex-cache "${jcache}" \
  "${opt[@]}" \
  --save        "${OUTDIR}/edges.tsv" \
  --save-names  "${OUTDIR}/edges_named.tsv" \
  --save-counts "${OUTDIR}/all_counts.tsv" \
  2>&1 | tee "${OUTDIR}/log_learn.txt"

log "完了: $(wc -l < "${OUTDIR}/edges.tsv") エッジ -> ${OUTDIR}/edges.tsv, edges_named.tsv, all_counts.tsv"

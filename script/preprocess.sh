#!/usr/bin/env bash
# =============================================================================
# preprocess.sh
#   preprocess_expr.py を環境変数から呼び出すラッパ (汎用版)。
#   バルク RNA 発現量データ (TSV/CSV/XLSX) -> 正規化 -> log 変換 -> フィルタ
#   -> 離散化 -> fast_bn 入力 TSV。
#
# 必須:
#   EXPR_INPUT      : 発現量ファイル (TSV/CSV/XLSX)
#
# 主な設定 (環境変数; 詳細は preprocess_expr.py --help)
#   ORIENTATION     : genes-in-rows (既定) | samples-in-rows
#   ID_COL          : 遺伝子 ID の列 (列名 or 0 始まり位置; 既定=先頭列)
#   NAME_COL        : 遺伝子シンボルの列
#   LENGTH_COL      : 遺伝子長の列 (NORMALIZE=tpm のとき必要)
#   DROP_COLS       : 無視する注釈列 (カンマ区切り)
#   SHEET           : Excel のシート名
#   HEADER_ROW      : ヘッダ行の位置 (0 始まり; 既定 0)
#   SAMPLE_META     : サンプル ID と群ラベルの表 (群別解析に必要)
#   NORMALIZE       : none (既定) | cpm | tpm
#   N_BINS          : 離散化の段階数 (既定 3)
#   DISC_METHOD     : quantile (既定) | uniform
#   TOP_VAR_GENES   : 分散上位 N 遺伝子 (既定 500; 0 で無効)
#   VAR_QUANTILE    : 分散の分位点で選択 (TOP_VAR_GENES より優先)
#   MIN_MEAN_LOG    : 低発現フィルタ (既定 0 = 無効)
#   MIN_DETECT_FRAC : 検出率フィルタ (既定 0 = 無効)
#   TARGET_FILE     : 必ず残す注目遺伝子リスト (既定 ./target_genes.txt; 無ければ無視)
#   PREPROCESS_OPTS : preprocess_expr.py への追加オプション (そのまま渡す)
#
# 出力:
#   INPUT    (既定 ./data/expr_disc.tsv) : fast_bn 入力
#   VARMAP   (既定 ./data/var_map.tsv)   : 列インデックス <-> 遺伝子 対応表
#   SAMPLES  (既定 ./data/samples.tsv)   : 行番号 <-> サンプル ID / 群ラベル
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BN_TAG="preprocess"; bn_tag="preprocess"

[[ -n "${EXPR_INPUT:-}" ]] || die "EXPR_INPUT (発現量ファイル) を指定してください。"
require_file "${EXPR_INPUT}"
mkdir -p "$(dirname "${INPUT}")"

opt=(--input "${EXPR_INPUT}"
     --orientation "${ORIENTATION:-genes-in-rows}"
     --normalize "${NORMALIZE:-none}"
     --n-bins "${N_BINS:-3}"
     --disc-method "${DISC_METHOD:-quantile}"
     --min-mean-log "${MIN_MEAN_LOG:-0.0}"
     --min-detect-frac "${MIN_DETECT_FRAC:-0.0}"
     --out "${INPUT}" --out-map "${VARMAP}" --out-samples "${SAMPLES}")

[[ -n "${FORMAT:-}"        ]] && opt+=(--format "${FORMAT}")
[[ -n "${SHEET:-}"         ]] && opt+=(--sheet "${SHEET}")
[[ -n "${HEADER_ROW:-}"    ]] && opt+=(--header-row "${HEADER_ROW}")
[[ -n "${ID_COL:-}"        ]] && opt+=(--id-col "${ID_COL}")
[[ -n "${NAME_COL:-}"      ]] && opt+=(--name-col "${NAME_COL}")
[[ -n "${LENGTH_COL:-}"    ]] && opt+=(--length-col "${LENGTH_COL}")
[[ -n "${DROP_COLS:-}"     ]] && opt+=(--drop-cols "${DROP_COLS}")
[[ -n "${SAMPLE_ID_COL:-}" ]] && opt+=(--sample-id-col "${SAMPLE_ID_COL}")
[[ -n "${PSEUDOCOUNT:-}"   ]] && opt+=(--pseudocount "${PSEUDOCOUNT}")
[[ -n "${DETECT_THRESHOLD:-}" ]] && opt+=(--detect-threshold "${DETECT_THRESHOLD}")
[[ -n "${GROUP_ORDER:-}"   ]] && opt+=(--group-order "${GROUP_ORDER}")
[[ "${LOG2:-1}" == "0"     ]] && opt+=(--no-log2)

# 分散フィルタ: VAR_QUANTILE があれば優先
if [[ -n "${VAR_QUANTILE:-}" ]]; then
  opt+=(--var-quantile "${VAR_QUANTILE}")
else
  opt+=(--top-var-genes "${TOP_VAR_GENES:-500}")
fi

# サンプルメタデータ (群ラベル)
if [[ -n "${SAMPLE_META:-}" ]]; then
  require_file "${SAMPLE_META}"
  opt+=(--sample-meta "${SAMPLE_META}")
  [[ -n "${META_SAMPLE_COL:-}" ]] && opt+=(--meta-sample-col "${META_SAMPLE_COL}")
  [[ -n "${META_GROUP_COL:-}"  ]] && opt+=(--meta-group-col "${META_GROUP_COL}")
fi

# 注目遺伝子ホワイトリスト (あれば使う)
if [[ -n "${TARGET_FILE:-}" && -f "${TARGET_FILE}" ]]; then
  opt+=(--keep-genes-file "${TARGET_FILE}")
  log "ホワイトリスト ${TARGET_FILE} を使用"
fi
[[ -n "${KEEP_GENES:-}" ]] && opt+=(--keep-genes "${KEEP_GENES}")

# 追加オプション (単語分割して渡す)
extra=()
if [[ -n "${PREPROCESS_OPTS:-}" ]]; then
  # shellcheck disable=SC2206
  extra=(${PREPROCESS_OPTS})
fi

"${PYTHON_BIN}" "$(py_tool preprocess_expr.py)" "${opt[@]}" "${extra[@]}"

log "完了: ${INPUT} (対応表: ${VARMAP}, サンプル表: ${SAMPLES})"

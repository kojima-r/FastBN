#!/usr/bin/env bash
# =============================================================================
# new_analysis.sh
#   解析ディレクトリを作り、inspect_matrix.py の推奨値を入れた config.sh を置く。
#   ユーザのデータは移動もコピーもしない (config.sh から絶対パスで参照する)。
#
#   使い方:
#     bash new_analysis.sh <解析ディレクトリ> --expr <発現量ファイル> \
#          [--meta <サンプル情報>] [--sheet <Excel シート>] [--targets <遺伝子リスト>]
#     bash new_analysis.sh <解析ディレクトリ>          # 空の雛形 (example_bulk/config.sh をコピー)
#
#   inspect_matrix.py に渡したい追加オプション (--orientation, --header-row など) は
#   -- のあとに書くとそのまま素通しする。
# =============================================================================
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
die() { echo "[new_analysis] エラー: $*" >&2; exit 1; }
log() { echo "[new_analysis] $*" >&2; }

[[ $# -ge 1 ]] || die "解析ディレクトリを指定してください (--help は new_analysis.sh の先頭を参照)"
target_dir="$1"; shift

expr_file=""; meta_file=""; targets_file=""; sheet=""; passthru=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --expr)    expr_file="$2"; shift 2 ;;
    --meta)    meta_file="$2"; shift 2 ;;
    --targets) targets_file="$2"; shift 2 ;;
    --sheet)   sheet="$2"; shift 2 ;;
    --)        shift; passthru+=("$@"); break ;;
    *)         passthru+=("$1"); shift ;;
  esac
done

# FastBN の場所を解決 (FASTBN_HOME / BN_SCRIPTS / PYTHON_BIN が入る)
eval "$(bash "${script_dir}/fastbn_env.sh")"

mkdir -p "${target_dir}"
target_dir="$(cd "${target_dir}" && pwd)"
config="${target_dir}/config.sh"

if [[ -f "${config}" ]]; then
  die "${config} が既にあります。上書きしたくないので中断します
  (新しい設定を試すなら別ディレクトリを作る、または RUNDIR を変えて実行する)"
fi

if [[ -n "${expr_file}" ]]; then
  [[ -f "${expr_file}" ]] || die "発現量ファイルがありません: ${expr_file}"
  opt=(--input "${expr_file}" --emit-config "${config}")
  [[ -n "${meta_file}" ]] && { [[ -f "${meta_file}" ]] || die "メタデータがありません: ${meta_file}"; opt+=(--meta "${meta_file}"); }
  [[ -n "${sheet}" ]] && opt+=(--sheet "${sheet}")
  "${PYTHON_BIN}" "${script_dir}/inspect_matrix.py" "${opt[@]}" "${passthru[@]+"${passthru[@]}"}"
else
  log "発現量ファイルが未指定なので example_bulk/config.sh を雛形としてコピーします"
  cp "${FASTBN_HOME}/example_bulk/config.sh" "${config}"
  log "${config} の EXPR_INPUT / SAMPLE_META / ID_COL / NAME_COL / NORMALIZE を書き換えてください"
  log "(DUMMY_* はダミーデータ生成用の項目なので無視してよい)"
fi

# 注目遺伝子リスト (あれば置く。図で赤く強調され、フィルタから免除される)
if [[ -n "${targets_file}" ]]; then
  [[ -f "${targets_file}" ]] || die "注目遺伝子リストがありません: ${targets_file}"
  cp "${targets_file}" "${target_dir}/target_genes.txt"
  log "注目遺伝子: $(wc -l < "${target_dir}/target_genes.txt") 件 -> ${target_dir}/target_genes.txt"
fi

cat >&2 <<EOF
[new_analysis] 完了: ${target_dir}
  次の手順:
    1. ${config} を読んで、列指定と正規化が実データと合っているか確認する
    2. cd ${target_dir} && source ./config.sh
    3. "\${BN_SCRIPTS}/preprocess.sh" から段階的に実行する
       (初回は TOP_VAR_GENES を小さく、BOOTSTRAP=3 SEEDS=2 で所要時間を測る)
EOF

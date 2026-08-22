#!/usr/bin/env bash
# =============================================================================
# fastbn_env.sh
#   FastBN の場所を解決し、バイナリと Python 依存を用意して、環境変数を
#   標準出力に export 文として書き出す。人向けのメッセージは標準エラーへ出す。
#
#   使い方:
#     eval "$(bash fastbn_env.sh)"     # 解決 + 必要ならビルド + export
#     bash fastbn_env.sh --check       # 状態の報告のみ (ビルドしない)
#     bash fastbn_env.sh --build       # 強制的に再ビルド
#
#   FASTBN_HOME の解決順:
#     1) 環境変数 FASTBN_HOME
#     2) このスクリプトから見たプラグイン/リポジトリのルート (../../..)
#     3) CLAUDE_PLUGIN_ROOT
#     4) カレントディレクトリから上に遡って探索
#     5) CLAUDE_PROJECT_DIR / よくある置き場 (~/FastBN, /opt/FastBN)
#   「FastBN のルート」の判定条件は fast_bn.cpp と script/common.sh があること。
# =============================================================================
set -euo pipefail

mode="ensure"
case "${1:-}" in
  --check) mode="check" ;;
  --build) mode="build" ;;
  -h|--help)
    sed -n '2,20p' "${BASH_SOURCE[0]}" >&2
    exit 0 ;;
  "") ;;
  *) echo "fastbn_env.sh: 不明な引数: $1 (--check | --build)" >&2; exit 2 ;;
esac

msg()  { echo "[fastbn_env] $*" >&2; }
die()  { echo "[fastbn_env] エラー: $*" >&2; exit 1; }

is_fastbn_home() {
  [[ -n "${1:-}" && -f "$1/fast_bn.cpp" && -f "$1/script/common.sh" ]]
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- 1) FASTBN_HOME の解決 ---------------------------------------------------
home=""
for cand in \
    "${FASTBN_HOME:-}" \
    "$(cd "${script_dir}/../../.." 2>/dev/null && pwd || true)" \
    "${CLAUDE_PLUGIN_ROOT:-}" ; do
  if is_fastbn_home "${cand}"; then home="${cand}"; break; fi
done

# カレントディレクトリから上へ (解析ディレクトリの中で呼ばれた場合)
if [[ -z "${home}" ]]; then
  d="${PWD}"
  while [[ "${d}" != "/" ]]; do
    if is_fastbn_home "${d}"; then home="${d}"; break; fi
    d="$(dirname "${d}")"
  done
fi

if [[ -z "${home}" ]]; then
  for cand in "${CLAUDE_PROJECT_DIR:-}" "${HOME}/FastBN" "/opt/FastBN"; do
    if is_fastbn_home "${cand}"; then home="${cand}"; break; fi
  done
fi

if [[ -z "${home}" ]]; then
  die "FastBN のリポジトリが見つかりません。
  次のいずれかで場所を教えてください:
    export FASTBN_HOME=/path/to/FastBN
  まだ持っていない場合:
    git clone https://github.com/kojima-r/FastBN.git && export FASTBN_HOME=\$PWD/FastBN"
fi

home="$(cd "${home}" && pwd)"
scripts="${home}/script"

# --- 2) fast_bn バイナリ -----------------------------------------------------
bin="${FASTBN_BIN:-${home}/fast_bn}"
need_build=0
reason=""
if [[ ! -x "${bin}" ]]; then
  need_build=1; reason="バイナリがありません: ${bin}"
elif [[ "${home}/fast_bn.cpp" -nt "${bin}" ]]; then
  need_build=1; reason="fast_bn.cpp がバイナリより新しい (再ビルドが必要)"
fi
[[ "${mode}" == "build" ]] && { need_build=1; reason="--build が指定された"; }

if [[ "${need_build}" -eq 1 && "${mode}" == "check" ]]; then
  msg "要ビルド: ${reason}"
elif [[ "${need_build}" -eq 1 ]]; then
  msg "ビルドします (${reason})"
  command -v g++ >/dev/null 2>&1 || die "g++ がありません (g++ 13 以上 / C++17 が必要)"
  ( cd "${home}" && g++ -O3 -march=native -std=c++17 fast_bn.cpp -o fast_bn ) >&2 \
    || die "ビルドに失敗しました。${home} で ./compile.sh を手で実行して出力を確認してください。"
  bin="${home}/fast_bn"
  msg "ビルド完了: ${bin}"
fi

# --- 3) Python 依存 ----------------------------------------------------------
py="${PYTHON_BIN:-python3}"
command -v "${py}" >/dev/null 2>&1 || die "${py} がありません"
missing_core=""
missing_viz=""
for m in numpy pandas; do
  "${py}" -c "import ${m}" 2>/dev/null || missing_core="${missing_core} ${m}"
done
for m in networkx matplotlib; do
  "${py}" -c "import ${m}" 2>/dev/null || missing_viz="${missing_viz} ${m}"
done
[[ -n "${missing_core}" ]] && msg "警告: 前処理に必要な Python パッケージが不足:${missing_core} (pip install${missing_core})"
[[ -n "${missing_viz}"  ]] && msg "警告: 可視化に必要な Python パッケージが不足:${missing_viz} (pip install${missing_viz})"

# --- 4) 報告と export --------------------------------------------------------
if [[ "${mode}" == "check" ]]; then
  msg "FASTBN_HOME = ${home}"
  msg "FASTBN_BIN  = ${bin} $( [[ -x "${bin}" ]] && echo '(実行可能)' || echo '(未ビルド)' )"
  msg "BN_SCRIPTS  = ${scripts}"
  msg "PYTHON_BIN  = ${py} ($("${py}" --version 2>&1))"
  [[ -x "${bin}" ]] && msg "動作確認: $(LANG=en "${bin}" --help 2>&1 | head -1)"
fi

cat <<EOF
export FASTBN_HOME="${home}"
export FASTBN_BIN="${bin}"
export BN_SCRIPTS="${scripts}"
export BS_PROB_PY="${home}/compute_bs_prob.py"
export PYTHON_BIN="${py}"
EOF

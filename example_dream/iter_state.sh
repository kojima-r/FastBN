#!/usr/bin/env bash
# =============================================================================
# iter_state.sh — 学習ログから「反復予算が律速しているか」を判定する
#
#   反復の最後の 10% でまだ最良スコア ([*] new best) が更新されていれば、
#   反復を増やせばまだ改善する = 反復予算が律速している と判定する。
#   Tabu サーチは局所最適でも動き続けるため [stop] はほとんど出ず、
#   「収束したか」の判定にはこちらの方が実用的。
#
#   iter_state.sh <log>            -> 人が読む文字列
#   iter_state.sh --flag <log>     -> 1 (律速) / 0 (十分)
# =============================================================================
set -uo pipefail
flag=0
if [[ "${1:-}" == "--flag" ]]; then flag=1; shift; fi
log="${1:-}"
if [[ ! -s "${log}" ]]; then
  [[ "${flag}" -eq 1 ]] && echo 0 || echo "ログなし"
  exit 0
fi
total=$(grep -c "^\[it " "${log}" 2>/dev/null); total="${total:-0}"
if [[ "${total}" -eq 0 ]]; then
  [[ "${flag}" -eq 1 ]] && echo 0 || echo "反復なし"
  exit 0
fi
# 最後に最良スコアが更新された反復番号
last_best=$(grep -n "^\[\*\] new best" "${log}" 2>/dev/null | tail -1 | cut -d: -f1)
last_best="${last_best:-0}"
# 行番号ではなく反復番号に直す (直前の [it N] を探す)
it_at_best=$(head -n "${last_best}" "${log}" 2>/dev/null | grep "^\[it " | tail -1 \
             | sed -n 's/^\[it \([0-9]*\)\].*/\1/p')
it_at_best="${it_at_best:-0}"
threshold=$(( total * 9 / 10 ))
if [[ "${it_at_best}" -ge "${threshold}" ]]; then
  [[ "${flag}" -eq 1 ]] && echo 1 || echo "**反復予算が律速** (最良更新 it=${it_at_best}/${total})"
else
  [[ "${flag}" -eq 1 ]] && echo 0 || echo "反復十分 (最良更新 it=${it_at_best}/${total})"
fi

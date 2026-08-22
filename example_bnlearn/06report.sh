#!/usr/bin/env bash
# =============================================================================
# 06report.sh — ステップ 6: HTML レポート
#   集約表・指標グラフ・正解 vs 学習ネットワークの比較図を 1 つの HTML に
#   まとめる。
#
# 出力: ${REPORT_HTML}
#
#   ./06report.sh            # 画像は相対リンク (軽量)
#   ./06report.sh --embed    # 画像を base64 埋め込み (1 ファイルで共有可)
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

[[ -s "${BENCHMARK}" ]] || { echo "[06report] エラー: ${BENCHMARK} がありません。先に ./03evaluate.sh" >&2; exit 1; }

# 対象ネットワークの一覧表 (BIF から直接読む)
netinfo="${RUNDIR}/networks.tsv"
mkdir -p "${RUNDIR}"
{
  printf 'network\tnodes\tarcs\tparameters\tmax_in_degree\tstate_space\n'
  for net in ${NETWORKS}; do
    python3 "${BN_SCRIPTS}/bif_io.py" info --bif "${NETDIR}/${net}.bif" \
      | awk -F'\t' -v n="${net}" '
          $1=="nodes"{a=$2} $1=="arcs"{b=$2} $1=="parameters"{c=$2}
          $1=="max_in_degree"{d=$2} $1=="state_space"{e=$2}
          END{printf "%s\t%s\t%s\t%s\t%s\t%s\n", n, a, b, c, d, e}'
  done
} > "${netinfo}"

python3 "${BN_SCRIPTS}/make_benchmark_report.py" \
  --benchmark "${BENCHMARK}" \
  --summary "${SUMMARY_TSV}" \
  --summary-overall "${SUMMARY_OVERALL_TSV}" \
  --networks "${netinfo}" \
  --plot "${SUMMARY_PNG}" \
  --compare-dir "${FIGDIR}" \
  --out "${REPORT_HTML}" \
  --title "bnlearn discrete-small 構造学習ベンチマーク" \
  --subtitle "ネットワーク: ${NETWORKS} / サンプル数: ${SAMPLE_SIZES} (各 ${REPLICATES} 反復) / スコア: ${SCORES}" \
  "$@"

echo "[06report] ブラウザで開いて確認してください: $(pwd)/${REPORT_HTML#./}"

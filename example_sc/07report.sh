#!/usr/bin/env bash
# =============================================================================
# 07report.sh — ステップ 7: HTML レポートの作成
#   図・エッジ重要度テーブル・データ要約を 1 つの HTML (${REPORT_HTML}) に
#   集約する。上部のボタンでメトリクス (dlogL/dBIC/...) を切り替えられる。
#
#   既定では図を相対リンクで参照する (軽量)。単体で共有・移動したい場合は
#   --embed で画像を base64 埋め込みする。
#     ./07report.sh --embed
#     ./07report.sh --metrics dlogL --top-edges 60
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh

"${BN_SCRIPTS}/make_report.sh" "$@"

echo "[07report] ブラウザで開いて確認してください: $(pwd)/${REPORT_HTML#./}"

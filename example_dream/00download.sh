#!/usr/bin/env bash
# =============================================================================
# 00download.sh — ステップ 0: DREAM データの取得
#
#   dream4 : GeneNetWeaver 配布の "DREAM4 in silico challenge.zip" (約 77 MB)
#            https://gnw.sourceforge.net/dreamchallenge.html
#   dream5 : Zenodo 17854236 の 1_Challenge_Data_Supplement.zip (約 38 MB)
#            https://zenodo.org/records/17854236  (CC-BY-ND-4.0)
#   hpn    : HPN-DREAM は Synapse (syn1720047) にあり、**認証が必要**なので
#            自動取得できません。下記の手順で手動配置してください。
#
# サーバが途中で接続を切ることがあるので、-C - で再開しながら取得します。
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
mkdir -p "${SRCDIR}"

fetch() {  # <出力パス> <URL> <説明>
  local out="$1" url="$2" what="$3"
  if [[ -s "${out}" ]] && python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1])" "${out}" 2>/dev/null; then
    echo "[00download] ${what}: 取得済み (スキップ)"
    return 0
  fi
  echo "[00download] ${what}: 取得中 ${url}"
  for attempt in 1 2 3 4 5; do
    curl -fsSL --retry 3 --retry-delay 2 -C - -o "${out}" "${url}" || true
    if python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1])" "${out}" 2>/dev/null; then
      echo "[00download] ${what}: 完了 ($(stat -c%s "${out}") bytes)"
      return 0
    fi
    echo "[00download] ${what}: 不完全 ($(stat -c%s "${out}" 2>/dev/null || echo 0) bytes) — 再開します (${attempt}/5)"
  done
  echo "[00download] エラー: ${what} を取得できませんでした" >&2
  exit 1
}

for ds in ${DATASETS}; do
  case "${ds}" in
    dream4) fetch "${SRCDIR}/dream4.zip" "${DREAM4_URL}" "DREAM4" ;;
    dream5) fetch "${SRCDIR}/dream5.zip" "${DREAM5_URL}" "DREAM5" ;;
    hpn)
      if [[ -s "${HPN_SRCDIR}/expression.tsv" && -s "${HPN_SRCDIR}/true_edges.tsv" ]]; then
        echo "[00download] HPN-DREAM: 手動配置されたデータを使います (${HPN_SRCDIR})"
      else
        cat >&2 <<'MSG'
[00download] HPN-DREAM: 自動取得できません。
    データは Synapse (https://www.synapse.org/#!Synapse:syn1720047) にあり、
    アカウント登録と利用規約への同意が必要です。匿名アクセスは 403 になります。

    手順:
      1. Synapse にログインし、HPN-DREAM Network Inference Challenge の
         Challenge Data からデータをダウンロードする
      2. 次の 2 ファイルを用意して配置する
           ${HPN_SRCDIR}/expression.tsv   行=サンプル, 列=変数, 1 行目=変数名
           ${HPN_SRCDIR}/true_edges.tsv   正解エッジ (u <TAB> v)
      3. DATASETS に hpn を含めて再実行する

    配置が無い場合、以降のステップは hpn を黙ってスキップします。
MSG
      fi
      ;;
    *) echo "[00download] 未知のデータセット: ${ds}" >&2; exit 1 ;;
  esac
done
echo "[00download] 完了"

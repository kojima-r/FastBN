---
description: FastBN の環境を確認する (バイナリのビルド・Python 依存・動作確認)
argument-hint: [--smoke-test]
---

引数: $ARGUMENTS

1. 環境を確認する:

   ```bash
   bash "${CLAUDE_PLUGIN_ROOT}/skills/expression-network-inference/scripts/fastbn_env.sh" --check
   ```

   `FASTBN_HOME` が解決できない場合は、FastBN の場所をユーザに尋ねる
   (未取得なら `git clone https://github.com/kojima-r/FastBN.git`)。

2. バイナリが無い / `fast_bn.cpp` より古い場合はビルドする:

   ```bash
   bash "${CLAUDE_PLUGIN_ROOT}/skills/expression-network-inference/scripts/fastbn_env.sh" --build
   ```

   g++ 13 以上 / C++17 が必要。Python は numpy・pandas (前処理)、
   networkx・matplotlib (可視化)、openpyxl (Excel 入力) を使う。

3. `--smoke-test` が指定された場合は、ダミーデータ生成から評価までの
   エンドツーエンド確認を実行する (1〜2 分):

   ```bash
   cd "${FASTBN_HOME}/example_bulk" && ./run_all.sh
   ```

   真の DAG が既知なので最後のステップが精度まで出す。ここが通れば環境は問題ない。

結果は簡潔に報告する (解決したパス、ビルドの有無、不足している依存、スモーク
テストの結果)。

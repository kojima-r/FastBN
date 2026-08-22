# CLI の癖・エラー対処

## fast_bn の CLI で事故になる点

* **未知のオプションは即エラー** (`Unknown option: ...`, exit 1)。寛容なフォールバックは
  一切ない。古いフラグ名を書くとその瞬間に死ぬ。
* `--bootstrap-include-zero` は `--help` に載っているが**パーサに分岐が無い**ので
  渡すと上のエラーになる。
* **ブートストラップのシードは `--seed`** (`--bootstrap-seed` は存在しない。
  `README.md` と `sample00_bs/` は古い)。
* ヘルプの言語は `LANG` で決まる (`en` で始まれば英語、既定は日本語)。
  `README.md` の `BN_LANG` は効かない。
* 詳細ログは既定で ON。`--verbose` と `--quiet` は最後に書いた方が勝つ。
* `--save-bootstrap-counts out/edges.tsv` は `out/edges_seed0004.tsv` を書く
  (シード mod 10000 を 4 桁ゼロ詰めで拡張子の前に挿入)。
* `--ess` は `bdeu` にしか効かない (他のスコアでは警告)。
* `--iters 0` は「探索せずに `--init` の構造のカウント表を作る」イディオム。
  コンセンサス網の `all_counts.tsv` はこれで作る (`bootstrap_stability.sh` の C)。
* `sample00/` と `sample00_bs/` は**古い**。削除された位置引数 CLI
  (`./fast_bn all_disc.tsv bdeu ...`) と `--bootstrap-seed` を使っているので
  使用例として読まない。

ドキュメント (`README.md`, `script/README.md`) は詳しいが一部古い。
食い違ったら `fast_bn.cpp` の引数パーサが正しい。

## よくあるエラーと対処

| 症状 | 原因 / 対処 |
| --- | --- |
| `fast_bn バイナリが見つかりません` | `${FASTBN_HOME}/compile.sh` を実行 (g++ 13+ / C++17)。別の場所なら `FASTBN_BIN` を指定 |
| `必須ファイルがありません: ./data/expr_disc.tsv` | 前処理をしていない、または CWD が解析ディレクトリでない |
| `EXPR_INPUT (発現量ファイル) を指定してください` | `source ./config.sh` を忘れている |
| 列順アライメント検証で中断 | `INPUT` の列順が対象網の学習入力と違う。前処理をやり直したなら、その学習に使った入力をそのまま渡す。古い `out/` を新しい行列に使い回していないか確認 |
| `コンセンサスエッジが 0 件でした` | `THRESHOLD_PROB` を下げる、または `BOOTSTRAP`/`SEEDS` を増やす |
| 辺が 0 本しか学習されない | `N_BINS` に対して N が足りない (`P_eff < 1`)、定数列だらけ、`TOPK` が小さすぎる |
| 辺が異常に多い (D の 3 倍以上) | `MAX_PARENTS` が大きすぎる、`SCORE=bdeu` の `ESS` が大きすぎる |
| メモリを食い尽くす | `q_i × r_i` の爆発。`MAX_PARENTS` と `N_BINS` を下げる。`REACH=dense` を使っているなら `lazy` に戻す |
| ブートストラップが終わらない | `BOOTSTRAP × SEEDS` 回の完全な学習。`MAX_JOBS` を CPU コア数に合わせ、`ITERS_BS` を `D × P_eff` 程度まで下げる |
| 図のラベルが読めない / 潰れる | `--top-n` を下げる、`--style-scale` で倍率を上げる、`--layout kamada` を試す |
| `report.html` に図が出ない | 図を作る前にレポートを作った。`viz*.sh` を先に実行 |

## 検算 (探索を疑う前に)

学習で報告されたスコアは、同じ構造をゼロから再計算したスコアと一致しなければ
ならない:

```bash
# 学習
"${BN_SCRIPTS}/learn_structure.sh"                      # out/log_learn.txt の最終スコア
# 同じ構造を探索なしで再スコア
"${FASTBN_BIN}" --input data/expr_disc.tsv --score bdeu --ess 10 \
  --iters 0 --init out/edges.tsv --save-counts /dev/null
```

一致しなければ探索側のバグ (過去に親桁順の不整合で実際に起きた)。
経緯と正確なコマンドは `${FASTBN_HOME}/example_bnlearn/README.md`。

同じ手で「真の構造のスコア」も計算できる (`--init <true edges>`)。学習網の
スコアが真の構造より高ければ、精度が出ないのは探索ではなくデータ量・離散化の
問題 (`interpretation.md` の最後)。

## 動作確認 (インストール直後 / 環境を疑うとき)

```bash
cd "${FASTBN_HOME}/example_bulk" && ./run_all.sh    # 1〜2 分。ダミーデータ生成から評価まで
```

真の DAG が既知なので step 08 が精度まで出す。ここが通れば環境は問題ない。

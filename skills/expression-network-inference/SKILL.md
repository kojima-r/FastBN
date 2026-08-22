---
name: expression-network-inference
description: This skill should be used when the user hands over an expression matrix (bulk RNA-seq counts/TPM, single-cell, proteomics, or any sample x feature table) and wants a network inferred, scored, stabilized or visualized — "estimate a gene network", "learn a Bayesian network / DAG from this data", "which edges are important", "bootstrap edge stability", "compare conditions", "draw the network", "発現データからネットワークを推定して", "遺伝子ネットワークを可視化", "エッジの重要度・安定性を調べたい", "群間で比較したい". It drives the FastBN pipeline — preprocess, Hill-Climb+Tabu structure learning, edge importance, bootstrap consensus, per-group comparison, figures, HTML report.
license: MIT
---

# 発現データからのベイジアンネットワーク推定・可視化 (FastBN)

`fast_bn` (C++ の Hill-Climb + Tabu 探索) と `script/` の汎用パイプラインを使い、
渡された発現量行列から DAG を推定し、エッジ重要度・ブートストラップ安定性・
群間比較・図・HTML レポートまで出す。

## 最初にやること (順番を守る)

### 0. 環境の解決 — 必ず最初

`SKILL_DIR` = **この `SKILL.md` があるディレクトリ** (ホストが提示したパスから決まる。
Claude Code のプラグインなら `${CLAUDE_PLUGIN_ROOT}/skills/expression-network-inference`、
Codex なら `<プラグイン>/skills/...` や `.codex/skills/...` / `${CODEX_HOME}/skills/...`、
FastBN のチェックアウト内なら `${FASTBN_HOME}/skills/...`)。

```bash
SKILL_DIR=<このファイルのあるディレクトリ>
eval "$(bash "${SKILL_DIR}/scripts/fastbn_env.sh")"
```

`FASTBN_HOME` / `FASTBN_BIN` / `BN_SCRIPTS` / `PYTHON_BIN` を export し、バイナリが
無ければビルドし、Python 依存 (numpy/pandas/networkx/matplotlib) を検査する。

FastBN 本体 (`fast_bn.cpp` と `script/`) が見つからないと言われた場合は、
`FASTBN_HOME` を教えてもらうか `git clone https://github.com/kojima-r/FastBN.git`
してから再実行する (詳細は `references/troubleshooting.md`)。


### 1. データを見てから設定を決める (推測で config を書かない)

```bash
python3 "${SKILL_DIR}/scripts/inspect_matrix.py" \
    --input <ユーザのファイル> [--meta <サンプル情報>] [--sheet <Excel シート>]
```

行列の向き・値の種類 (生カウント / 正規化済み / 既に離散化済み)・サンプル数 N・
変数数 D・欠損・重複を判定し、**そのデータに合った推奨パラメータ**
(`NORMALIZE` / `N_BINS` / `TOP_VAR_GENES` / `MAX_PARENTS` / `ITERS` / ブートストラップ回数)
を根拠つきで出す。推奨値の導出は `references/parameter-sizing.md`。

判定できない点 (どの列が遺伝子 ID か、生カウントか正規化済みか、群ラベルの列名) は
**ユーザに確認する**。ここを間違えると以降すべてが無意味になる。

### 2. 解析ディレクトリを作る

```bash
bash "${SKILL_DIR}/scripts/new_analysis.sh" \
    <解析ディレクトリ> --expr <ユーザのファイル> [--meta <サンプル情報>]
```

`<解析ディレクトリ>/config.sh` (推奨値入り) を生成する。ユーザのデータは**移動も
コピーもしない** — config.sh に絶対パスで参照させる。生成後に `config.sh` を読み、
`EXPR_INPUT` / `ID_COL` / `NAME_COL` / `DROP_COLS` / `NORMALIZE` / `SAMPLE_META` が
実データと合っているか自分で確認してから走らせる。

### 3. 実行

```bash
cd <解析ディレクトリ>
source ./config.sh
"${BN_SCRIPTS}/preprocess.sh"          # 1) 正規化 -> log -> フィルタ -> 離散化
"${BN_SCRIPTS}/learn_structure.sh"     # 2) 構造学習
"${BN_SCRIPTS}/edge_importance.sh"     # 3) エッジ重要度
"${BN_SCRIPTS}/bootstrap_stability.sh" # 4) 安定性 -> コンセンサス網 (最も重い)
"${BN_SCRIPTS}/importance_groups.sh"   # 5) 群別重要度 (群ラベルがある場合)
"${BN_SCRIPTS}/viz.sh" --metrics "${VIZ_METRICS}" --top-n "${VIZ_TOP_N}"   # 6) 図
"${BN_SCRIPTS}/make_report.sh"         # 7) report.html
```

`script/viz*.sh` は `VIZ_METRICS` / `VIZ_TOP_N` を**読まない** (環境変数を読むのは
`example_*/06visualize.sh` の側)。直接呼ぶときは上のように `--metrics` /
`--top-n` を渡す。渡さないと 4 メトリクス全部を描いて何倍も時間がかかる。

**段階ごとに実行する** (`run_pipeline.sh` の一括実行は、所要時間が読めている
再実行のときだけ)。理由: 2 で辺が 0 本や過剰に出たら 3 以降は無駄になり、4 は
`BOOTSTRAP × SEEDS` 回の学習で全体の実行時間を支配するため。

初回は必ず**小さく試す**。`TOP_VAR_GENES` を 100〜300、`BOOTSTRAP=3 SEEDS=2` で
1 周させて時間を測り、実測値からユーザに本番設定の所要時間を伝えて合意を取る。
数時間〜数日かかる規模 (D が数千以上) では特に必須。

各ステージの環境変数・出力・コストは `references/pipeline-stages.md`。

### 4. 結果を読んで報告する

図とテーブルを貼るだけで終わらせない。`references/interpretation.md` に従って、
エッジ重要度・ブートストラップ確率・群間差を解釈し、**マルコフ同値による向きの
不確かさ**を必ず添える。既知の正解構造やパスウェイと比べる依頼なら
`network-structure-evaluation` スキルへ。

## 絶対に守る 2 つの約束 (パイプラインの前提)

1. **カレントディレクトリ = 解析ディレクトリ。** `script/*.sh` は自分の場所へ
   `cd` しない。`./data`, `./out`, `./bs`, `./groups`, `./figures` はすべて CWD 基準。
   別ディレクトリから呼ぶと出力が散らばる。
2. **設定は環境変数だけ。** `script/` 以下のスクリプトを編集しない。解析ディレクトリの
   `config.sh` を `source` してから呼ぶ。1 回だけ変えたいなら
   `SCORE=bic ../script/learn_structure.sh` のように前置きする。
   `RUNDIR` を変えれば設定違いの実験を並べて比較できる。

## 落とし穴 (事故が起きる順)

* **ノード番号 = 入力 TSV の列位置。** `fast_bn` は遺伝子名を保存しない。
  前処理をやり直す (`TOP_VAR_GENES` を変える等) と番号が全部ずれるので、
  古い `out/edges.tsv` / `all_counts.tsv` を新しい行列に使うと結果は無意味になる。
  `edges.tsv` と `edges_named.tsv` は**行が 1:1 対応**しており、名前はここから復元する。
* **未知のオプションは即エラー** (`Unknown option: ...`, exit 1)。フォールバックは無い。
* **ブートストラップのシード指定は `--seed`** (`--bootstrap-seed` は存在しない)。
* `--iters 0` は no-op ではなく「`--init` で与えた構造のカウント表を再計算する」
  イディオム。コンセンサス網の `all_counts.tsv` はこれで作る。
* 群別解析の前に列順アライメント検証が走り、不一致なら中断する (これは正しい挙動)。

その他の既知の癖・エラー対処は `references/troubleshooting.md`。
リポジトリ内ドキュメント (`README.md` / `script/README.md`) は一部古い。
食い違ったら `fast_bn.cpp` の引数パーサが正しい。

## 参照ファイル

| ファイル | 中身 |
| --- | --- |
| `references/pipeline-stages.md` | 各ステージの環境変数・出力・計算コスト・実行順の分岐 |
| `references/preprocessing.md` | 入力形式の判定、正規化と離散化の選び方、離散化済みデータの扱い、注目遺伝子 |
| `references/parameter-sizing.md` | `ITERS` / `MAX_PARENTS` / `N_BINS` / ブートストラップ回数の決め方 (D, N, r からの導出) |
| `references/interpretation.md` | 重要度・安定性・群間比較・図の読み方と報告の書き方 |
| `references/troubleshooting.md` | CLI の癖、よくあるエラー、探索を疑う前に確認すること |

実データでの完全な実行例はリポジトリ内にある: `example_bulk/` (生カウントからの
バルク RNA、ダミーデータで真の構造つき)、`example_sc/` (離散化済み単一細胞)。
迷ったら該当する例の `config.sh` と番号付きスクリプトを読む。

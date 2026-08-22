# example_sc — 単一細胞発現データでのベイジアンネットワーク解析

**離散化済みの単一細胞遺伝子発現データ** (マウス; Tabula Muris Senis 由来) を使って、
`fast_bn` による構造学習からエッジ重要度・ブートストラップ安定性・組織別比較・
HTML レポートまでを一通り実行する例題です。

[`../example_bulk/`](../example_bulk/) と同じ構成 (`config.sh` + 番号付きステップ +
`run_all.sh`) になっていて、中身の解析ステップも共通の [`../script/`](../script/)
を呼んでいます。違いは入力データだけです。

| | example_bulk | example_sc (これ) |
| --- | --- | --- |
| データ | 生成したダミーのバルク RNA カウント | 実データ (単一細胞; **離散化済み**) |
| ステップ 1 | 前処理 (正規化 → log → フィルタ → 離散化) | 準備のみ (離散化済みなので変換しない) |
| 群 (条件) | 投薬群 (Control / TreatA / ...) | **組織** (Aorta / Liver / ...) |
| 最終評価 | 真の DAG と比較 (P/R/F1, SHD) | 対数尤度の比較 (正解構造が無いため) |

## 前提

* このディレクトリをカレントディレクトリにして実行します
* `fast_bn` をビルド済みであること (リポジトリ直下で `./compile.sh`)
* Python 3 + `networkx` / `matplotlib` (可視化に使用)

## クイックスタート

```bash
NVARS=100 ./run_all.sh              # まずはこちら: 100 遺伝子で全ステップ (数分)
./run_all.sh                        # 既定: bbknn / 2 値 / 全 2488 遺伝子 (2 時間程度)
DATASET=ss NVARS=1000 ./run_all.sh  # Smart-seq2 データ (1000 遺伝子)
```

結果は `run_<DATASET>_<DISC>_<NVARS>/` 以下 (既定では `run_bbknn_bin_all/`) に
まとまり、その中の `report.html` をブラウザで開くと図と表を一覧できます。

既定は**全遺伝子** (`NVARS=all`) です。まず流れを確認したいときは `NVARS=100` を
指定してください (下の「データの切り替え」参照)。

### 全実験の一括実行 (`bootstrap.sh`)

`run_all.sh` は 1 設定ぶんです。**時間をかけて全条件を回す**ときは
`bootstrap.sh` を使います。中身は 1 行 1 実験を並べただけです。

```bash
DATASET=bbknn DISC=bin ITERS=22000 ITERS_BS=7500  BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh
DATASET=bbknn DISC=tri ITERS=15000 ITERS_BS=5000  BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh
DATASET=ss    DISC=bin ITERS=62000 ITERS_BS=20000 BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh
DATASET=ss    DISC=tri ITERS=62000 ITERS_BS=20000 BOOTSTRAP=5 SEEDS=20 MAX_JOBS=20 ./run_all.sh
```

`nohup ./bootstrap.sh > bootstrap.log 2>&1 &` で流してください。出力は
`run_<DATASET>_<DISC>_all/` に分かれます。

#### 反復数 (`ITERS` / `ITERS_BS`) とリサンプル数の決め方

反復数は実行時間ではなく**変数数 D とサンプル数 N から**決めています。

1. Hill-Climb は成長段階では 1 反復に 1 辺しか足せないので、必要な反復数の
   下限は学習される辺の本数 `|E|` そのものです。
2. `|E| <= D x P_eff`。`P_eff = min(MAX_PARENTS, floor(log_r(N / 10)))` は
   1 ノードが実際に持てる親の数です (`r` = 1 変数の状態数 = 離散化の段階数)。
   親を `P` 個持つノードの CPT は `r^P` 通りの親設定を持ち、1 設定あたり
   10 サンプル程度は無いと罰則付きスコアはその親を保持しないため、
   **`P_eff` はサンプル数 `N` が決めます**。
3. 成長後の修正 (REMOVE / REVERSE) と Tabu が局所最適から抜ける分に成長段階の
   2 倍を見込み、`ITERS = 3 x D x P_eff`。
4. ブートストラップは `out/edges.tsv` からの warm start なので成長段階が要らず、
   `ITERS_BS = D x P_eff` (= `ITERS` の 1/3)。
5. リサンプル総数 `B = BOOTSTRAP x SEEDS` はエッジ出現確率の標準誤差
   `sqrt(p(1-p)/B) <= 0.5/sqrt(B)` から決めます。

4 ケースに当てはめた値:

| ケース | `D` | `N` | `r` | `P_eff` | `ITERS` = 3·D·P_eff | `ITERS_BS` = D·P_eff |
| --- | --- | --- | --- | --- | --- | --- |
| `bbknn` / `bin` | 2,488 | 240 | 2 | 3 | 22,392 → 22000 | 7,464 → 7500 |
| `bbknn` / `tri` | 2,488 | 240 | 3 | 2 | 14,928 → 15000 | 4,976 → 5000 |
| `ss` / `bin` | 6,862 | 4,500 | 2 | 3 | 61,758 → 62000 | 20,586 → 20000 |
| `ss` / `tri` | 6,862 | 4,500 | 3 | 3 | 61,758 → 62000 | 20,586 → 20000 |

`bbknn` / `tri` だけ反復が少ないのはサンプル数のためです。`N=240` で 3 値だと
親 3 個は `3^3 = 27` 通りの親設定に対し 8.9 サンプル/設定しかなく、スコアが
3 個目の親を保持しません (2 値なら `2^3 = 8` 通りで 30 サンプル/設定あり、親 3 個を
支えられます)。逆に `ss` は `N=4500` と多いので、離散化段階によらず
`MAX_PARENTS=3` を支えられます。リサンプル総数は `B = 5 x 20 = 100` (`SE <= 0.05`)。

> **注意**: `ss` の 2 行は上の規則どおりの反復数なので、遺伝子数が多いぶん
> 実行は数週間規模になります (1 反復のコストは辺が増えるほど上がるため)。
> 先に骨格だけ見たい場合は、`bootstrap.sh` 末尾にコメントで置いてある縮小版
> (`ITERS=6862` = `D`、平均入次数 1 まで成長させる。`ITERS_BS` はその 1/3) と
> 入れ替えてください。

## データの切り替え

`config.sh` の先頭にある 3 つのスイッチで対象データを選びます。環境変数で
上書きできるので、`config.sh` を編集せずに切り替えられます。

| 変数 | 選択肢 | 既定 | 説明 |
| --- | --- | --- | --- |
| `DATASET` | `bbknn` / `ss` | `bbknn` | どのデータセットを使うか (下表) |
| `DISC` | `bin` / `tri` | `bin` | 2 値離散化 / 3 値離散化 |
| `NVARS` | `10` / `100` / `1000` / `all` | `all` | 使用する遺伝子数 (先頭 N 列) |

| `DATASET` | ディレクトリ | 内容 | 規模 |
| --- | --- | --- | --- |
| `bbknn` | `data_bbknn_r_tissue_disc/` | 複数プロトコルを統合 (BBKNN でバッチ補正) | 240 サンプル x 2,488 遺伝子 / 24 組織 |
| `ss` | `data_ss_r_tissue_disc/` | Smart-seq2 のみ (前処理で遺伝子を絞り込み済み) | 4,500 サンプル x 6,862 遺伝子 / 45 組織 |


```bash
DATASET=ss ./run_all.sh                    # Smart-seq2 データで一式
DISC=tri NVARS=1000 ./run_all.sh           # 3 値離散化・1000 遺伝子
DATASET=ss DISC=tri NVARS=10 ./02learn.sh  # 特定のステップだけ再実行
```

**出力は組み合わせごとに別ディレクトリ** (`run_<DATASET>_<DISC>_<NVARS>/`) に出ます。
設定を変えて何度実行しても互いに上書きしません。

```
run_bbknn_bin_all/
├── data/        expr_disc.tsv (入力へのリンク), var_map.tsv, samples.tsv
├── out/         edges*.tsv, all_counts.tsv, edge_importance*.tsv, integ_*
├── bs/          ブートストラップのシード別出力
├── groups/      組織ごとの評価データ
├── figures/     学習網の図 (+ subsets/ に組織別図)
├── figures_bs/  コンセンサス網の図
└── report.html  まとめ
```

### データセットごとの規模と計算量

`NVARS=all` の意味はデータセットごとに違います。

| データセット | `NVARS=all` の変数数 | サンプル数 | 組織 |
| --- | --- | --- | --- |
| `bbknn` | 2,488 | 240 | 24 |
| `ss` | 6,862 | 4,500 | 45 |

`ss` はサンプルが 10 倍多いぶん 1 手の評価コストが高く、候補親の前処理
(相互情報量の計算) も `O(変数数 x mi-budget x サンプル数)` で効きます。実測値は
下の「実行時間の目安」を参照してください。

**反復あたりのコストは反復が進むほど増えます。** エッジが増えると REMOVE /
REVERSE の候補手が増えるためで、総時間は概ね反復数の 2 乗で伸びます
(`ss` 1,000 遺伝子・450 サンプルの旧データでは 300 反復時点 0.13 秒/反復 →
2,000 反復時点 2.2 秒/反復)。`ITERS` を倍にすると時間は 4 倍近くになるので、
大きい入力では控えめにしてください。

### 実行時間の目安

`NVARS` を上げると変数が増え、構造学習・ブートストラップ・描画の
時間が大きく伸びます。既定の `bbknn` / `bin` での実測値:

| ステップ | `NVARS=100` (100 遺伝子) | `NVARS=all` (2488 遺伝子) |
| --- | --- | --- |
| 1 準備 | 1 秒未満 | 数秒 |
| 2 構造学習 | 15 秒 | 候補親の前処理 11 秒 + **0.25 秒/反復** (既定 `ITERS=5000` で約 21 分) |
| 3 重要度 | 1 秒未満 | 1 秒未満 |
| 4 ブートストラップ | 53 秒 | **約 1〜2 時間** (`BOOTSTRAP` x `SEEDS` = 50 回) |
| 5 組織別重要度 | 30 秒 | 30 秒 |
| 6 可視化 | 1 分 | 4 分 |
| 7 レポート / 8 スコア | 1 秒未満 | 1 秒未満 |
| **合計 (`run_all.sh`)** | **約 6 分** | **2 時間前後** |

**反復あたりのコストは反復が進むほど増えます。** エッジが増えると REMOVE /
REVERSE の候補手が増えるためで、`ss` / 1000 遺伝子では 300 反復時点で
0.13 秒/反復、2000 反復時点で 2.2 秒/反復まで悪化しました。総時間は概ね
反復数の 2 乗で増えるので、`ITERS` を倍にすると時間は 4 倍近くになります。

そのため大きい入力では `ITERS` を控えめにしてください。`ss` / 1000 遺伝子を
既定の `ITERS=5000` + `BOOTSTRAP=10` x `SEEDS=5` で回すと 30 時間規模になります。

短時間で試すには反復数と本数を絞ります。

```bash
NVARS=100 ./run_all.sh                                  # 一番速い
ITERS=2000 BOOTSTRAP=3 SEEDS=2 ITERS_BS=500 ./run_all.sh # 変数を減らさず軽く
```

## 保存済みの実行例

このディレクトリには、確認済みの実行結果が 3 組入っています (どれも
`report.html` つき)。`.gitignore` されているので手元で再生成する前提です。

| ディレクトリ | データセット | 変数 | サンプル | 組織 | 図 | 所要 |
| --- | --- | --- | --- | --- | --- | --- |
| `run_bbknn_bin100/` | `bbknn` 先頭 100 遺伝子 | 100 | 240 | 24 | 63 枚 | 6 分 |
| `run_bbknn_binall/` | `bbknn` 全遺伝子 | 2,488 | 240 | 24 | 63 枚 | 7 分 |
| `run_ss_binall/` | `ss` 全遺伝子 | 6,862 | 4,500 | 45 | 105 枚 | 3 時間 50 分 |

**`run_bbknn_binall` と `run_ss_binall` は同じ縮小パラメータで回してあります**
(直接比較できるようにするため)。

```bash
# bbknn (全 2,488 遺伝子)
ITERS=300 BOOTSTRAP=2 SEEDS=2 ITERS_BS=100 MAX_JOBS=2 \
  VIZ_METRICS=dlogL VIZ_TOP_N=40 ./run_all.sh

# ss (全 6,862 遺伝子)
DATASET=ss ITERS=300 BOOTSTRAP=2 SEEDS=2 ITERS_BS=100 MAX_JOBS=2 \
  VIZ_METRICS=dlogL VIZ_TOP_N=40 ./run_all.sh
```

結果 (ステップ 8 の対数尤度):

| | エッジ数 | logL / サンプル |
| --- | --- | --- |
| `bbknn` 学習網 | 300 | -112.38 |
| `bbknn` コンセンサス網 | 318 | -111.47 |
| `ss` 学習網 | 300 | -2703.97 |
| `ss` コンセンサス網 | 358 | -2682.74 |

`ss` の logL がずっと小さいのは変数が 2.8 倍あるためで、`bbknn` との直接比較は
できません (logL は変数数とサンプル数に比例)。同じデータセット内での
学習網 vs コンセンサス網の比較として読んでください。

組織が 45 件あるので、`ss` の `figures/subsets/` と `figures_bs/subsets/` は
`bbknn` (24 組織) より図が多くなります (105 枚 vs 63 枚)。

> `ITERS=300` は既定の `ITERS=5000` よりかなり小さく、学習網のエッジ数が
> ちょうど 300 = 反復数になっている (= 全反復が ADD で、まだ辺を足し続けている)
> ことからも分かるように、**反復予算に律速された結果**です。構造の質を上げるには
> `ITERS` を増やしてください (時間は反復数の 2 乗で増えます)。

## ステップ

各スクリプトは単独で実行できます (`./03importance.sh` だけ、など)。

| # | スクリプト | 内容 | 主な出力 |
| --- | --- | --- | --- |
| 0 | `00download.sh` | データの取得・展開 (既にあればスキップ) | `data_*_r_tissue_disc/` |
| 1 | `01prepare.sh` | 行列の選択と `script/` 用の付随ファイル作成 | `data/expr_disc.tsv`, `var_map.tsv`, `samples.tsv` |
| 2 | `02learn.sh` | Hill-Climb + Tabu で構造学習 | `out/edges.tsv`, `edges_named.tsv`, `all_counts.tsv` |
| 3 | `03importance.sh` | エッジ重要度 (全サンプル) | `out/edge_importance.tsv` |
| 4 | `04bootstrap.sh` | ブートストラップ → コンセンサス網 → 重要度 | `bs/edges_seed*.tsv`, `out/integ_*` |
| 5 | `05importance_groups.sh` | **組織別**のエッジ重要度 | `out/edge_importance_g*_*.tsv` |
| 6 | `06visualize.sh` | ネットワークの描画 | `figures/`, `figures_bs/` |
| 7 | `07report.sh` | HTML レポート | `report.html` |
| 8 | `08score_check.sh` | 対数尤度の確認 (学習網 vs コンセンサス網) | `out/score_hc.tsv`, `score_bs.tsv` |

### 0. データの取得

```bash
./00download.sh          # config.sh の DATASET のぶんだけ
./00download.sh --all    # 両方のデータセット
```

`${DATA_BASE_URL}/<ディレクトリ名>.tar.gz` から取得して展開します。既にディレクトリが
ある場合は何もしません。手元にアーカイブがある場合は、このディレクトリで `tar xf`
するだけでかまいません。

各データセットの中身:

| ファイル | 内容 |
| --- | --- |
| `all_disc.tsv` / `all_disc10.tsv` / `all_disc100.tsv` / `all_disc1000.tsv` | 2 値離散化した行列 (数字は先頭 N 遺伝子) |
| `all_disc_tri*.tsv` | 同じデータを 3 値離散化したもの |
| `tissue/` , `tissue_tri/` | 組織ごとに分けた同じデータ (群ラベルの元) |

いずれも単純な TSV (1 行目=遺伝子名、以降 1 行=1 サンプル、値=離散コード) なので
`less` などで中身を確認できます。

### 1. 準備 (前処理ではない)

このデータは**すでに離散化済み**なので、`example_bulk` のような前処理
(正規化 → log → フィルタ → 離散化) は行いません。`01prepare.sh` は選んだ行列を
`script/` の汎用スクリプトが扱える形に整えるだけです。

* `data/expr_disc.tsv` — 選んだ `all_disc*.tsv` へのシンボリックリンク
* `data/var_map.tsv` — 列インデックス ↔ 遺伝子名 / 分散 / 出現水準数
* `data/samples.tsv` — 行番号 ↔ サンプル ID / **組織名**

組織名は `tissue/` 以下のファイル名から取ります。`all_disc*.tsv` の行順が
「組織別ファイルをファイル名順に連結したもの」と一致することを検証したうえで
群ラベルにしているので、ステップ 5 でそのまま組織別解析ができます。

図で強調したい遺伝子がある場合は `TARGET_GENES` に指定します
(データに無い名前は警告して無視されます)。

```bash
TARGET_GENES="Ugt1a1,Rims1" ./01prepare.sh   # -> 04_targets_highlight.png など
```

### 2〜4. 学習・重要度・ブートストラップ

`../script/` の汎用スクリプトをそのまま呼びます。パラメータは `config.sh` の
`SCORE` / `ESS` / `MAX_PARENTS` / `ITERS` / `TABU` / `TOPK` などで変更できます。

ステップ 4 が最も時間がかかります (`BOOTSTRAP` x `SEEDS` 回の構造学習)。まず
小さく試してから増やしてください。

```bash
BOOTSTRAP=3 SEEDS=2 ITERS_BS=500 ./04bootstrap.sh
```

ブートストラップは初期構造としてステップ 2 の結果 (`out/edges.tsv`) を使います
(warm start)。de novo で学習させたい場合は `WARM_START=0`。

### 5. 組織別のエッジ重要度

構造を固定したまま、評価データ (`--score-dataset`) を組織ごとに差し替えて
「どのエッジがどの組織で効いているか」を比較します。学習網とコンセンサス網の
両方について計算します。

> 1 組織あたり 10 サンプル程度と少ないため ΔlogL のノイズは大きく、絶対値も
> 全サンプル版より小さくなります (logL はサンプル数に比例)。**組織間の相対比較**
> として解釈してください。

### 6〜7. 可視化とレポート

```bash
VIZ_METRICS=dlogL,dBIC,dK2,dBDeu ./06visualize.sh
./06visualize.sh --layout kamada --top-n 80
./07report.sh --embed              # 画像を base64 埋め込みして 1 ファイルに
```

既定の `NVARS=all` では図のノードが数千個になり、全体図はほぼ塗りつぶしになります
(描画自体は通る: bbknn/all で約 4 分)。`--top-n` を絞るか、`TARGET_GENES` を指定した
図 (`04_targets_highlight.png`, `05_target_ego_*.png`) や重要度上位の部分グラフ
(`03_importance_top_*.png`) を見てください。

### 8. スコア (対数尤度) の確認

学習網とコンセンサス網の当てはまりを対数尤度で比較します。既定では**学習に使った
のと同じデータ**を評価データにしています。別のデータで評価する場合:

```bash
EVAL_INPUT=./run_bbknn_binall/groups/expr_g13_Liver.tsv ./08score_check.sh
```

> 評価データの**列 (遺伝子) の順序**は学習時の入力と一致している必要があります
> (`fast_bn` のノード番号 = 列位置)。`groups/` 以下のファイルは `data/expr_disc.tsv`
> を行方向に切り出したものなので、そのまま使えます。

## つまずきやすい点

* **未知のオプションはエラーで停止します** (`Unknown option: ...`)。古いスクリプトを
  流用するときは注意してください。なお `fast_bn` のログは既定で詳細に出ます
  (`VERBOSE=1` で `--verbose` を明示指定、静かにするなら `--quiet`)。
* **ノード番号 = 入力の列位置** です。`NVARS` や `DATASET` を変えたら、その設定で
  学習し直した `out/` を使ってください (`RUNDIR` が分かれているので通常は自動的に
  そうなります)。ステップ 5 は実行前に列順を検証し、不一致なら中断します。
* **組織別ファイルは全遺伝子ぶんの列を持っています** (`NVARS` で切られていません)。
  そのため `--score-dataset` に `tissue/*.tsv` を直接渡さないでください。
  ステップ 5 が `data/expr_disc.tsv` から切り出した `groups/*.tsv` を使います。

## 素の `fast_bn` の使い方

このディレクトリのスクリプトは `../script/` 経由で `fast_bn` を呼んでいます。
オプションの一覧や `fast_bn` を直接叩く方法は、リポジトリ直下の
[`../README.md`](../README.md) と `../fast_bn --help` を参照してください。

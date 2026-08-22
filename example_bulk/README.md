# example_bulk — バルク RNA 発現データ BN 解析の使用例 (ダミーデータ)

`../script/` の汎用ツールを使って、**バルク RNA 発現量データ (生カウント行列)**
からベイジアンネットワークを推定し、エッジ重要度・安定性・群間比較・可視化・
HTML レポートまで一通り行う実行可能な例です。データはこのディレクトリ内で
生成する**ダミーデータ**なのでダウンロード不要で、**真のネットワーク構造が既知**
なので学習精度も測れます。

`../example/` が `fast_bn` 単体の素の使い方 (既に離散化済みのデータを使う) を
示すのに対し、こちらは**生の発現量から始まる実務的なパイプライン全体**を示します。

## 前提

* このディレクトリがカレントディレクトリであること
* `fast_bn` がビルド済みであること (無ければリポジトリ直下で `./compile.sh`)
* Python 3 + `numpy` / `pandas` / `networkx` / `matplotlib`

## クイックスタート

```bash
cd FastBN/example_bulk
./run_all.sh            # 00 〜 08 を順番に実行 (手元では合計 1〜2 分)
```

終わったら `report.html` をブラウザで開くと、図・重要度テーブル・データ要約を
一覧できます (上部のボタンで dlogL / dBIC を切り替え)。

各ステップは単独で実行できます。パラメータを変えて試すときは `config.sh` を編集し、
該当ステップだけ再実行してください。

```bash
./00make_data.sh          # ダミーデータ生成
./01preprocess.sh         # 前処理 (正規化 → log → フィルタ → 離散化)
./02learn.sh              # 構造学習
./03importance.sh         # エッジ重要度
./04bootstrap.sh          # ブートストラップ → コンセンサス網
./05importance_groups.sh  # 群 (条件) 別のエッジ重要度
./06visualize.sh          # 可視化
./07report.sh             # HTML レポート
./08evaluate.sh           # 真の構造との比較 (ダミーデータのみ)
```

### 全実験の一括実行 (`bootstrap.sh`)

`run_all.sh` は 1 設定ぶんです。**時間をかけて全条件を回す**ときは
`bootstrap.sh` を使います。共通設定を `export` したうえで、データの乱数シード
(5 通り) x スコア (3) x 離散化段階 (2 / 3) x 最大親数 (2 / 3) = 60 ケースを
1 行 1 実験で並べただけのスクリプトです。

```bash
export ITERS=5000
export ITERS_BS=1500
export BOOTSTRAP=10; export SEEDS=20; export MAX_JOBS=20
...
RUNDIR=./run_s1_bdeu_b3_p2  DUMMY_SEED=1 SCORE=bdeu N_BINS=3 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./run_s1_bdeu_b3_p3  DUMMY_SEED=1 SCORE=bdeu N_BINS=3 MAX_PARENTS=3 ./run_all.sh
```

出力は `run_s<seed>_<score>_b<段階>_p<最大親数>/` に分かれます。ダミーデータは
真の DAG が既知なので、各ケースの `out/eval_hc.tsv` (学習網) と
`out/eval_bs.tsv` (コンセンサス網) を並べれば設定の良し悪しをそのまま比較できます。

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

このデータは `D` = `TOP_VAR_GENES` = 60 変数、`N` = 4 群 x 30 反復 = 120 サンプルです。

| 離散化 | `r` | `P_eff` | 必要な `ITERS` | 必要な `ITERS_BS` |
| --- | --- | --- | --- | --- |
| `N_BINS=2` | 2 | `min(3, floor(log2 12))` = 3 | 540 | 180 |
| `N_BINS=3` | 3 | `min(3, floor(log3 12))` = 2 | 360 | 120 |

真の DAG の辺が 60 本前後 (`DUMMY_EDGE_PROB=0.06` x `DUMMY_MAX_PARENTS=2`) なので
妥当な桁です。既定の `ITERS=5000` / `ITERS_BS=1500` は必要量の 9 倍 / 8 倍あり、
そのまま使っています。3 段階の方が必要反復が少ないのはサンプル数のためで、
`N=120` で 3 値だと親 3 個は `3^3 = 27` 通りの親設定に対し 4.4 サンプル/設定しか
なく、スコアが 3 個目の親を保持しません (2 値なら `2^3 = 8` 通りで 15 サンプル/設定)。

リサンプル総数は `B = 10 x 20 = 200`。`SE <= 0.035` なので `THRESHOLD_PROB=0.3` の
採否が ±0.07 (2SE) の精度で決まります。

`RUNDIR` を指定すると `data` / `out` / `bs` / `groups` / `figures` /
`report.html` がまるごとその下に出ます (既定は `.` = 従来どおりの配置)。
`config.sh` の変数はすべて環境変数で上書きできます。

## ディレクトリ構成

```text
.
├── config.sh              # 全設定 (環境変数)。自分のデータではここだけ書き換える
├── 00make_data.sh 〜 08evaluate.sh
├── run_all.sh
├── target_genes.txt       # 注目遺伝子 (00 が生成) #generated
├── data/                  #generated
│   ├── counts.tsv         # 生カウント行列 (行=遺伝子, 列=サンプル)
│   ├── sample_meta.tsv    # sample_id / group / replicate / library_size
│   ├── true_edges.tsv     # 真の DAG (評価用)
│   ├── expr_disc.tsv      # 離散化済み fast_bn 入力
│   ├── var_map.tsv        # 列インデックス ↔ 遺伝子 対応表
│   └── samples.tsv        # 行番号 ↔ サンプル / 群
├── out/                   #generated 学習結果・重要度・評価
├── bs/                    #generated ブートストラップの生出力
├── groups/                #generated 群ごとの評価データ
├── figures/               #generated 学習網の図 (+ subsets/)
├── figures_bs/            #generated コンセンサス網の図 (+ subsets/)
└── report.html            #generated
```

## ダミーデータの中身 (`00make_data.sh`)

| 項目 | 値 (`config.sh` の `DUMMY_*` で変更可) |
| --- | --- |
| 遺伝子 | 90 = 構造に参加する 60 (`G001`…) + 無相関ノイズ 30 (`NOISE001`…) |
| サンプル | 120 = 4 群 (`Control` / `TreatA` / `TreatB` / `Combo`) × 30 反復 |
| 真の DAG | 57 エッジ (最大親数 2)、親子相関 ≈ 0.8 (`DUMMY_SIGNAL_FRAC=0.65`) |
| 群効果 | 遺伝子の 25% が条件依存の発現シフトを持ち、影響は DAG 下流にも伝播 |
| 値 | 負の二項分布による**生リードカウント** (総リード数はサンプルごとにばらつく) |

ノイズ遺伝子は変動が小さいので、前処理の分散フィルタ (`TOP_VAR_GENES=60`) で
ちょうど落ちます。実際、真の構造に含まれる遺伝子は 100% 解析対象に残ります
(`08evaluate.sh` のログで確認できます)。

**実データで解析する場合**は 00 を飛ばし、`config.sh` の `EXPR_INPUT` /
`SAMPLE_META` / `ID_COL` / `NAME_COL` / `NORMALIZE` を自分のファイルに合わせて
`./01preprocess.sh` から始めてください。

## 各ステップの内容と出力

### 01 前処理

`counts.tsv` → CPM 正規化 → `log2(x+1)` → 検出率フィルタ → 分散上位 60 遺伝子 →
3 段階 (低/中/高) に等頻度離散化 → `data/expr_disc.tsv` (120 行 × 60 列)。

サンプルは `sample_meta.tsv` の群順に並べ替えられ、対応表が `data/samples.tsv` に
出ます (この表が 05 の群分割に使われます)。**ノード番号 = `expr_disc.tsv` の列位置**
なので、以後この列順が基準になります。

### 02 構造学習

Hill-Climb + Tabu で DAG を学習 (`SCORE=bdeu`, `MAX_PARENTS=2`, `ITERS=5000`)。
→ `out/edges.tsv` (62 エッジ), `out/edges_named.tsv`, `out/all_counts.tsv`。

`edges.tsv` と `edges_named.tsv` は行が 1:1 対応しており、以降の可視化・評価は
この対応から遺伝子名を復元します。

### 03 エッジ重要度

各エッジを除いたときのスコア変化 → `out/edge_importance.tsv`
(`u v ΔlogL ΔBIC ΔK2 ΔBDeu meanΔlogL_per_sample stdΔlogL_per_sample`)。
値が大きいエッジほど、そのエッジがデータの説明に効いています。

### 04 ブートストラップ (安定性)

復元抽出で 10 × 5 シード = 50 回学習し、エッジ出現頻度 (ブートストラップ確率) が
0.3 以上のエッジでコンセンサス網を構成 → `out/integ_*`。
`out/integ_edges_score.tsv` の 4 列目が各エッジの出現確率です。

### 05 群 (条件) 別のエッジ重要度

構造を固定し、評価データだけを群ごとに差し替えて重要度を再計算
→ `out/edge_importance_g{1..4}_<label>.tsv` と、コンセンサス網版
`out/integ_edge_importance_g{1..4}_<label>.tsv`。
「どのエッジがどの条件で効いているか」を比較できます。

実行前に**列順アライメント検証**が走ります (ノード番号=列位置なので、対象網の
学習入力と列順が違うと中断)。

> ΔlogL はサンプル数に比例するため、群別 (各 30 サンプル) の値は全サンプル版
> (120 サンプル) より小さくなります。**群間の相対比較**として読んでください。

### 06 可視化

`figures/` (学習網) と `figures_bs/` (コンセンサス網) に図を出力します。

| 図 | 内容 |
| --- | --- |
| `01_structure_full.png` | 全体構造 (ノードサイズ=次数, 赤=注目遺伝子) |
| `02_importance_full_<metric>.png` | エッジの色・太さ = 重要度 |
| `03_importance_top_<metric>.png` | 重要度上位エッジのみの部分グラフ |
| `04_targets_highlight.png` | 注目遺伝子 (赤) と近傍 (緑) を強調 |
| `05_target_ego_<metric>.png` | 注目遺伝子の直接の親子のみ |
| `06_bootstrap_prob.png` | エッジのブートストラップ確率 (= 安定性) |
| `subsets/subsets_grid_<metric>.png` | 4 群を共通レイアウトで並べた比較図 (群ごとに色) |
| `subsets/subsets_overlay_<metric>.png` | 全群を 1 枚に重ねた統合図 (色相=群, 濃さ=重要度) |
| `subsets/subsets_multichannel_<metric>.png` | + 線の太さ = ブートストラップ確率 |

ノード・線・文字の大きさはノード数から自動調整されます (`--style-scale` で上書き可)。

### 07 HTML レポート

図・重要度テーブル (上位 40)・データ要約を `report.html` に集約します。既定は
相対リンク参照なので軽量 (数十 KB)。単体で共有したいときは `./07report.sh --embed`。

### 08 真の構造との比較 (ダミーデータのみ)

```text
 n_true_edges                  : 57      ← 両端が解析対象に残った真のエッジ
 n_learned_edges               : 66
 undirected_precision          : 0.47    ← 骨格 (向きを無視) の一致
 undirected_recall             : 0.54
 undirected_f1                 : 0.50
 directed_f1                   : 0.23    ← 向きも含めた一致
 reversed_edges                : 17      ← 骨格は当たっているが向きが逆
```

骨格の半分程度を復元できています。有向の一致が低いのは、**観測データだけでは
同じ独立性を表す DAG (マルコフ同値類) の向きを決められない**ためで、実際に
誤りの多くは「向きが逆」(`reversed_edges`) です。`out/eval_hc_edges.tsv` に
エッジ単位の判定 (TP / FP / FP_reversed / FN) が出ます。

> 数値はサンプル数・変数数・離散化の段階数・スコア関数に強く依存します。
> `config.sh` の `DUMMY_REPLICATES` を減らす (= サンプルを減らす) と精度が
> どう落ちるか試すと、実データに必要なサンプル数の感覚がつかめます。

> **この設定ではスコアが真の構造で最大になりません。** 120 サンプルに対し
> 60 変数・3 段階離散化と情報が乏しいため、学習結果の BDeu (-7642.5) は
> 真の構造そのもの (-7871.2) より高く出ます。つまり探索が完璧でも正解には
> 戻らない領域です。指標が伸びないときは、探索の問題ではなくサンプル数・
> 離散化・スコア設定の問題を先に疑ってください。真の構造がスコア最大に
> なる条件での探索精度を見たい場合は
> [`../example_bnlearn/`](../example_bnlearn/) を参照してください。

## 主な設定 (`config.sh` 抜粋)

| 変数 | 既定 | 効果 |
| --- | --- | --- |
| `NORMALIZE` | `cpm` | 生カウントなら `cpm`。TPM/FPKM 済みなら `none` |
| `TOP_VAR_GENES` | 60 | 使用する高分散遺伝子数。多いほど網が大きく重い |
| `N_BINS` | 3 | 離散化の段階数。サンプルが少ないなら 3 程度に抑える |
| `SCORE` / `ESS` | `bdeu` / 10 | スコア関数 / BDeu 等価サンプルサイズ |
| `MAX_PARENTS` | 2 | 各ノードの最大親数。**少サンプルでは小さく** |
| `ITERS` | 5000 | Hill-Climb 反復数 |
| `BOOTSTRAP` / `SEEDS` / `MAX_JOBS` | 10 / 5 / 5 | 総リサンプル数 = 積。計算時間に直結 |
| `THRESHOLD_PROB` | 0.3 | コンセンサス採用のブートストラップ確率閾値 |
| `VIZ_METRICS` | `dlogL,dBIC` | 図を作るメトリクス (`dK2`,`dBDeu` も指定可) |

環境変数はコマンドラインからも上書きできます。

```bash
TOP_VAR_GENES=120 N_BINS=5 ./01preprocess.sh && ./02learn.sh
BOOTSTRAP=3 SEEDS=2 ./04bootstrap.sh          # 短時間で試す
VIZ_METRICS=dlogL,dBIC,dK2,dBDeu ./06visualize.sh
```

## 自分のデータに適用する

1. このディレクトリをコピー (または新しいディレクトリを作って `config.sh` だけコピー)
2. 発現量ファイルと、サンプル ID + 群ラベルの表を用意する

   ```text
   counts.tsv           sample_meta.tsv
   gene_id  gene_name  S1   S2  ...      sample_id  group
   ENSG...  ACTB       1523 987 ...      S1         Control
                                         S2         Treated
   ```
3. `config.sh` の `EXPR_INPUT` / `SAMPLE_META` / `ID_COL` / `NAME_COL` /
   `DROP_COLS` / `NORMALIZE` を合わせる (Excel なら `SHEET`, `HEADER_ROW` も)
4. `00make_data.sh` と `08evaluate.sh` は不要 (`08` は既知パスウェイを正解として
   与えれば使えます)
5. `source ./config.sh && ../script/run_pipeline.sh` で一括実行

オプションの詳細は [`../script/README.md`](../script/README.md) を参照してください。

# example_dream — DREAM チャレンジのデータによるベンチマーク

DREAM (Dialogue on Reverse Engineering Assessment and Methods) チャレンジの
ネットワーク推定データを使って、`fast_bn` の構造推定精度を測るベンチマークです。
構成は [`../example_bnlearn/`](../example_bnlearn/) と同じ (`config.sh` + 番号付き
ステップ + `run_all.sh`) で、評価も同じ `../script/evaluate_structure.py` を使います。

## 対象データセット

| `DATASETS` の値 | 内容 | 出典 | 自動取得 |
| --- | --- | --- | --- |
| `dream4` | DREAM4 in silico network challenge。人工の遺伝子制御ネットワーク size10 x 5 と size100 x 5 | [GeneNetWeaver](https://gnw.sourceforge.net/dreamchallenge.html) (約 77 MB) | ○ |
| `dream5` | DREAM5 network inference challenge。in silico / E. coli / S. cerevisiae | [Zenodo 17854236](https://zenodo.org/records/17854236) (約 38 MB, CC-BY-ND-4.0) | ○ |
| `hpn` | HPN-DREAM breast cancer network inference challenge | [Synapse syn1720047](https://www.synapse.org/#!Synapse:syn1720047) | **×** (下記) |

既定は `DATASETS="dream4 dream5"` です。

参考文献:

* DREAM4 / GeneNetWeaver: Marbach et al., *Revealing strengths and weaknesses of
  methods for gene network inference*, PNAS 107:6286–6291 (2010)
  ([PMC2963605](https://pmc.ncbi.nlm.nih.gov/articles/PMC2963605/))
* DREAM5: Marbach et al., *Wisdom of crowds for robust gene network inference*,
  Nature Methods 9:796–804 (2012)
* HPN-DREAM: Hill et al., *Inferring causal molecular networks: empirical
  assessment through a community-based effort*, Nature Methods 13:310–318 (2016)

### HPN-DREAM は自動取得できません

HPN-DREAM のデータは Synapse にあり、**アカウント登録と利用規約への同意が必要**です。
匿名アクセスはファイル本体で 403 を返します (プロジェクトの一覧表示だけは可能)。
そのため `00download.sh` は取得を試みず、手順を表示するだけにしています。

使う場合は Synapse から取得して次の 2 ファイルを用意し、`source/hpn/` に置いて
`DATASETS` に `hpn` を含めてください。無い場合は以降のステップが黙ってスキップします。

```
source/hpn/expression.tsv   行 = サンプル, 列 = 変数, 1 行目 = 変数名 (数値行列)
source/hpn/true_edges.tsv   正解エッジ (u <TAB> v)
```

> HPN-DREAM の本来の評価は「時系列＋阻害剤摂動から因果的な子孫集合を当てる」
> もので、完全な DAG を正解とする本ベンチマークの枠組みとは目的が異なります。
> 上記の形に整形できる範囲で扱う、という位置づけです。

## 評価指標と、このデータ特有の注意

`../script/evaluate_structure.py` が SHD と Directed / Skeleton の
Precision / Recall / F1 を計算します。定義は
[`../example_bnlearn/README.md`](../example_bnlearn/README.md) を参照。

このデータでは次の 2 点が bnlearn の例と大きく違います。

**1. 正解が DAG ではない** — 遺伝子制御ネットワークはフィードバックループを
含みます (DREAM4 size100 で 2〜7 組、DREAM5 で 5〜16 組の相互作用)。SID は DAG
同士でしか定義できないため自動的に **NA** になります (`true_is_dag=0` に記録)。
KL も真の CPT が無いため計算しません。

**2. gold standard が一部のペアしか判定していない** — 特に DREAM5 の gold
standard は「TF × 遺伝子」の一部ペアについてのみ 1/0 を与えており、載っていない
ペアは「不明」であって「エッジが無い」ではありません。そのまま Precision を
計算すると、判定対象外のペアを予測しただけで不正解扱いになり過小評価されます。

そこで `01prepare.sh` が `<network>_pairs.tsv` (gold standard が判定したペア) を
書き出し、`03evaluate.sh` が `--eval-pairs` でそのペアだけに絞って指標を計算します。
除外された学習エッジの本数は `n_pred_not_evaluable` 列に、絞る前の本数は
`n_pred_edges_raw` / `n_true_edges_raw` 列に残ります。

```
dream5_net3_ecoli bic   学習 529 本 -> 判定対象 153 本 (対象外 376 本)  正解 133 本
```

## 使い方

```bash
./run_all.sh                       # 既定: DREAM4 (10 網) + DREAM5 (3 網, 全遺伝子) x 3 スコア = 39 実行
DATASETS=dream4 ./run_all.sh       # DREAM4 だけ (数分)
D4_SIZES=10 ./run_all.sh           # size10 だけ
D5_MAX_VARS=300 ./run_all.sh       # DREAM5 の変数を 300 個に絞って手早く試す
SKIP_DOWNLOAD=1 ./run_all.sh
```

> **既定は DREAM5 を全遺伝子 (1643 / 4511 / 5950) で回します。**反復数も変数数に
> 応じて自動スケールするため、**DREAM5 だけで約 5.9 日**かかります (実測レートからの
> 見積もり; 下の「実測結果」参照)。まず流れを確認したいときは `DATASETS=dream4`
> (5 分) か `D5_MAX_VARS=300` を指定してください。

結果は `results/` にまとまり、`results/report.html` に集約表・グラフ・
比較図が 1 枚で出ます。

| 変数 | 既定 | 説明 |
| --- | --- | --- |
| `DATASETS` | `dream4 dream5` | 対象データセット |
| `D4_SIZES` | `10 100` | DREAM4 のネットワークサイズ |
| `D4_PARTS` | multifactorial knockouts knockdowns timeseries wildtype | 縦に結合する実験の種類 |
| `D5_NETWORKS` | `1 3 4` | DREAM5 の対象網 (2 = S. aureus は採点に未使用のため除外) |
| `D5_MAX_VARS` | `0` (全遺伝子) | DREAM5 / HPN で使う変数の上限。手早く試すなら `300` などに絞る (TF を優先し残りは分散上位) |
| `ITERS_PER_VAR` / `ITERS_MIN` | `10` / `2000` | 反復数 = `max(ITERS_MIN, ITERS_PER_VAR x 変数数)`。`ITERS` を明示すると固定値で上書き |
| `BINS` / `DISC_METHOD` | `3` / `quantile` | 離散化 |
| `SCORES` | `bic bdeu k2` | スコア関数 |
| `MAX_SID_NODES` | `150` | これを超えると SID を省略 (SID は O(p²) の d 分離判定) |
| `VIZ_MAX_NODES` | `100` | 比較図のノード数上限。超える網は次数上位のハブ部分グラフとして描く (0 で全ノード) |
| `FASTBN_BIN` | `../fast_bn` | 評価するバイナリ |

### 全実験の一括実行 (`bootstrap.sh`)

`run_all.sh` は 1 設定ぶんです。**時間をかけて全条件を回す**ときは
`bootstrap.sh` を使います。中身は 1 行 1 パスの実行を並べただけです。

```bash
export SCORES="bic bdeu k2"; export BINS=3; export DISC_METHOD="quantile"

RUNDIR=./results_dream4      DATASETS=dream4 D4_SIZES="10 100"                   ./run_all.sh
RUNDIR=./results_dream5_300  DATASETS=dream5 D5_NETWORKS="1 3 4" D5_MAX_VARS=300 ./run_all.sh
RUNDIR=./results_dream5_full DATASETS=dream5 D5_NETWORKS="1 3 4" D5_MAX_VARS=0   ./run_all.sh
#RUNDIR=./results_hpn        DATASETS=hpn    D5_MAX_VARS=0                       ./run_all.sh
```

| パス | 内容 | 反復数 (自動) | 規則からの必要量 |
| --- | --- | --- | --- |
| `results_dream4` | DREAM4 10 網 x 3 スコア = 30 実行 | 2000 (size10 / size100 とも下限) | 60 / 900 |
| `results_dream5_300` | DREAM5 を 300 変数に絞った参照値 (9 実行) | 3000 | 2700 |
| `results_dream5_full` | DREAM5 全遺伝子 (9 実行) | 16,430 / 45,110 / 59,500 | 14,787 / 40,599 / 53,550 |
| `results_hpn` | HPN-DREAM (`source/hpn` を手動配置した場合のみ) | 自動 | — |

#### 反復数 (`ITERS_PER_VAR` / `ITERS_MIN`) の決め方

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

この例題は `D` が 10 から 5950 まで 3 桁ちがうので、`ITERS` を固定せず
`ITERS = max(ITERS_MIN, ITERS_PER_VAR x D)` と自動スケールさせます。上の 3. は
`D` に比例する形なので、規則はそのまま `ITERS_PER_VAR = 3 x P_eff` と
読み替えられます。`BINS=3` で各データのサンプル数は 136 / 411 / 805 / 536 なので
`P_eff` は 2〜3、つまり `3 x P_eff = 6〜9` となり、**既定の
`ITERS_PER_VAR=10` がそのまま妥当な値**になります。そのため `bootstrap.sh` の
各行では反復数を上書きしていません (`ITERS` を渡すと自動スケールが固定値で
潰れてしまうので渡してはいけません)。

上の表のとおり全パスで必要量を満たします。実際に満たせたかは `benchmark.tsv` の
`budget_binding` 列で確認できます。`dream4` でこれが 1 になる実行がありますが、
必要量に対して反復は 2 倍以上あるので、Tabu の彷徨いによる終盤の微小更新です。

重い `dream5_full` を外したいときは、その行をコメントアウトしてください。

### 反復数は変数数に応じて自動スケールします

Hill-Climb は成長段階では 1 反復あたり 1 辺しか追加できません。そのため反復数を
固定すると、大きい網では「反復上限で打ち切られた構造」を評価してしまいます。
実際 DREAM5 の全遺伝子 (Network1 = 1643 遺伝子) を `ITERS=2000` で回すと、
2000 反復すべてが ADD で REMOVE が 1 度も起きず、正解 4012 辺に対して最大
2000 辺しか置けません。

そこで反復数はネットワークごとに次で決めます。

```
ITERS = max(ITERS_MIN, ITERS_PER_VAR x 変数数)     # 既定 max(2000, 10 x D)
```

DREAM4 は変数が 10 / 100 なので `ITERS_MIN=2000` に張り付き、従来と同じです。
DREAM5 の全遺伝子では 16,430 / 45,110 / 59,500 反復になります。

**反復予算が足りているかは指標表で確認できます。** `03evaluate.sh` が各実行について

* `iters_used` — 実際に回した反復数
* `budget_binding` — 最後の 10% の反復でまだ最良スコアが更新されたか (1 = 律速)

を記録します (`iter_state.sh` が判定)。`budget_binding=1` の行の指標は「探索性能」
ではなく「反復予算」を反映しているので、`ITERS_PER_VAR` を上げて再実行してください。

DREAM4 では全 30 実行が `budget_binding=0` で、最良スコアの更新は size100 でも
90 反復目あたりで終わっています。つまり DREAM4 の結果は反復数に律速されて
いません。

### データの作り方

* **DREAM4** — `DREAM4 training data/` の各実験 (multifactorial / knockouts /
  knockdowns / timeseries / wildtype) を**縦に結合**して 1 つの行列にします。
  size10 で 136 サンプル、size100 で 411 サンプルになります。
  時系列は `Time` 列を落として各時点を 1 サンプルとして扱います。
* **DREAM5** — 発現量行列をそのまま使います。既定 (`D5_MAX_VARS=0`) では
  全遺伝子が対象で、Network1 は 1643 遺伝子 / 805 サンプル、Network3 は
  4511 / 805、Network4 は 5950 / 536 です。`D5_MAX_VARS` に正の値を指定すると
  TF を優先しつつ残りを分散上位で絞ります。
* 連続値はいずれも既定で 3 段階の等頻度離散化 (`../script/discretize_matrix.py`
  と同じロジック) にします。値が 1 種類しかない列は落とします。

> **摂動データを観測として扱っています**: knockout / knockdown は本来「介入」
> ですが、`fast_bn` は観測データ用の構造学習器なので、ここでは全実験を単に
> 観測とみなして結合しています。介入情報を使えば向きの決定精度は上がるはずで、
> 下の Directed F1 はその情報を使わない場合の値です。

## ステップ

| # | スクリプト | 内容 |
| --- | --- | --- |
| 0 | `00download.sh` | DREAM4 / DREAM5 の zip を取得 (途中で切れるので再開しながら) |
| 1 | `01prepare.sh` → `prepare_dream.py` | 発現量の結合・離散化、正解エッジと評価対象ペアの書き出し |
| 2 | `02learn.sh` | `../script/learn_structure.sh` で構造学習 |
| 3 | `03evaluate.sh` | `../script/evaluate_structure.py` で評価 (`--eval-pairs` つき) |
| 4 | `04summarize.sh` | 集約表とグラフ |
| 5 | `05compare.sh` | 正解 vs 学習の比較図。`VIZ_MAX_NODES` を超える網は次数上位のハブ部分グラフ |
| 6 | `06report.sh` | HTML レポート |

## 実測結果

### DREAM4 (既定設定で実測)

| データセット | 変数数 | 実行 | SHD | P (dir) | R (dir) | F1 (dir) | P (skel) | R (skel) | F1 (skel) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| size10 | 10 | 15 | 14.7 ± 4.0 | 0.365 ± 0.209 | 0.210 ± 0.128 | 0.240 ± 0.117 | 0.647 ± 0.251 | 0.412 ± 0.186 | 0.453 ± 0.118 |
| size100 | 100 | 15 | 290.0 ± 62.5 | 0.099 ± 0.030 | 0.053 ± 0.025 | 0.064 ± 0.017 | 0.183 ± 0.043 | 0.103 ± 0.049 | 0.122 ± 0.033 |

全 30 実行が `budget_binding=0` (反復予算は十分)。学習は 187 秒、評価・作図・
レポートまで含めて 5 分程度です。

### DREAM5 (全遺伝子 = 既定) は未計測です

既定の全遺伝子 + 自動スケールでは**計算量が数日規模**になるため、この README の
時点では計測していません。net1 で実測した反復レートからの見積もり:

| 網 | 変数数 | サンプル | ITERS (自動) | 秒/反復 | 1 実行 | 3 スコア |
| --- | --- | --- | --- | --- | --- | --- |
| Network1 (in silico) | 1,643 | 805 | 16,430 | 0.48 | 2.2 h | 6.6 h |
| Network3 (E. coli) | 4,511 | 805 | 45,110 | 1.32 (推定) | 16.6 h | 49.7 h |
| Network4 (S. cerevisiae) | 5,950 | 536 | 59,500 | 1.74 (推定) | 28.8 h | 86.4 h |
| **合計** | | | | | | **約 5.9 日** |

秒/反復は net1 の実測 (1288 反復 / 620 秒) を変数数に比例させた値です。実行する
場合は次のようにバックグラウンドで回してください。

```bash
nohup ./02learn.sh > learn.log 2>&1 &        # 数日かかる
./03evaluate.sh && ./04summarize.sh && ./05compare.sh && ./06report.sh
```

短時間で済ませたい場合の選択肢:

```bash
D5_NETWORKS=1 ./run_all.sh              # in silico だけ (3 スコアで約 7 時間)
SCORES=bic ./run_all.sh                 # スコア比較をやめて 1/3 の時間に
ITERS_PER_VAR=2 ./run_all.sh            # 反復を 1/5 に (budget_binding=1 になる見込み)
D5_MAX_VARS=300 ./run_all.sh            # 変数を 300 個に絞る (下の参考値と同条件)
```

### 参考: DREAM5 を 300 変数に絞った場合 (`D5_MAX_VARS=300`)

全遺伝子版の前に計測した値です (この条件では `ITERS=3000`)。

| 網 | SHD | P (skel) | R (skel) | F1 (skel) | F1 (dir) |
| --- | --- | --- | --- | --- | --- |
| Network1 (in silico) | 709.3 ± 137.9 | 0.328 ± 0.086 | 0.216 ± 0.031 | 0.253 ± 0.012 | 0.155 ± 0.002 |
| Network3 (E. coli) | 270.3 ± 39.3 | 0.066 ± 0.006 | 0.090 ± 0.014 | 0.075 ± 0.002 | 0.066 ± 0.004 |
| Network4 (S. cerevisiae) | 500.0 ± 52.8 | 0.039 ± 0.005 | 0.028 ± 0.009 | 0.032 ± 0.007 | 0.009 ± 0.004 |

読み方:

* **規模が大きくなるほど急激に難しくなります。** size10 (骨格 F1 0.45) →
  size100 (0.12) → DREAM5 300 変数 (0.12)。変数が増えると候補エッジ数は p² で
  増える一方サンプル数はほぼ一定なので、当然の傾向です。
* **in silico > 実測** という並びも DREAM チャレンジ全体の傾向と一致します
  (Network1 の骨格 F1 0.25 に対し E. coli 0.075、S. cerevisiae 0.032)。実測データの
  gold standard は既知の制御関係しか含まず不完全で、かつ発現量から制御関係を
  読み取ること自体が難しいためです。DREAM5 では参加した全手法でも実測ネットワーク
  のスコアは低く、これ自体が挑戦的な課題であることが報告されています。
* **絶対値の低さは想定内**です。この例題の目的は「fast_bn が既知の難問で
  どの程度の値を出すか」を再現可能な形で記録することと、改良の前後を比較できる
  基準を用意することにあります。真の分布が既知で探索精度そのものを測りたい場合は
  [`../example_bnlearn/`](../example_bnlearn/) を使ってください。

## 補足

* **ダウンロードが途中で切れる件** — GeneNetWeaver のサーバは大きいファイルで
  接続を切ることがあります。`00download.sh` は zip として開けるようになるまで
  `curl -C -` で最大 5 回まで再開します。
* **DREAM5 の Network2 (S. aureus)** は gold standard がチャレンジの採点に
  使われなかったため、既定の `D5_NETWORKS` から外しています
  (`D5_NETWORKS="1 2 3 4"` で含められます)。
* **ライセンス** — DREAM5 の Zenodo レコードは CC-BY-ND-4.0 です。本例題は
  データを取得して評価するだけで、改変物の再配布は行いません。

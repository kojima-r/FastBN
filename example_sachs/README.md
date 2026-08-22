# example_sachs — Sachs タンパク質シグナル伝達データによるベンチマーク

Sachs らのフローサイトメトリー・データ (単一細胞レベルの 11 タンパク質・リン脂質の
測定値) を使って、`fast_bn` の構造推定を**実験的に検証済みのシグナル伝達経路**と
比較するベンチマークです。

構成は [`../example_bnlearn/`](../example_bnlearn/) と同じ (`config.sh` + 番号付き
ステップ + `run_all.sh`) で、評価も同じ `../script/evaluate_structure.py` を使います。
違いは「正解が実測に基づく生物学的経路であること」と「データが連続値なので
離散化が要ること」です。

| | example_bnlearn | example_sachs (これ) |
| --- | --- | --- |
| 正解 | 人工の BIF ネットワーク (CPT つき) | 実験で検証された経路 (エッジのみ) |
| データ | 正解 CPT からのサンプリング | 実測のフローサイトメトリー値 (連続) |
| 正解は DAG か | はい | **いいえ** (PKA ↔ PIP3 の相互作用を含む) |
| 計算できる指標 | SHD, P/R/F1 x2, SID, KL | SHD, P/R/F1 x2 (SID・KL は不可) |

## データ

* 出典: [Zenodo 7681811](https://zenodo.org/records/7681811) "Sachs: Protein and
  Phospholipids Expressions" (CC-BY-4.0)
* 元論文: Sachs et al., *Causal Protein-Signaling Networks Derived from
  Multiparameter Single-Cell Data*, Science 308:523–529 (2005)
* 変数: `Raf, Mek, Plcg, PIP2, PIP3, Erk, Akt, PKA, PKC, P38, Jnk` の 11 個
* 実験条件: 14 種 (一般刺激 + 各種阻害剤)、合計 11,672 細胞
* 正解: `GroundTruth.csv` の 20 エッジ

`00download.sh` が Zenodo から取得して `source/` に展開します (再取得はスキップ)。

### 条件セット (`PRESETS`)

| 値 | 内容 | 細胞数 |
| --- | --- | --- |
| `obs` | 一般刺激のみ (`cd3cd28`)。介入なしの観測データに相当 | 853 |
| `int` | 阻害剤などの介入条件のみ (13 条件) | 10,819 |
| `all` | 14 条件すべてを結合 | 11,672 |

> **介入データの扱いについて**: `int` / `all` は阻害剤で特定のタンパク質を
> 固定した条件を含みます。本来こうしたデータは「介入」として明示的にモデル化
> すべきですが、`fast_bn` は観測データ用の構造学習器なので、ここでは全条件を
> 単に観測とみなして結合しています。原論文が介入情報を使って向きを決めたのに
> 対し、この設定では向きの決定に使える情報が少ない点に注意してください。

## 評価指標

`../script/evaluate_structure.py` が計算します。詳しい定義は
[`../example_bnlearn/README.md`](../example_bnlearn/README.md) を参照。

* **SHD** — 欠損 + 余分 + 向き違いの合計 (小さいほど良い)
* **Directed / Skeleton の Precision / Recall / F1**
* **SID** — 正解が DAG でないため **NA**。SID は DAG 同士でしか定義されません
  (評価スクリプトが自動判定してスキップし、`true_is_dag=0` を記録します)
* **KL** — 真の CPT が与えられていないため計算しません

## 使い方

```bash
./run_all.sh                          # 既定: 18 実行 (1 分弱)
PRESETS=obs BINS_LIST=3 ./run_all.sh  # 観測データのみ・3 段階だけ
SCORES=bic ./02learn.sh               # 特定のステップだけ再実行
SKIP_DOWNLOAD=1 ./run_all.sh
```

結果は `results/` にまとまり、`results/report.html` に集約表・グラフ・
正解経路との比較図が 1 枚で出ます。

| 変数 | 既定 | 説明 |
| --- | --- | --- |
| `PRESETS` | `obs int all` | 使う条件セット |
| `BINS_LIST` | `2 3` | 離散化の段階数 (原論文は 3 段階) |
| `SCORES` | `bic bdeu k2` | スコア関数 |
| `DISC_METHOD` | `quantile` | `quantile` (等頻度) / `uniform` (等幅) |
| `USE_LOG2` | `1` | `log2(x+1)` 変換してから離散化 (蛍光強度は裾が長い) |
| `MAX_PARENTS` | `3` | 正解の最大入次数は 3 |
| `FASTBN_BIN` | `../fast_bn` | 評価するバイナリ |

### 全実験の一括実行 (`bootstrap.sh`)

`run_all.sh` は 1 設定ぶんです。**時間をかけて全条件を回す**ときは
`bootstrap.sh` を使います。離散化を 2 / 3 / 4 段階に広げた共通設定を `export`
したうえで、離散化方法 (quantile / uniform) x 最大親数 (2 / 3 / 4) を
1 行 1 実験で並べただけのスクリプトです。

```bash
export PRESETS="obs int all"
export BINS_LIST="2 3 4"
export SCORES="bic bdeu k2"
export ITERS=3000
...
RUNDIR=./results_quantile_p3 DISC_METHOD=quantile MAX_PARENTS=3 ./run_all.sh
RUNDIR=./results_uniform_p3  DISC_METHOD=uniform  MAX_PARENTS=3 ./run_all.sh
```

1 掃引 27 実行 x 6 掃引 = 162 実行。結果は
`results_<離散化方法>_p<最大親数>/` に分かれて出ます。

#### 反復数 (`ITERS`) の決め方

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

Sachs は `D = 11` で固定、サンプルは最少の `obs` でも 853 細胞あります。
一番厳しい組み合わせ (4 段階離散化 x `obs`) でも
`P_eff = min(4, floor(log4(85.3))) = 3` なので、必要な反復は
`3 x 11 x 3 = 99` です (取りうる有向辺 110 本・正解 20 辺・最大入次数 3 と整合)。
既定の `ITERS=3000` はその 30 倍あるため、そのまま使っています。

## ステップ

| # | スクリプト | 内容 |
| --- | --- | --- |
| 0 | `00download.sh` | Zenodo から `sachs.zip` を取得して `source/` に展開 |
| 1 | `01prepare.sh` → `prepare_sachs.py` | 条件を結合 → 離散化 → fast_bn 入力 + 正解エッジ |
| 2 | `02learn.sh` | `../script/learn_structure.sh` で構造学習 |
| 3 | `03evaluate.sh` | `../script/evaluate_structure.py` で評価 |
| 4 | `04summarize.sh` | 集約表とグラフ |
| 5 | `05compare.sh` | 正解経路 vs 学習結果の比較図 |
| 6 | `06report.sh` | HTML レポート |

## 実測結果

既定設定 (18 実行) での条件セットごとの平均 ± 標準偏差 (スコアと段階数にわたる):

| 条件セット | 細胞数 | SHD | P (dir) | R (dir) | F1 (dir) | P (skel) | R (skel) | F1 (skel) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `obs` | 853 | 17.50 ± 1.38 | 0.413 ± 0.126 | 0.167 ± 0.041 | 0.236 ± 0.060 | 0.786 ± 0.064 | 0.342 ± 0.044 | 0.474 ± 0.040 |
| `int` | 10,819 | 28.33 ± 1.97 | 0.203 ± 0.037 | 0.267 ± 0.052 | 0.231 ± 0.043 | 0.439 ± 0.045 | 0.605 ± 0.073 | 0.509 ± 0.055 |
| `all` | 11,672 | 28.33 ± 1.86 | 0.209 ± 0.049 | 0.275 ± 0.052 | 0.237 ± 0.051 | 0.441 ± 0.017 | 0.614 ± 0.027 | 0.513 ± 0.017 |

読み方:

* **`obs` は Precision が高く Recall が低い** (骨格 0.79 / 0.34)。853 細胞では
  BIC/BDeu のペナルティが効いて 7〜10 本しかエッジを引かず、引いた辺はよく当たる、
  という保守的な挙動です。SHD が小さいのも単にエッジ数が少ないためで、
  「よく復元できている」わけではありません。
* **`int` / `all` は逆に Recall が高く Precision が低い** (骨格 0.61 / 0.44)。
  1 万細胞あると 25〜30 本引けるようになり、正解 20 本のうち 12 本前後を拾う
  代わりに余分な辺も増えます。
* **Directed F1 はどの設定でも 0.23 前後**で頭打ちです。観測データだけでは
  マルコフ同値な DAG を区別できないうえ、上で述べたとおり介入情報を使っていない
  ためで、この設定では想定どおりの結果です。原論文が向きまで決められたのは
  介入を明示的にモデル化したからです。

## 補足

* **フローサイトメトリー値の分布** — 蛍光強度は右に大きく裾を引くので、既定で
  `log2(x+1)` 変換してから等頻度離散化しています (`USE_LOG2=0` で無効化)。
  等頻度離散化なので変換自体は順序を変えず結果に影響しませんが、`uniform`
  (等幅) と組み合わせるときは効きます。
* **正解経路の網羅性** — `GroundTruth.csv` は「実験的に確立した経路」であって
  完全な因果グラフではありません。学習結果の「余分なエッジ」の一部は、
  正解に含まれていないだけで実在する関係かもしれない点に注意してください。
* **他の例題との違い** — 真の分布が既知の条件で探索精度そのものを測りたい場合は
  [`../example_bnlearn/`](../example_bnlearn/) を、遺伝子発現データに対する
  一般的なパイプラインは [`../example_bulk/`](../example_bulk/) /
  [`../example_sc/`](../example_sc/) を参照してください。

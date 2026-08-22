# example_bnlearn — bnlearn ベンチマークによる構造推定精度の評価

[bnlearn Bayesian Network Repository](https://www.bnlearn.com/bnrepository/discrete-small.html)
の **discrete-small (ノード数 20 未満) 全 5 ネットワーク**を正解として、`fast_bn` の
構造推定精度を測るベンチマークです。

正解ネットワークの CPT からデータを生成 → `fast_bn` で構造学習 → 正解 DAG と比較、
という流れを「ネットワーク x サンプル数 x 反復 x スコア関数」の全組み合わせで回し、
結果を 1 枚の表とグラフにまとめます。

| ネットワーク | ノード | 辺 | パラメータ | 最大入次数 | 状態空間 |
| --- | --- | --- | --- | --- | --- |
| `asia` | 8 | 8 | 18 | 2 | 256 |
| `cancer` | 5 | 4 | 10 | 2 | 32 |
| `earthquake` | 5 | 4 | 10 | 2 | 32 |
| `sachs` | 11 | 17 | 178 | 3 | 177,147 |
| `survey` | 6 | 6 | 21 | 2 | 144 |

## 評価指標

`../script/evaluate_structure.py` が以下の 5 指標を計算します
(正解 DAG が分かっていれば `fast_bn` 以外の結果にも使えます)。

| 指標 | 意味 | 良い方向 |
| --- | --- | --- |
| **Structural Hamming Distance (SHD)** | 学習 DAG を正解 DAG にするのに必要な編集数。<br>`欠損辺 + 余分な辺 + 向きだけ違う辺` (反転は 1 と数える) | 小 |
| **Directed Edge P / R / F1** | 順序つきペア `u -> v` の一致 | 大 |
| **Skeleton Edge P / R / F1** | 向きを無視した無向ペア `{u, v}` の一致 | 大 |
| **Structural Intervention Distance (SID)** | 学習 DAG の親調整で介入分布 `p(x_j \| do(x_i))` を<br>誤って計算する順序つきペア数 (0 〜 p(p-1)) | 小 |
| **KL divergence** | `KL(P_true \|\| P_learned)` を全状態列挙で厳密計算 (nat) | 小 |

### 指標の読み方

* **Directed F1 は原理的に低く出ます。** 観測データだけではマルコフ同値な DAG を
  区別できないため、骨格が完全に正しくても向きは反転しがちです。Skeleton F1 と
  並べて見てください。
* **SID は「構造の使い道」に近い評価です。** SHD が小さくても、ハブの親を取り違えて
  いると因果効果の推定を大きく誤るため SID は悪化します。SID は非対称で、
  `SID(G, G) = 0` です。
* **KL は構造とパラメータを合わせた総合的な当てはまりです。** `P_learned` は
  「学習した構造 + 学習に使ったデータから推定した CPT」で、CPT の推定には
  Dirichlet 平滑化 (`KL_ALPHA`, 既定 1.0) を使います。真の構造を与えれば
  `n -> 無限大` で KL は 0 に収束し、辺を落としすぎると 0 に収束しません
  (例: asia で空グラフの KL は total correlation の 0.743 に張り付く)。

SID の実装は Peters & Bühlmann (2015) の R パッケージ `SID` の判定規則に合わせて
あり、同規則を経路列挙で書き下した独立実装とランダム DAG 4800 ケースで一致を
確認しています。

## 使い方

```bash
./run_all.sh                                   # 既定: 300 実行 (数分)
NETWORKS=asia SAMPLE_SIZES=1000 ./run_all.sh   # 1 ネットワークだけ
SCORES=bic REPLICATES=2 ./run_all.sh           # 軽く試す
SKIP_DOWNLOAD=1 ./run_all.sh                   # 取得済みの .bif を使う
```

結果は `results/` にまとまります。**`results/report.html` をブラウザで開くと、
集約表・推移グラフ・正解ネットワークとの比較図が 1 枚で見られます。**

```
results/
├── data/                 <net>_n<N>_r<R>.tsv (生成データ), <net>_true_edges.tsv (正解)
├── out/                  <net>_n<N>_r<R>_<score>/ 学習結果 (edges.tsv ほか)
├── eval/                 実行ごとの全指標
├── figures/<net>/*.png   正解 vs 学習ネットワークの比較図
├── benchmark.tsv         全実行を 1 行ずつ (生データ)
├── summary.tsv / .md     (network, n, score) ごとの平均 ± 標準偏差
├── summary_overall.tsv   サンプル数だけで集約した全体傾向
├── summary.png           指標 x ネットワークの折れ線グラフ (横軸 = サンプル数)
└── report.html           上記をまとめた HTML レポート
```

### 全実験の一括実行 (`bootstrap.sh`)

`run_all.sh` は 1 設定ぶんです。**時間をかけて全条件を回す**ときは
`bootstrap.sh` を使います。共通設定を `export` したうえで、最大親数 2 / 3 / 4 を
1 行 1 実験で並べただけの簡単なスクリプトなので、軸を変えたいときは行を
足し引きしてください。

```bash
export NETWORKS="asia cancer earthquake sachs survey"
export SAMPLE_SIZES="100 200 500 1000 2000 5000 10000"
export REPLICATES=10
export SCORES="bic bdeu k2"
export ITERS=2000
...
RUNDIR=./results_p2 MAX_PARENTS=2 ./run_all.sh
RUNDIR=./results_p3 MAX_PARENTS=3 ./run_all.sh
RUNDIR=./results_p4 MAX_PARENTS=4 ./run_all.sh
```

1 掃引 1050 実行 x 3 掃引 = 3150 実行。結果は `results_p2/` `results_p3/`
`results_p4/` に分かれて出ます (それぞれ `report.html` / `summary.md` /
`benchmark.tsv` つき)。

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

このベンチマークで一番大きいのは sachs の `D = 11`、掃引する最大親数は 4、
サンプル数は最大 10000 なので `P_eff = 4`。必要な反復は
`3 x 11 x 4 = 132` です (取りうる有向辺が 110 本、正解の辺が最大 17 本なので
妥当な桁)。既定の `ITERS=2000` はその 15 倍あるため、そのまま使っています。

### 正解ネットワークとの比較図

`05compare.sh` が、左に正解 DAG・右に学習 DAG を**同じノード配置**で並べた図を
作ります。ノード位置は正解 DAG から決めるので、スコア関数やサンプル数を変えた
図どうしもそのまま見比べられます。学習側のエッジは判定で色分けされます。

| 色 | 意味 |
| --- | --- |
| 緑 実線 | 向きまで一致 (matched) |
| 橙 実線 | 骨格は合っているが向きが逆 (reversed) |
| 赤 実線 | 正解に無い余分なエッジ (extra) |
| 灰 破線 | 学習が見落としたエッジ (missing) |

図の副題にはその実行の SHD / F1 / SID / KL が入ります。既定では
「各ネットワーク x 各スコア x 最大サンプル数 x 反復 1」の 15 枚を描きます。

```bash
VIZ_N="100 5000" ./05compare.sh   # 少数データと多数データを見比べる
VIZ_N=all ./05compare.sh          # SAMPLE_SIZES 全部 (枚数に注意)
VIZ_SCORES=bic VIZ_REP=2 ./05compare.sh
./06report.sh --embed             # 画像を埋め込んで HTML 1 ファイルにする
```

### 切り替え

すべて環境変数です (`config.sh` を編集しなくてよい)。

| 変数 | 既定 | 説明 |
| --- | --- | --- |
| `NETWORKS` | `asia cancer earthquake sachs survey` | 対象ネットワーク |
| `SAMPLE_SIZES` | `100 500 1000 5000` | 生成するサンプル数 |
| `REPLICATES` | `5` | サンプル数ごとの反復 (乱数シード 1..R) |
| `SCORES` | `bic bdeu k2` | 比較するスコア関数 |
| `MAX_PARENTS` / `ITERS` / `TABU` / `TOPK` / `ESS` | 3 / 2000 / 10 / 20 / 1 | 学習パラメータ |
| `KL_ALPHA` | `1.0` | KL 用 CPT 推定の平滑化 |
| `SUMMARY_METRICS` | SHD, P/R/F1 (dir), P/R/F1 (skel), SID, KL | 表に載せる指標 |
| `SUMMARY_PLOT_METRICS` | SHD, F1 x2, SID, KL | グラフに描く指標 (絞る) |
| `VIZ_N` / `VIZ_SCORES` / `VIZ_REP` | 最大値 / 全部 / 1 | 比較図を描く対象 |
| `FASTBN_BIN` | `../fast_bn` | 評価するバイナリ (差し替えて比較できる) |

総実行回数 = `|NETWORKS| x |SAMPLE_SIZES| x REPLICATES x |SCORES|` (既定 300)。

## ステップ

| # | スクリプト | 内容 |
| --- | --- | --- |
| 0 | `00download.sh` | bnlearn から BIF をダウンロード (`networks/<net>.bif`) |
| 1 | `01sample.sh` | 正解 CPT から祖先サンプリングでデータ生成 + 正解エッジ出力 |
| 2 | `02learn.sh` | 全データ x 全スコアで `../script/learn_structure.sh` を実行 |
| 3 | `03evaluate.sh` | 全実行を `../script/evaluate_structure.py` で評価し 1 行ずつ追記 |
| 4 | `04summarize.sh` | `../script/summarize_benchmark.py` で集約し表とグラフを出力 |
| 5 | `05compare.sh` | `../script/plot_dag_comparison.py` で正解 vs 学習の比較図を描画 |
| 6 | `06report.sh` | `../script/make_benchmark_report.py` で HTML レポートに集約 |

`00download.sh` が取得するのは**正解のネットワーク定義 (構造 + CPT)** であって
データではありません。データは `01sample.sh` がその CPT からサンプリングして作る
ので、真の構造が既知の状態で精度を測れます。

## 評価スクリプトを単体で使う

`../script/evaluate_structure.py` は本例題に依存しません。正解 DAG が分かって
いれば、任意の学習結果を評価できます。

```bash
# 正解が BIF の場合 (5 指標すべて)
python3 ../script/evaluate_structure.py \
    --true-bif networks/asia.bif \
    --pred-edges results/out/asia_n1000_r1_bic/edges.tsv \
    --input results/data/asia_n1000_r1.tsv

# 正解がエッジ表だけの場合 (KL 以外)
python3 ../script/evaluate_structure.py \
    --true-edges results/data/asia_true_edges.tsv \
    --pred-edges out/edges.tsv --input data/expr_disc.tsv
```

エッジファイルは 1 行 1 エッジの `u<TAB>v`。整数なら列インデックス、文字列なら
変数名として自動判定します (`fast_bn` の `edges.tsv` / `edges_named.tsv` どちらでも可)。

## 実測結果

既定設定 (300 実行) を実行したときの、サンプル数ごとの平均 ± 標準偏差です。

Directed (向きまで一致):

| n | runs | SHD | Precision | Recall | F1 |
| --- | --- | --- | --- | --- | --- |
| 100 | 75 | 7.240 ± 3.931 | 0.256 ± 0.254 | 0.213 ± 0.208 | 0.227 ± 0.219 |
| 500 | 75 | 5.013 ± 3.160 | 0.486 ± 0.290 | 0.402 ± 0.269 | 0.425 ± 0.265 |
| 1000 | 75 | 4.400 ± 2.686 | 0.520 ± 0.254 | 0.479 ± 0.268 | 0.487 ± 0.251 |
| 5000 | 75 | 3.307 ± 2.124 | 0.613 ± 0.272 | 0.581 ± 0.263 | 0.592 ± 0.265 |

Skeleton (向きを無視):

| n | Precision | Recall | F1 | SID (正規化) | KL |
| --- | --- | --- | --- | --- | --- |
| 100 | 0.670 ± 0.333 | 0.498 ± 0.285 | 0.544 ± 0.278 | 0.541 ± 0.176 | 0.233 ± 0.293 |
| 500 | 0.912 ± 0.128 | 0.706 ± 0.274 | 0.755 ± 0.216 | 0.455 ± 0.223 | 0.066 ± 0.104 |
| 1000 | 0.906 ± 0.124 | 0.820 ± 0.209 | 0.839 ± 0.136 | 0.411 ± 0.238 | 0.031 ± 0.046 |
| 5000 | 0.951 ± 0.095 | 0.909 ± 0.136 | 0.920 ± 0.087 | 0.357 ± 0.258 | 0.006 ± 0.007 |

サンプル数に対してすべての指標が単調に改善します。骨格は n=5000 で
Precision 0.951 / Recall 0.909 とほぼ復元できている一方、**Directed F1 は 0.592
どまり**です。これは実装の限界ではなく、観測データだけではマルコフ同値な DAG を
区別できないためで、この種のベンチマークでは想定どおりの結果です
(`cancer` / `earthquake` のように 5 ノード中 4 辺しかないネットワークほど、
向きの決まらない辺の比率が高くなります)。

## このベンチマークで見つかり、修正した不具合 (増分スコアの破綻)

このベンチマークを最初に流した時点で `fast_bn.cpp` の不具合が見つかりました。
**現在は修正済み**です (上の実測結果は修正後のもの)。記録として経緯を残します。

**症状** — `fast_bn` が出力する `# learned_score` が、保存された構造の本当のスコアと
一致しませんでした。しかも理論上ありえない値になります。

```
asia (n=5000, bic):  報告 -8709.01   実際の BIC -11767.02   真の構造の BIC -11213.59
sachs (n=5000, bic): 報告 -32904.03  実際の BIC -37710.15   真の構造の BIC -36573.30
```

`asia` の n=5000 では、飽和モデル (経験分布そのもの) の対数尤度が -11112.7 なので
`-8709` はどのベイジアンネットワークでも到達不可能です。

**原因** — `deltaAdd_andBuildNewCounts()` (fast_bn.cpp:887) は追加する親 `u` の桁を
常に**最下位桁**として counts を作ります (`j2 = j * ru + xu`)。一方 `DAG::parents()` は
親を**添字の昇順**で返し、`computeCountsForNode_full()` / `JIndexCache` /
`deltaRemove_andBuildNewCounts()` はすべて昇順を前提に桁を解釈します。そのため
既存の親より小さい添字の親を追加すると、以後 `nodeCounts[v]` の桁順が正準順と
ずれ、REMOVE が別の親の桁をマージしてしまいます。誤差が `totalNow` に蓄積し、
探索は壊れたスコアを最大化することになります。

`--iters 0` (探索なし) では増分更新が走らないためスコアは正しく、実際
`--iters 0 --init <正解構造>` は正しい BIC (-11213.593792) を返します。

**修正案** — 昇順の正しい位置に桁を挿入する (REMOVE 側のマージ写像の逆写像):

```cpp
// right_u = u より添字が大きい既存の親の基数の総乗
int right_u = 1;
for (int p : g.parents(v)) if (p > u) right_u *= ds.r[p];
const int period_u = right_u * ru;
// ループ内:  const int j2 = j * ru + xu;   を次に置き換える
const int j2 = (j / right_u) * period_u + xu * right_u + (j % right_u);
```

**検証** — 修正後は、探索が報告するスコアと、同じ構造を `--iters 0` で
ゼロから計算し直したスコアが完全に一致します。5 ネットワーク x サンプル数 3 種 x
反復 2 種 x スコア 3 種 x `--max-parents` 3 種 x `--tabu` 2 種 x `--iters` 2 種 =
**1080 構成すべてで差 0.0** (修正前は同じテストの 90 構成中 58 件が不一致、
最大差 5120)。厳密 BIC との突き合わせでも全ケース一致し、探索結果は真の構造と
同等以上のスコアに到達します。

修正前 → **修正後** のベンチマーク:

| n | SHD | Precision (skel) | Recall (skel) | F1 (skel) | KL |
| --- | --- | --- | --- | --- | --- |
| 100 | 8.00 → **7.24** | 0.584 → **0.670** | 0.478 → **0.498** | 0.506 → **0.544** | 0.298 → **0.233** |
| 500 | 7.08 → **5.01** | 0.756 → **0.912** | 0.681 → **0.706** | 0.683 → **0.755** | 0.139 → **0.066** |
| 1000 | 7.29 → **4.40** | 0.736 → **0.906** | 0.767 → **0.820** | 0.734 → **0.839** | 0.095 → **0.031** |
| 5000 | 7.15 → **3.31** | 0.740 → **0.951** | 0.848 → **0.909** | 0.783 → **0.920** | 0.071 → **0.006** |

頭打ちだった **Precision (骨格) が n=5000 で 0.740 → 0.951 に回復**しており、
「壊れたスコアを最大化した結果、余分なエッジが減らなかった」という診断と
整合します。SID (正規化) も 0.397 → 0.357、Directed F1 も 0.461 → 0.592 に
改善しました。

**回帰テストとしての使い方** — 探索まわりを触ったら、この一致確認が最も手軽な
検証になります。

```bash
# 探索の報告スコア == --iters 0 での再計算スコア か確認する
../fast_bn --input results/data/asia_n5000_r1.tsv --score bic --iters 2000 \
  --tabu 10 --topk 20 --max-parents 3 --save /tmp/e.tsv --quiet | grep learned_score
../fast_bn --input results/data/asia_n5000_r1.tsv --score bic --iters 0 \
  --topk 20 --init /tmp/e.tsv --quiet | grep learned_score
```

別バイナリとの比較は `FASTBN_BIN` で差し替えられます。

```bash
FASTBN_BIN=/path/to/other_fast_bn RUNDIR=./results_other ./run_all.sh
```

## 補足

* **決定的な CPT** — `asia` の `either` は `lung OR tub` の決定的な関数です。
  決定性は忠実性 (faithfulness) を破るため、どの構造学習アルゴリズムでも
  余分な辺が出やすくなります。`asia` の SHD が他より大きいのはこの影響もあります。
* **稀な状態** — `asia` の `asia=yes` は確率 0.01 なので、n=100 では 1 度も
  出現しないことがあります (`01sample.sh` が警告します)。`fast_bn` 側の基数は
  データの最大値+1 になりますが、評価スクリプトは正解ネットワークの基数を使うので
  KL は正しく計算されます。
* **`--iters` の打ち切り** — 結果が反復数に律速していないかは
  `results/out/<run>/log_learn.txt` の `[stop] no improving move.` の有無で確認
  できます。既定の `ITERS=2000` では 300 実行中 97 実行がこの停止条件に達しました
  (残りは Tabu サーチが動き続けて打ち切られたもので、保存されるのは探索中の
  最良構造です)。

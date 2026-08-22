---
name: network-structure-evaluation
description: This skill should be used when a learned Bayesian network / DAG must be scored against a reference structure — a known pathway, a gold-standard network, a BIF file, a simulation's true DAG, or another run's output. Triggers include "how accurate is this network", "compare with the known pathway", "SHD / precision / recall / F1 / SID / KL", "benchmark the learner", "did we recover the true edges", "既知パスウェイと比較", "正解構造との一致を評価", "ベンチマークしたい", "この推定は当たっているのか". It covers the FastBN evaluation tools (compare_edges.py, evaluate_structure.py, plot_dag_comparison.py, summarize_benchmark.py) and how to read the metrics honestly.
license: MIT
---

# 推定ネットワークの評価 (正解構造との比較)

学習した DAG を、既知パスウェイ・正解ネットワーク (BIF)・シミュレーションの真の
構造・別の設定での推定結果と突き合わせて数値化する。

## 環境

```bash
eval "$(bash "${CLAUDE_PLUGIN_ROOT}/skills/expression-network-inference/scripts/fastbn_env.sh")"
```

> **スクリプトの場所**: プラグインとして導入した場合は
> `${CLAUDE_PLUGIN_ROOT}/skills/expression-network-inference/scripts/`、`skills/` を単独でコピーした
> 場合はこの `SKILL.md` と同じディレクトリの `scripts/` にある。どちらでもない場合は
> FastBN のチェックアウト内 (`${FASTBN_HOME}/skills/.../scripts/`) を使う。


## どのツールを使うか

| 状況 | ツール | 突き合わせの単位 |
| --- | --- | --- |
| 既知パスウェイ / 遺伝子名で書かれた正解 | `${BN_SCRIPTS}/compare_edges.py` | **遺伝子名** |
| 正解が `u v` のインデックス表、または BIF | `${BN_SCRIPTS}/evaluate_structure.py` | ノード番号 (列位置) |
| 正解と推定を並べて描く | `${BN_SCRIPTS}/plot_dag_comparison.py` | — |
| 条件を変えた多数の実行を集約 | `${BN_SCRIPTS}/summarize_benchmark.py` → `make_benchmark_report.py` | — |

### 遺伝子名で比較する (実データで最も多いケース)

```bash
"${PYTHON_BIN}" "${BN_SCRIPTS}/compare_edges.py" \
  --true pathway_edges.tsv \
  --edges out/edges.tsv --edges-named out/edges_named.tsv \
  --input data/expr_disc.tsv --restrict-to-analyzed \
  --out out/eval.tsv --out-edges out/eval_edges.tsv
```

* `--true` は 1 行 1 エッジ (タブ区切りの遺伝子名 2 列)。
* **`--restrict-to-analyzed` を必ず付ける。** 前処理のフィルタで落ちた遺伝子の
  エッジは原理的に復元できないので、分母から除く方が学習器の精度に近い。
  付けない値も併記すると「解析対象の選び方でどれだけ損したか」が見える。
* 出力: 有向/無向の precision・recall・F1、`reversed_edges`、SHD 相当。
  `--out-edges` はエッジ単位の判定 (TP / FP / FP_reversed / FN)。

### インデックス表 / BIF で比較する (ベンチマーク)

```bash
"${PYTHON_BIN}" "${BN_SCRIPTS}/evaluate_structure.py" \
  --true-edges true_edges.tsv   `# または --true-bif net.bif` \
  --pred-edges out/edges.tsv \
  --input data/expr_disc.tsv \
  --out eval/run1.tsv --append benchmark.tsv --extra "score=bdeu" --quiet
```

* `--true-bif` を渡すと真の CPT が分かるので **KL(P_true || P_learned)** も出る
  (全状態列挙。`--max-states` が上限、`--alpha` は学習構造の CPT 推定の平滑化。
  0 だと未観測の親設定で確率 0 が出て KL が発散しうる)。
* `--eval-pairs <pairs.tsv>` — **正解が部分的なとき (判定済みのペアだけが分かる
  gold standard) は必須。** これを渡さないと未判定のペアが偽陽性として数えられ、
  Precision が不当に低く出る。判定対象外だった推定エッジの本数は
  `n_pred_not_evaluable` 列に出る。
* `--skip-sid` / `--max-sid-nodes` — SID は重い。ノード数が上限を超える場合と
  正解が DAG でない (フィードバックループを含む) 場合は NA になる。
* `--append` で 1 行ずつ追記し、`--extra KEY=VALUE` で条件 (スコア関数・サンプル数・
  離散化・反復数など) を列として残す。この積み上げが集約の入力になる。

### 図と集約

```bash
"${PYTHON_BIN}" "${BN_SCRIPTS}/plot_dag_comparison.py" \
  --true-edges true_edges.tsv --pred-edges out/edges.tsv \
  --input data/expr_disc.tsv --out figures/compare.png --title "true vs learned"

"${PYTHON_BIN}" "${BN_SCRIPTS}/summarize_benchmark.py" \
  --input benchmark.tsv --group-by network,n,score \
  --metrics shd,f1_directed,f1_skeleton,sid_normalized,kl_divergence \
  --out summary.tsv --markdown summary.md --plot summary.png
```

`plot_dag_comparison.py` は正解と推定を**同じノード配置**で並べ、一致・逆向き・
余分・見落としを色分けする。数値表より誤りの構造 (どこがハブごと落ちているか) が
分かる。多数の実行をまとめた HTML は `make_benchmark_report.py`。

## 指標の読み方 (ここを外すと誤った結論になる)

* **向きの一致 (directed) は低く出るのが正常。** 観測データだけではマルコフ同値類の
  中の向きを決められない。骨格 (skeleton) の P/R/F1 と `reversed_edges` を必ず
  併記し、「骨格は当たっているが向きが逆」と「そもそも辺が無い」を区別する。
* **SHD は辺の本数に敏感。** 推定が疎なら SHD は小さく出る。単独で使わず
  precision/recall と一緒に読む。
* **SID** は介入分布の違いを測る (向きの誤りに厳しい)。正解が DAG でないと定義
  できない。
* **KL** は構造 + CPT の総合評価。真の CPT (BIF) が無ければ計算できない。
* **部分的な gold standard** (DREAM5 など) では `--eval-pairs` 無しの Precision は
  無意味。
* **反復予算が律速していないか**を必ず添える。
  `"${FASTBN_HOME}/example_dream/iter_state.sh" out/log_learn.txt` が
  「反復予算が律速」と言う場合、その指標は探索性能ではなく `ITERS` の不足を
  測っている。

## 指標が悪いときの切り分け (この順で)

1. **そもそも正解がスコア最大か。** 真の構造のスコアを計算して学習網と比べる:

   ```bash
   "${FASTBN_BIN}" --input data/expr_disc.tsv --score bdeu --ess 10 \
     --iters 0 --init true_edges.tsv --save-counts /dev/null
   ```

   学習網のスコアの方が高ければ、**探索は成功していて正解に戻れない**
   (サンプル数・離散化・スコア設定の問題)。`ITERS` を増やしても無駄。
   `example_bulk` はまさにこの領域にある例。
2. **反復予算** (`iter_state.sh`)。
3. **サンプル数と離散化** — `P_eff = floor(log_r(N/10))` が 1 なら親 1 個の構造しか
   学習できない (`expression-network-inference` スキルの
   `references/parameter-sizing.md`)。
4. **前処理で正解の遺伝子が落ちていないか** (`--restrict-to-analyzed` 付きと
   無しの差を見る)。
5. 最後に探索パラメータ (`SCORE` / `ESS` / `MAX_PARENTS` / `TOPK`)。

## 既製のベンチマーク (実行可能な参照実装)

| ディレクトリ | 対象 | 特徴 |
| --- | --- | --- |
| `example_bnlearn/` | bnlearn の discrete-small (asia/cancer/earthquake/sachs/survey) | 真のネットワークからサンプリングするので**正解がスコア最大**。探索性能そのものを見たいときはここ |
| `example_sachs/` | Sachs フローサイトメトリ (実データ) | 正解が PKA↔PIP3 の 2-cycle を含むので SID は NA |
| `example_dream/` | DREAM4 / DREAM5 | 正解が**巡回かつ部分的**なので `--eval-pairs` で判定済みペアに限定。反復数を `max(ITERS_MIN, ITERS_PER_VAR × D)` で自動スケール |
| `example_bulk/` | ダミーデータ (真の DAG 既知) | 前処理からの全体パイプラインの精度 |

新しいベンチマークを作るなら、これらの `config.sh` + `00`〜`06` の番号付き
スクリプトの構成をそのまま踏襲する (設定は環境変数、CWD = 解析ディレクトリ)。
BIF からのサンプリングは `${BN_SCRIPTS}/bif_io.py sample --bif net.bif --n 1000
--seed 1 --out data.tsv --out-edges true_edges.tsv`。

## 報告の書き方

正解の性質 (完全な DAG か / 巡回か / 部分的か)、評価の範囲
(`--restrict-to-analyzed` / `--eval-pairs`)、向きの不確かさ、反復予算の状態を
必ず添える。この 4 つを書かない精度の数字は解釈できない。

# パイプライン各ステージ (環境変数・出力・コスト)

すべて「CWD = 解析ディレクトリ」で `source ./config.sh` 済みを前提とする。
`${BN_SCRIPTS}` は FastBN の `script/` (= `${FASTBN_HOME}/script`)。

権威ある一覧は `${FASTBN_HOME}/script/README.md` (日本語)。ここはその運用要約。

## 実行順と分岐

```
preprocess.sh → learn_structure.sh → edge_importance.sh → bootstrap_stability.sh
              → importance_groups.sh → viz.sh / viz_bs.sh / viz_subsets.sh → make_report.sh
```

* 既に離散化済みの行列を渡された → `preprocess.sh` は使わない (`references/preprocessing.md`)。
* 群ラベル (`SAMPLE_META`) が無い → `importance_groups.sh` と `viz_subsets.sh` は不可。
* 時間が無い / 探索的な一次解析 → `bootstrap_stability.sh` を後回しにする。
  ただしコンセンサス網 (`out/integ_*`) が無いと `viz_bs.sh` と安定性の議論はできない。
* 一括ドライバ `run_pipeline.sh` はステップ ON/OFF を持つ:
  `DO_PREPROCESS` / `DO_LEARN` / `DO_IMPORTANCE` / `DO_BOOTSTRAP` / `DO_GROUPS` /
  `DO_VIZ` / `DO_REPORT` (既定すべて 1)。

## ディレクトリ (すべて環境変数で移動可)

`DATADIR=./data` `OUTDIR=./out` `BSDIR=./bs` `GROUPDIR=./groups`
`FIGDIR=./figures` `FIGDIR_BS=./figures_bs` `REPORT_HTML=./report.html`
`INPUT=${DATADIR}/expr_disc.tsv` `VARMAP=${DATADIR}/var_map.tsv`
`SAMPLES=${DATADIR}/samples.tsv` `TARGET_FILE=./target_genes.txt`

`RUNDIR` をまとめて切り替える書き方は `example_bulk/config.sh` / `example_sc/config.sh`
を参照 (設定違いの実験を並列に置ける)。

## 1. preprocess.sh

発現量行列 → 正規化 → log → フィルタ → 離散化 → `fast_bn` 入力。

必須: `EXPR_INPUT`。主な変数:
`ORIENTATION` (genes-in-rows | samples-in-rows) / `ID_COL` / `NAME_COL` /
`LENGTH_COL` (`NORMALIZE=tpm` のとき必須) / `DROP_COLS` / `SHEET` / `HEADER_ROW` /
`SAMPLE_META` / `META_SAMPLE_COL` / `META_GROUP_COL` / `GROUP_ORDER` /
`NORMALIZE` (none|cpm|tpm) / `LOG2` / `PSEUDOCOUNT` /
`MIN_DETECT_FRAC` / `DETECT_THRESHOLD` / `MIN_MEAN_LOG` /
`TOP_VAR_GENES` / `VAR_QUANTILE` / `N_BINS` / `DISC_METHOD` (quantile|uniform) /
`TARGET_FILE` / `KEEP_GENES` / `PREPROCESS_OPTS`。

出力: `data/expr_disc.tsv` (行=サンプル, 列=遺伝子, 値=離散コード。**ノード番号=列位置**)、
`data/var_map.tsv` (index / column_name / gene_id / gene_name / variance /
detected_frac / used_levels / whitelisted)、`data/samples.tsv` (row_index /
sample_id / group)。

コスト: 数秒〜数分 (I/O 律速)。`data/var_map.tsv` の `used_levels` が 1 の列
(定数列) が多い、または `detected_frac` が低い列が残っているならフィルタが緩い。

## 2. learn_structure.sh

Hill-Climb + Tabu で DAG を学習。

`SCORE` (bic|k2|bdeu, 既定 bdeu) / `ESS` (bdeu のみ) / `ALPHA` /
`MAX_PARENTS` / `MAX_CHILDREN` / `ITERS` / `TABU` /
`TOPK` / `CAND_METRIC` (mi|chi2) / `REACH` (lazy|dense) / `JINDEX_CACHE` /
`INIT` (warm start) / `VERBOSE`。

出力: `out/edges.tsv` (u v = 列インデックス)、`out/edges_named.tsv` (同じ行順で遺伝子名)、
`out/all_counts.tsv` (CPT 推定と重要度評価に使うカウント表)、`out/log_learn.txt`。

コスト: 支配項は `ITERS × D`。MI による候補親の事前選択も D が大きいと重い
(`example_sc` の 6,862 遺伝子 × 4,500 サンプルで MI だけで約 1,000 秒)。
反復あたりのコストは |E| とともに増える (REMOVE/REVERSE 候補が増える) ので、
総時間は概ね `ITERS²` で伸びる。

学習後の確認:

```bash
wc -l out/edges.tsv                             # 辺が 0 本 / D の 3 倍以上なら設定を疑う
grep "learned_score" out/log_learn.txt          # 最終スコア
grep -c "new best" out/log_learn.txt            # スコアが更新された回数
"${FASTBN_HOME}/example_dream/iter_state.sh" out/log_learn.txt  # 反復予算が律速していないか
```

## 3. edge_importance.sh

各エッジを 1 本除いたときのスコア変化。

`INPUT` (評価データ) / `INIT` / `COUNTS` / `OUT_IMP` / `SCORE_IMP` (既定 bic)。

出力列: `u v ΔlogL ΔBIC ΔK2 ΔBDeu meanΔlogL_per_sample stdΔlogL_per_sample`。
値が大きいほどそのエッジがデータの説明に効いている。

コスト: 辺の本数に比例、通常は学習より軽い。

## 4. bootstrap_stability.sh (最も重い)

A) `fast_bn --bootstrap` を `SEEDS` 個並列 (`MAX_JOBS` で制限) →
B) `compute_bs_prob.py` で統合 (`--remove-cycle` 済み) →
C) `--iters 0` でコンセンサス構造のカウント表を再計算 →
D) コンセンサス網のエッジ重要度。

`BOOTSTRAP` / `SEEDS` (総リサンプル数 = 積) / `MAX_JOBS` / `ITERS_BS` /
`THRESHOLD_PROB` / `THRESHOLD_COUNT` / `WARM_START` (既定 1: `out/edges.tsv` から開始)。

出力: `bs/edges_seed####.tsv` (u v count prob)、`out/integ_edges.tsv` (コンセンサス構造)、
`out/integ_edges_score.tsv` (u v count prob)、`out/integ_edges2.tsv` +
`out/integ_edges_named.tsv` + `out/integ_all_counts.tsv`、`out/integ_edge_importance.tsv`。

コスト: **`BOOTSTRAP × SEEDS` 回の完全な学習**。ここが総実行時間を決める。
必ず `BOOTSTRAP=3 SEEDS=2` で 1 回計測してから本番の回数を決める。
回数の決め方は `parameter-sizing.md`。

`--save-bootstrap-counts out/edges.tsv` は `out/edges_seed0004.tsv` を書く
(シードを mod 10000・4 桁ゼロ詰めで拡張子の前に挿入)。

## 5. importance_groups.sh

構造とカウント表を**固定したまま**、評価データだけを群ごとに差し替えて重要度を
再計算する (どのエッジがどの条件で効くか)。

群の決まり方 (この優先順): `SAMPLES` (既定 `./data/samples.tsv` の group 列) →
`GROUP_META` → `GROUP_LABELS` + `GROUP_SIZES`。
`MIN_GROUP_SAMPLES` (既定 2) 未満の群はスキップ。

対象網の切り替え (コンセンサス網でやる場合):

```bash
INIT=out/integ_edges.tsv COUNTS=out/integ_all_counts.tsv \
REF_EDGES=out/integ_edges2.tsv REF_NAMED=out/integ_edges_named.tsv \
OUT_PREFIX=integ_edge_importance "${BN_SCRIPTS}/importance_groups.sh"
```

出力: `groups/expr_g{N}_<label>.tsv`、`out/<OUT_PREFIX>_g{N}_<label>.tsv`。

実行前に `check_column_alignment.py` が `INPUT` の列順と対象網の参照を照合し、
不一致なら中断する。中断したら前処理をやり直した入力を渡していないか確認する。

## 6. viz*.sh

* `viz.sh` — 学習網 → `${FIGDIR}`
* `viz_bs.sh` — コンセンサス網 → `${FIGDIR_BS}` (線の太さにブートストラップ確率)
* `viz_subsets.sh` / `viz_bs_subsets.sh` — 群別比較図 → `${FIGDIR}/subsets`

引数はそのまま `visualize.py` / `viz_subsets.py` に渡る:
`--metrics dlogL,dBIC,dK2,dBDeu` / `--top-n` / `--cmap` /
`--layout spring|kamada|circular` / `--seed` / `--hub-labels` / `--style-scale`。

**`script/viz*.sh` は `VIZ_METRICS` / `VIZ_TOP_N` を読まない。** これらを読むのは
`example_*/06visualize.sh` の側で、`--metrics` / `--top-n` に変換して渡している
(`run_pipeline.sh` も渡していないので 4 メトリクス全部を描く)。直接呼ぶときは

```bash
"${BN_SCRIPTS}/viz.sh" --metrics "${VIZ_METRICS}" --top-n "${VIZ_TOP_N}"
```

メトリクス 1 つあたり全体図・上位図・群別図が増えるので、4 つ描くと 4 倍の時間に
なる (D=90 の小さな例でも 4 メトリクスで約 70 秒)。議論に使うものだけ描く。

図の一覧と読み方は `interpretation.md`。

## 7. make_report.sh

図・重要度テーブル・データ要約を 1 枚の `report.html` に集約 (存在する成果物のみ
取り込むので途中でも作れる)。`--embed` で画像を base64 埋め込み (単一ファイルで
共有可, ただしメトリクスが多いと数百 MB になるので `--metrics` で絞る)。
`--top-edges` / `--title` も指定可。

## fast_bn の 4 つのモード (フラグで排他的に決まる)

1. **Bootstrap** — `--bootstrap B` (+ `--save-bootstrap-counts`)。**学習結果は出ない。**
2. **Edge importance** — `--edge-importance` + `--score-dataset` + `--init` + `--counts`
3. **Scoring** — `--score-dataset` 単独 (既存網の対数尤度を新しいデータで評価)
4. **Search** (既定) — `--input` から Hill-Climb + Tabu

優先順はこの順。モードを間違えると「どのフラグが効くか」が黙って変わる。

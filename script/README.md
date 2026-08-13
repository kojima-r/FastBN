# script/ — バルク RNA 発現データ BN 解析の汎用ツール群

`fast_bn` を使って**一般のバルク RNA 発現量データ**からベイジアンネットワークを
推定・評価・可視化するための汎用スクリプト集です。`gssg_analysis/` (特定データ用の
実例) を任意のデータセットに適用できるよう一般化したものです。

動く使用例は [`../example_bulk/`](../example_bulk/) にあります (ダミーデータ生成つき)。

## 設計 (2 つの約束)

1. **カレントディレクトリ = 解析ディレクトリ**
   スクリプトは自分の場所へ `cd` しません。データセットごとにディレクトリを作り、
   そこから `../script/xxx.sh` を呼びます。相対パス (`./data`, `./out`, `./figures`)
   はすべて解析ディレクトリ基準です。
2. **設定は環境変数**
   コードを編集せず環境変数で振る舞いを変えます。解析ディレクトリに `config.sh`
   を置き、`source ./config.sh` してから呼ぶ運用を推奨します
   (`../example_bulk/config.sh` がテンプレート)。

```
my_analysis/            ← 解析ディレクトリ (どこに作ってもよい)
├── config.sh           ← 設定 (example_bulk/config.sh をコピーして編集)
├── counts.tsv          ← 自分の発現量データ
├── sample_meta.tsv     ← サンプル ID と群ラベル
├── target_genes.txt    ← 注目遺伝子 (任意)
├── data/  out/  bs/  groups/  figures/  figures_bs/   ← 自動生成
└── report.html                                        ← 自動生成
```

## ツール一覧

| スクリプト | 種別 | 役割 | 主な出力 |
| --- | --- | --- | --- |
| `preprocess.sh` → `preprocess_expr.py` | Bash+Py | 発現量行列 → 正規化 → log → フィルタ → 離散化 | `data/expr_disc.tsv`, `data/var_map.tsv`, `data/samples.tsv` |
| `learn_structure.sh` | Bash | Hill-Climb + Tabu で DAG 構造学習 | `out/edges.tsv`, `out/edges_named.tsv`, `out/all_counts.tsv` |
| `edge_importance.sh` | Bash | エッジ除去によるスコア変化 (重要度) | `out/edge_importance.tsv` |
| `bootstrap_stability.sh` | Bash | ブートストラップ → コンセンサス網 → 重要度 | `bs/edges_seed*.tsv`, `out/integ_*` |
| `importance_groups.sh` → `make_groups.py`, `check_column_alignment.py` | Bash+Py | 構造を固定して群別に重要度を計算 | `groups/expr_g*_*.tsv`, `out/<prefix>_g*_*.tsv` |
| `viz.sh` / `viz_bs.sh` → `visualize.py` | Bash+Py | 学習網 / コンセンサス網の描画 | `figures/*.png`, `figures_bs/*.png` |
| `viz_subsets.sh` / `viz_bs_subsets.sh` → `viz_subsets.py` | Bash+Py | 群別重要度の比較図 (共通レイアウト・群ごとの色) | `figures*/subsets/*.png` |
| `make_report.sh` → `make_report.py` | Bash+Py | 図・表・要約を 1 枚の HTML に集約 | `report.html` |
| `run_pipeline.sh` | Bash | 上記を一括実行するドライバ | 一式 |
| `make_dummy_expr.py` | Py | 真の DAG が既知のダミーデータ生成 | `counts.tsv`, `sample_meta.tsv`, `true_edges.tsv` |
| `compare_edges.py` | Py | 学習結果を既知の正解構造と比較 (P/R/F1, SHD) | 評価 TSV |
| `common.sh` | Bash | 共通設定・パス解決・ログ (各 `*.sh` が source) | — |

## 依存

* `fast_bn` バイナリ (リポジトリ直下。無ければ `../compile.sh` でビルド)。
  別の場所にある場合は `FASTBN_BIN=/path/to/fast_bn`。
* Python 3 + `numpy` / `pandas` (前処理)、`networkx` / `matplotlib` (可視化)、
  `openpyxl` (Excel 入力を使う場合)。

## 典型的な実行順

```bash
cd my_analysis
source ./config.sh

../script/preprocess.sh           # 1) 前処理
../script/learn_structure.sh      # 2) 構造学習
../script/edge_importance.sh      # 3) エッジ重要度
../script/bootstrap_stability.sh  # 4) 安定性 (コンセンサス網)
../script/importance_groups.sh    # 5) 群別重要度 (学習網)
../script/viz.sh                  # 6) 可視化
../script/viz_bs.sh
../script/viz_subsets.sh
../script/make_report.sh          # 7) HTML レポート

# または一括
../script/run_pipeline.sh
```

コンセンサス網に対する群別重要度・可視化は、対象ファイルを環境変数で差し替えます。

```bash
INIT=out/integ_edges.tsv COUNTS=out/integ_all_counts.tsv \
REF_EDGES=out/integ_edges2.tsv REF_NAMED=out/integ_edges_named.tsv \
OUT_PREFIX=integ_edge_importance ../script/importance_groups.sh
../script/viz_bs_subsets.sh
```

---

## 1. 前処理 (`preprocess.sh` / `preprocess_expr.py`)

発現量行列を fast_bn の入力形式 (ヘッダ=変数名, 各行=サンプル, 値=離散コード) に
変換します。

処理順: **読み込み → サンプル整列 → 正規化 → log 変換 → フィルタ → 離散化 → 出力**

対応する入力:

* TSV / CSV / Excel (`--format auto` が拡張子で判定)
* 行=遺伝子・列=サンプル (既定) / 行=サンプル・列=遺伝子 (`ORIENTATION`)
* 遺伝子 ID 列・シンボル列・遺伝子長列・その他の注釈列の混在 (`ID_COL` 等で指定)
* ヘッダ行の前に余分な行がある場合 (`HEADER_ROW`)

| 環境変数 | 既定 | 説明 |
| --- | --- | --- |
| `EXPR_INPUT` | (必須) | 発現量ファイル |
| `ORIENTATION` | `genes-in-rows` | 行と列の向き |
| `ID_COL` / `NAME_COL` | 先頭列 / なし | 遺伝子 ID 列 / シンボル列 (列名 or 0 始まり位置) |
| `LENGTH_COL` | なし | 遺伝子長の列 (`NORMALIZE=tpm` のとき必須) |
| `DROP_COLS` | なし | 無視する注釈列 (カンマ区切り) |
| `SHEET` / `HEADER_ROW` | 先頭 / 0 | Excel シート名 / ヘッダ行位置 |
| `SAMPLE_META` | なし | サンプル ID と群ラベルの表 (**群別解析に必要**) |
| `META_SAMPLE_COL` / `META_GROUP_COL` | 先頭列 / `group` | メタデータの列名 |
| `GROUP_ORDER` | メタデータ登場順 | 群の並び順 (カンマ区切り) |
| `NORMALIZE` | `none` | `none` / `cpm` / `tpm`。**生カウントなら `cpm`**、TPM 済みなら `none` |
| `LOG2` / `PSEUDOCOUNT` | 1 / 1.0 | `log2(x + pseudocount)` 変換 |
| `MIN_DETECT_FRAC` / `DETECT_THRESHOLD` | 0 / 0 | 検出率フィルタ (低発現遺伝子の除去) |
| `MIN_MEAN_LOG` | 0 | 平均発現によるフィルタ |
| `TOP_VAR_GENES` | 500 | 分散上位 N 遺伝子のみ使用 (0 で無効) |
| `VAR_QUANTILE` | なし | 分散の分位点で選択 (`TOP_VAR_GENES` より優先) |
| `N_BINS` / `DISC_METHOD` | 3 / `quantile` | 離散化の段階数 / `quantile`(等頻度) or `uniform`(等幅) |
| `TARGET_FILE` | `./target_genes.txt` | フィルタを免除して必ず残す注目遺伝子 |
| `PREPROCESS_OPTS` | なし | `preprocess_expr.py` への追加オプション |

出力 3 点:

* `data/expr_disc.tsv` — fast_bn 入力。**ノード番号 = この列の位置**。
* `data/var_map.tsv` — `index / column_name / gene_id / gene_name / variance /
  detected_frac / used_levels / whitelisted`。
* `data/samples.tsv` — `row_index / sample_id / group`。群別解析 (`make_groups.py`)
  が参照します。

**変数の個数と離散化の段階数の目安**: パラメータ数は
`(段階数 - 1) × 段階数^親数` で増えます。サンプル数が少ないほど
`N_BINS` と `MAX_PARENTS` を小さくしてください (例: 12 サンプルなら
`N_BINS=3, MAX_PARENTS=2` 程度)。

## 2. 構造学習 (`learn_structure.sh`)

| 環境変数 | 既定 | 説明 |
| --- | --- | --- |
| `SCORE` | `bdeu` | `bic` / `k2` / `bdeu` |
| `ESS` | 10 | BDeu の等価サンプルサイズ |
| `MAX_PARENTS` | 3 | 各ノードの最大親数 |
| `MAX_CHILDREN` | 0 (無制限) | 各ノードの最大子数 |
| `ITERS` / `TABU` | 10000 / 30 | Hill-Climb 反復数 / Tabu 禁制期間 |
| `TOPK` / `CAND_METRIC` | 20 / `mi` | 候補親の上位 K / 関連度指標 (`mi`/`chi2`) |
| `REACH` / `JINDEX_CACHE` | `lazy` / 1024 | 到達可能性チェック / 親配置キャッシュ |
| `INIT` | なし | 初期構造 (warm start) |
| `VERBOSE` | 0 | 1 で `--verbose` |

`out/edges.tsv` (インデックス) と `out/edges_named.tsv` (遺伝子名) は**行が 1:1 対応**
します。可視化・評価はこの行対応から `idx → 遺伝子名` を復元するので、
`var_map.tsv` を作り直しても図の名前がずれません。

## 3. エッジ重要度 (`edge_importance.sh`)

各エッジを除いたときのスコア変化を計算します
(`u v ΔlogL ΔBIC ΔK2 ΔBDeu meanΔlogL_per_sample stdΔlogL_per_sample`)。

| 環境変数 | 既定 | 説明 |
| --- | --- | --- |
| `INPUT` | `./data/expr_disc.tsv` | 評価データ (`--score-dataset`) |
| `INIT` / `COUNTS` | `out/edges.tsv` / `out/all_counts.tsv` | 対象ネットワーク |
| `OUT_IMP` | `out/edge_importance.tsv` | 出力 |
| `SCORE_IMP` | `bic` | 重要度評価に用いるスコア |

## 4. ブートストラップ安定性 (`bootstrap_stability.sh`)

サンプルを復元抽出して構造学習を繰り返し、エッジ出現頻度 (ブートストラップ確率)
で安定なエッジのみを残したコンセンサス網を作り、その重要度まで計算します
(A: リサンプリング → B: 統合 → C: カウント再計算 → D: 重要度)。

| 環境変数 | 既定 | 説明 |
| --- | --- | --- |
| `BOOTSTRAP` / `SEEDS` | 10 / 5 | 総リサンプル数 = 積。**計算時間に直結** |
| `MAX_JOBS` | 0 (全並列) | 同時実行プロセス数。CPU コア数に合わせる |
| `ITERS_BS` | `ITERS` or 1000 | リサンプルごとの反復数 |
| `THRESHOLD_PROB` / `THRESHOLD_COUNT` | 0.2 / 2 | コンセンサス採用の閾値 |
| `WARM_START` | 1 | `out/edges.tsv` を初期構造に使う |

## 5. 群 (条件) 別のエッジ重要度 (`importance_groups.sh`)

構造を固定したまま `--score-dataset` を群ごとに差し替え、「どのエッジがどの条件で
効くか」を比較します。群の定義は次の優先順で決まります。

1. `SAMPLES` (既定 `./data/samples.tsv`) の `group` 列 ← **通常はこれ**
2. `GROUP_META` (サンプル ID + 群の表。行順が `INPUT` と一致している前提)
3. `GROUP_LABELS` + `GROUP_SIZES` (先頭から n 件ずつ、手動指定)

| 環境変数 | 既定 | 説明 |
| --- | --- | --- |
| `INIT` / `COUNTS` | 学習網 | 対象ネットワーク |
| `REF_EDGES` / `REF_NAMED` | `out/edges.tsv` / `out/edges_named.tsv` | 列順検証の参照 |
| `OUT_PREFIX` | `edge_importance` | 出力ファイル名の接頭辞 |
| `MIN_GROUP_SAMPLES` | 2 | この件数未満の群はスキップ |

> **列順アライメント**: fast_bn のノード番号は**列の位置**です。`INPUT` の遺伝子
> 列順が対象網の学習入力と一致しないと結果が無意味になるため、実行前に
> `check_column_alignment.py` で検証し、不一致なら中断します。前処理の
> パラメータ (例 `TOP_VAR_GENES`) を変えて `out/` を作り直したら、その学習に
> 使った入力をそのまま `INPUT` に渡してください。

> **群サイズ**: ΔlogL はサンプル数に比例するため、群別の値は全サンプル版より
> 小さくなります。絶対値ではなく**群間の相対比較**として解釈してください。

## 6. 可視化 (`viz.sh` / `viz_bs.sh` / `viz_subsets.sh` / `viz_bs_subsets.sh`)

| 図 | 内容 |
| --- | --- |
| `01_structure_full.png` | 全体構造 (ノードサイズ=次数, 赤=注目遺伝子, ラベル=ハブ+注目遺伝子) |
| `02_importance_full_<metric>.png` | 全体構造。エッジの色・太さ = 重要度 |
| `03_importance_top_<metric>.png` | 重要度上位エッジの部分グラフ (全ノードにラベル) |
| `04_targets_highlight.png` | 注目遺伝子 (赤) と近傍 (緑)・関連エッジ (橙) を強調 |
| `05_target_ego_<metric>.png` | 注目遺伝子とその直接の親子のみ |
| `06_bootstrap_prob.png` | エッジのブートストラップ確率 (コンセンサス網のみ) |
| `subsets/subset_<label>_<metric>.png` | 群ごとの重要エッジ (背景に全体網を薄く) |
| `subsets/subsets_grid_<metric>.png` | 全群を共通レイアウトで並べた比較図 |
| `subsets/subsets_overlay_<metric>.png` | 全群を 1 枚に重ねた統合図 (色相=群, 濃さ=重要度) |
| `subsets/subsets_multichannel_<metric>.png` | 3 チャンネル統合図 (+ 線の太さ=ブートストラップ確率) |
| `figures*/edge_importance_named_<metric>.tsv` | 重要度降順の名前付きエッジ表 |

主なオプション (`viz.sh` 以降の引数はそのまま Python に渡ります)。

| オプション | 既定 | 説明 |
| --- | --- | --- |
| `--metrics` | `dlogL,dBIC,dK2,dBDeu` | 図を作るメトリクス (カンマ区切り) |
| `--top-n` | 60 / 40 | 強調するエッジ本数 |
| `--cmap` | `Blues` | エッジの単色グラデーション |
| `--layout` | `spring` | `spring` / `kamada` / `circular` |
| `--seed` | 42 | レイアウト乱数シード (再現性) |
| `--hub-labels` | 15 | 全体図でラベルを付けるハブ数 |
| `--style-scale` | 自動 | ノード・線・文字の倍率。既定はノード数から自動決定 (ノード数 300 で 1.0、少ないほど大きく描く) |

環境変数 `FIGDIR` (既定 `./figures`) / `FIGDIR_BS` (既定 `./figures_bs`) で出力先を
変えられます。

## 7. HTML レポート (`make_report.sh`)

図・重要度テーブル・データ要約を 1 つの `report.html` に集約します。上部のボタンで
メトリクスを切り替えられます (メトリクス非依存の図は常時表示)。存在する成果物のみ
取り込むので、途中まででもレポートは作れます。

```bash
../script/make_report.sh                       # 相対リンク参照 (軽量, 数十 KB)
../script/make_report.sh --embed               # 画像を base64 埋め込み (単一ファイル)
../script/make_report.sh --metrics dlogL,dBIC  # 含めるメトリクスを限定
../script/make_report.sh --top-edges 60 --title "私の解析"
```

`--embed` はメトリクスが多いと数百 MB になり得るので、`--metrics` で絞ってください。

## 8. ダミーデータと精度評価 (`make_dummy_expr.py` / `compare_edges.py`)

`make_dummy_expr.py` は、疎な DAG を作り、各遺伝子の変動を
「親からの寄与 (割合 `--signal-frac`) + 独立ノイズ」で構成し、負の二項分布から
リードカウントを生成します。群 (条件) 効果を持つ遺伝子も混ぜられます。
ノイズ遺伝子は変動を小さくしてあるので、分散フィルタで落ちる挙動を確認できます。

`compare_edges.py` は学習結果を正解構造と**遺伝子名**で突き合わせ、有向/無向の
precision / recall / F1、逆向きエッジ数、SHD 相当を出します。既知パスウェイを
正解として与えれば実データでも使えます。

> **向きの一致は低く出ます**: 観測データだけでは同じ独立性を表す DAG (マルコフ
> 同値類) を区別できないため、骨格 (無向) が合っていても向きは反転しがちです。
> `reversed_edges` の本数を見て判断してください。

## gssg_analysis との対応

| gssg_analysis (特定データ用) | script/ (汎用) |
| --- | --- |
| `preprocess.py` (Excel 固定フォーマット) | `preprocess_expr.py` (TSV/CSV/Excel, 列指定・正規化を一般化) |
| `split_groups.py` (4 群 x 3 反復を固定) | `make_groups.py` (メタデータ / サイズ指定) |
| `run_learn.sh`, `run_importance.sh` | `learn_structure.sh`, `edge_importance.sh` |
| `bootstrap.sh` (ドライバ) | `run_pipeline.sh` |
| `bootstrap_bs.sh` | `bootstrap_stability.sh` |
| `run_importance_groups.sh` | `importance_groups.sh` + `check_column_alignment.py` |
| `viz/*` , `make_report.py` | `visualize.py`, `viz_subsets.py`, `make_report.py` (パスを引数化) |

`fast_bn` 自体のオプションはリポジトリ直下の `../README.md`、素の使い方は
`../example/` を参照してください。

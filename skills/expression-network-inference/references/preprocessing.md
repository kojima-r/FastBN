# 入力データの判定と前処理の決め方

`scripts/inspect_matrix.py` が機械的に判定できることはそれに任せ、ここは
判断の根拠と、スクリプトでは決められない部分の扱いを書く。

## 受け付ける形

| 項目 | 対応 | 指定する変数 |
| --- | --- | --- |
| ファイル形式 | TSV / CSV / Excel (拡張子で自動判定) | `FORMAT`, `SHEET`, `HEADER_ROW` |
| 向き | 行=遺伝子・列=サンプル (既定) / 行=サンプル・列=遺伝子 | `ORIENTATION=genes-in-rows\|samples-in-rows` |
| 注釈列の混在 | 遺伝子 ID・シンボル・遺伝子長・その他 | `ID_COL`, `NAME_COL`, `LENGTH_COL`, `DROP_COLS` |
| ヘッダ前の余分な行 | あり | `HEADER_ROW` (0 始まり) |
| 群ラベル | 別表 (sample_id + group) | `SAMPLE_META`, `META_SAMPLE_COL`, `META_GROUP_COL`, `GROUP_ORDER` |

`ID_COL` / `NAME_COL` は列名でも 0 始まりの位置でも指定できる。`NAME_COL` は図の
ラベルに使われるので、シンボル列があるなら必ず指定する (ID だけの図は読めない)。

## 正規化 (`NORMALIZE`)

| データ | 設定 |
| --- | --- |
| 生リードカウント (整数, 最大値が大きい) | `NORMALIZE=cpm`, `LOG2=1` |
| TPM / FPKM / RPKM / CPM 済み | `NORMALIZE=none`, `LOG2=1` |
| 既に log 変換済み (値が概ね 0〜20 の実数) | `NORMALIZE=none`, `LOG2=0` |
| 生カウント + 遺伝子長があり長補正したい | `NORMALIZE=tpm` + `LENGTH_COL` |
| 単一細胞の正規化済み行列 (Seurat/scanpy の出力) | `NORMALIZE=none`, `LOG2=0` (既に log1p 済みが普通) |

**二重 log 変換は静かに結果を壊す** (分散が潰れて離散化の境界がずれる)。値の範囲が
0〜20 程度の実数なら既に log 済みと考え、`LOG2=0` を疑う。

ライブラリサイズがサンプル間で桁違いなら生カウントを `none` で通してはいけない
(発現量ではなく深さで離散化される)。

## フィルタ

| 変数 | 役割 | 目安 |
| --- | --- | --- |
| `MIN_DETECT_FRAC` / `DETECT_THRESHOLD` | 低発現遺伝子の除去 | バルクなら 0.5 前後、単一細胞なら 0.1〜0.3 |
| `MIN_MEAN_LOG` | 平均発現での除去 | 併用は任意 |
| `TOP_VAR_GENES` | 分散上位 N 遺伝子のみ使用 | ここで D (変数数) が決まる。計算時間と精度の主要な操作点 |
| `VAR_QUANTILE` | 分散の分位点で選択 (`TOP_VAR_GENES` より優先) | D をデータ依存で決めたいとき |
| `TARGET_FILE` / `KEEP_GENES` | フィルタを免除して必ず残す | 注目遺伝子は必ず入れる |

**フィルタで落ちた遺伝子のエッジは原理的に復元できない。** ユーザが特定の遺伝子に
関心があるなら `TARGET_FILE` (1 行 1 遺伝子名) に必ず書く。図では赤で強調され、
`05_target_ego_*.png` / `04_targets_highlight.png` の対象になる。

## 離散化 (`N_BINS` / `DISC_METHOD`)

* `DISC_METHOD=quantile` (等頻度, 既定) — 発現分布が歪んでいる普通のケース。
* `DISC_METHOD=uniform` (等幅) — 値のスケール自体に意味があるとき。
* `N_BINS` は 2 (低/高) または 3 (低/中/高) が実用範囲。**大きくすると
  パラメータ数が `(N_BINS-1) × N_BINS^親数` で爆発する。**

サンプル数が少ないほど `N_BINS` と `MAX_PARENTS` を小さくする。目安と導出は
`parameter-sizing.md`。

離散化後は `data/var_map.tsv` の `used_levels` を確認する。1 の列は定数
(情報が無い) で、多いならフィルタか離散化が不適切。

## 既に離散化済みの行列を渡された場合

`preprocess.sh` は使わない。`data/expr_disc.tsv` (行=サンプル, 列=変数, 値=0 始まりの
整数コード, ヘッダ必須, タブ区切り) をそのまま用意し、`var_map.tsv` と `samples.tsv`
を自分で作る。`example_sc/prepare_data.py` と `example_sc/01prepare.sh` がその実例
(組織名を群ラベルにしている)。

連続値行列を離散化するだけなら `script/discretize_matrix.py` が使える
(`--bins` / `--method` / `--log2` / `--max-vars` / `--keep-vars` / `--drop-constant` /
`--transpose`)。ただし `var_map.tsv` / `samples.tsv` の体裁は `preprocess_expr.py` の
出力に合わせる必要がある (群別解析とレポートがこれらを読む)。

## 前処理後に必ず確認すること

```bash
head -2 data/expr_disc.tsv | cut -c1-200      # 向きと値の形
awk -F'\t' 'NR==1{print NF" 変数"}' data/expr_disc.tsv
echo "$(( $(wc -l < data/expr_disc.tsv) - 1 )) サンプル"
cut -f7 data/var_map.tsv | sort | uniq -c      # used_levels の分布
cut -f3 data/samples.tsv | sort | uniq -c      # 群ごとのサンプル数
```

* D (変数数) と N (サンプル数) が想定どおりか。
* 群ごとのサンプル数が `MIN_GROUP_SAMPLES` 以上あるか (少ない群は群別解析で落ちる)。
* 注目遺伝子が残っているか (`grep -w <gene> data/var_map.tsv`)。

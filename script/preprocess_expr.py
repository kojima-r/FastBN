#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess_expr.py
==================
**一般のバルク RNA 発現量データ**を FastBN/`fast_bn` の入力 TSV に変換する
汎用前処理スクリプト。gssg_analysis/preprocess.py を任意のデータセットに
適用できるよう一般化したもの。

対応する入力:
  * TSV / CSV / Excel (.xlsx)
  * 行=遺伝子・列=サンプル (バルク RNA-seq のカウント行列で一般的; 既定)
    または 行=サンプル・列=遺伝子 (--orientation samples-in-rows)
  * 遺伝子 ID 列 / 遺伝子シンボル列 / 遺伝子長列 / その他注釈列の混在
  * ヘッダ行の前に余分な行がある場合 (--header-row)

処理の流れ:
  1. 読み込み        : 発現行列 (genes x samples) と遺伝子 ID/名前を取り出す
  2. サンプル整列    : --sample-meta があれば群ごとにサンプルを並べ替え/絞り込み
  3. 正規化 (任意)   : none | cpm | tpm (--length-col が必要)
  4. log 変換        : log2(x + pseudocount)
  5. フィルタ        : 低発現 / 低検出率 / 低分散 (ホワイトリストは免除)
  6. 離散化          : 遺伝子ごとに quantile (等頻度) または uniform (等幅)
  7. 出力            : 行=サンプル・列=遺伝子の TSV (fast_bn 入力) + 対応表

出力:
  --out          : fast_bn 入力 TSV (ヘッダ=遺伝子名, 各行=サンプル, 値=0..n_bins-1)
  --out-map      : 列インデックス <-> gene id / name / 統計量 の対応表 TSV
  --out-samples  : 行番号 <-> サンプル ID / 群ラベル の対応表 TSV
                   (make_groups.py が群分割に使用)

使用例:
  # 生カウント (行=遺伝子, 1列目=gene_id, 2列目=symbol) を CPM 正規化して 3 値化
  python3 preprocess_expr.py --input counts.tsv --id-col 0 --name-col 1 \
      --normalize cpm --n-bins 3 --top-var-genes 500 \
      --sample-meta sample_meta.tsv \
      --out data/expr_disc.tsv --out-map data/var_map.tsv \
      --out-samples data/samples.tsv

  # 既に TPM 正規化済みの Excel (ヘッダが 2 行目) を 5 値化
  python3 preprocess_expr.py --input expr.xlsx --sheet TPM --header-row 1 \
      --normalize none --n-bins 5 --out data/expr_disc.tsv
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# ログ
# ---------------------------------------------------------------------------
QUIET = False


def log(*args):
    if not QUIET:
        print("[preprocess]", *args, file=sys.stderr)


def warn(*args):
    print("[preprocess] 警告:", *args, file=sys.stderr)


def die(msg):
    sys.exit(f"[preprocess] エラー: {msg}")


# ---------------------------------------------------------------------------
# 読み込み
# ---------------------------------------------------------------------------
def read_table(path, fmt, sheet, header_row):
    """発現量ファイルを DataFrame (ヘッダ付き, index はデフォルト整数) で読む。"""
    if fmt == "auto":
        ext = os.path.splitext(path)[1].lower()
        if ext in (".xlsx", ".xls", ".xlsm"):
            fmt = "excel"
        elif ext == ".csv":
            fmt = "csv"
        else:
            fmt = "tsv"
    log(f"入力形式 = {fmt} ({path})")

    if fmt == "excel":
        sheet_arg = sheet if sheet is not None else 0
        df = pd.read_excel(path, sheet_name=sheet_arg, header=header_row)
    else:
        sep = "," if fmt == "csv" else "\t"
        df = pd.read_csv(path, sep=sep, header=header_row)
    df.columns = [str(c).strip() for c in df.columns]
    log(f"読み込み shape = {df.shape} (ヘッダ行={header_row})")
    return df


def resolve_col(df, spec, what):
    """列指定 (列名 または 0 始まりの位置) を実際の列名に解決する。"""
    if spec is None or spec == "":
        return None
    if spec in df.columns:
        return spec
    try:
        i = int(spec)
    except (TypeError, ValueError):
        die(f"{what} に指定した列 '{spec}' が見つかりません。"
            f" 実際の列: {list(df.columns)[:8]} ...")
    if not (-len(df.columns) <= i < len(df.columns)):
        die(f"{what} の列位置 {i} が範囲外です (列数={len(df.columns)})")
    return df.columns[i]


def to_numeric_matrix(df, sample_cols):
    """指定列を数値行列に変換する。全く数値を含まない列は除外する。"""
    keep_cols, dropped = [], []
    data = {}
    for c in sample_cols:
        col = pd.to_numeric(df[c], errors="coerce")
        if col.notna().any():
            keep_cols.append(c)
            data[c] = col
        else:
            dropped.append(c)
    if dropped:
        log(f"数値を含まない列を除外: {dropped}")
    if not keep_cols:
        die("数値のサンプル列が 1 つも見つかりませんでした。"
            " --id-col / --name-col / --orientation の指定を確認してください。")
    mat = pd.DataFrame(data, columns=keep_cols)
    n_nan = int(mat.isna().sum().sum())
    if n_nan:
        log(f"欠損値 {n_nan} 個を 0 で補完")
    return mat.fillna(0.0)


def load_expression(args):
    """(gene_id, gene_name, expr[genes x samples], gene_length or None) を返す。

    expr は DataFrame (index=0..n_genes-1, columns=サンプル ID)。
    """
    df = read_table(args.input, args.format, args.sheet, args.header_row)

    if args.orientation == "samples-in-rows":
        # 行=サンプル, 列=遺伝子 -> 転置して genes x samples に揃える
        sid_col = resolve_col(df, args.sample_id_col, "--sample-id-col")
        if sid_col is None:
            sid_col = df.columns[0]
        sample_ids = df[sid_col].astype(str).str.strip().tolist()
        gene_cols = [c for c in df.columns if c != sid_col]
        drop = {resolve_col(df, s, "--drop-cols") for s in args.drop_cols}
        gene_cols = [c for c in gene_cols if c not in drop]
        mat = to_numeric_matrix(df, gene_cols).T          # genes x samples
        mat.columns = sample_ids
        gene_names = pd.Series([str(c) for c in mat.index])
        expr = mat.reset_index(drop=True)
        log(f"転置しました: 遺伝子 {expr.shape[0]}, サンプル {expr.shape[1]}")
        return gene_names.copy(), gene_names, expr, None

    # --- 行=遺伝子, 列=サンプル (既定) ---------------------------------------
    id_col = resolve_col(df, args.id_col, "--id-col")
    name_col = resolve_col(df, args.name_col, "--name-col")
    len_col = resolve_col(df, args.length_col, "--length-col")
    if id_col is None:
        id_col = df.columns[0]

    drop = set()
    for s in args.drop_cols:
        c = resolve_col(df, s, "--drop-cols")
        if c is not None:
            drop.add(c)
    meta_cols = {c for c in (id_col, name_col, len_col) if c is not None} | drop
    sample_cols = [c for c in df.columns if c not in meta_cols]

    gene_id = df[id_col].astype(str).str.strip()
    gene_name = (df[name_col].astype(str).str.strip() if name_col is not None
                 else gene_id.copy())
    gene_len = (pd.to_numeric(df[len_col], errors="coerce")
                if len_col is not None else None)

    expr = to_numeric_matrix(df, sample_cols)
    log(f"遺伝子 {expr.shape[0]}, サンプル {expr.shape[1]}")
    log(f"サンプル列: {list(expr.columns)[:10]}"
        f"{' ...' if expr.shape[1] > 10 else ''}")
    return gene_id, gene_name, expr, gene_len


# ---------------------------------------------------------------------------
# サンプルメタデータ (群ラベル) による整列
# ---------------------------------------------------------------------------
def load_sample_meta(path, sample_col, group_col):
    """サンプルメタデータを読み、[(sample_id, group), ...] (ファイル順) を返す。"""
    sep = "," if os.path.splitext(path)[1].lower() == ".csv" else "\t"
    meta = pd.read_csv(path, sep=sep)
    meta.columns = [str(c).strip() for c in meta.columns]
    scol = resolve_col(meta, sample_col, "--meta-sample-col")
    if scol is None:
        scol = meta.columns[0]
    if group_col is None or group_col == "":
        gcol = "group" if "group" in meta.columns else (
            meta.columns[1] if len(meta.columns) > 1 else None)
    else:
        gcol = resolve_col(meta, group_col, "--meta-group-col")
    if gcol is None:
        die(f"サンプルメタデータ {path} に群ラベル列が見つかりません "
            f"(--meta-group-col で指定してください)")
    log(f"サンプルメタデータ: {path} (サンプル列='{scol}', 群列='{gcol}')")
    rows = [(str(s).strip(), sanitize_label(str(g).strip()))
            for s, g in zip(meta[scol], meta[gcol])]
    return rows


def sanitize_label(label):
    """群ラベルをファイル名に使える形に正規化する (空白・記号 -> '_')。"""
    out = []
    for ch in label:
        out.append(ch if (ch.isalnum() or ch in "-.") else "_")
    s = "".join(out).strip("_")
    return s or "NA"


def align_samples(expr, meta_rows, group_order):
    """メタデータに従ってサンプル列を並べ替え、群ラベルの配列を返す。

    群ごとにサンプルが連続して並ぶよう並べ替える (make_groups.py が
    「先頭から n 件ずつ」で分割できるようにするため)。メタデータに無い
    サンプルは除外する。
    """
    cols = list(expr.columns)
    meta_map = dict(meta_rows)
    missing = [c for c in cols if c not in meta_map]
    if missing:
        warn(f"メタデータに無いサンプル {len(missing)} 件を除外: {missing[:8]}"
             f"{' ...' if len(missing) > 8 else ''}")
    not_in_data = [s for s, _ in meta_rows if s not in set(cols)]
    if not_in_data:
        warn(f"データに無いメタデータ行 {len(not_in_data)} 件を無視: "
             f"{not_in_data[:8]}{' ...' if len(not_in_data) > 8 else ''}")

    # 群の順序: --group-order 指定があればそれ、無ければメタデータ登場順
    if group_order:
        order = [sanitize_label(g.strip()) for g in group_order.split(",") if g.strip()]
    else:
        order, seen = [], set()
        for _, g in meta_rows:
            if g not in seen:
                seen.add(g)
                order.append(g)

    new_cols, groups = [], []
    for g in order:
        for c in cols:
            if meta_map.get(c) == g:
                new_cols.append(c)
                groups.append(g)
    dangling = [c for c in cols if c in meta_map and meta_map[c] not in set(order)]
    if dangling:
        warn(f"--group-order に無い群のサンプル {len(dangling)} 件を除外: {dangling[:8]}")
    if not new_cols:
        die("メタデータと一致するサンプルがありません。サンプル ID の表記を確認してください。")

    counts = {g: groups.count(g) for g in order if g in groups}
    log(f"サンプルを群順に整列: {counts} (合計 {len(new_cols)})")
    return expr[new_cols], groups


# ---------------------------------------------------------------------------
# 正規化・離散化
# ---------------------------------------------------------------------------
def normalize(expr, method, gene_len, scale):
    """ライブラリサイズ等の正規化。expr は genes x samples。"""
    if method == "none":
        log("正規化なし (既に TPM/FPKM/CPM 等に正規化済みの入力を想定)")
        return expr
    if method == "cpm":
        total = expr.sum(axis=0)
        total = total.replace(0, np.nan)
        out = expr.div(total, axis=1) * scale
        log(f"CPM 正規化 (列和 -> {scale:g}); ライブラリサイズ "
            f"min={expr.sum(axis=0).min():.4g}, max={expr.sum(axis=0).max():.4g}")
        return out.fillna(0.0)
    if method == "tpm":
        if gene_len is None:
            die("--normalize tpm には --length-col (遺伝子長列) が必要です")
        length_kb = (gene_len.to_numpy(dtype=float) / 1000.0)
        length_kb = np.where(np.isfinite(length_kb) & (length_kb > 0), length_kb, np.nan)
        rate = expr.div(pd.Series(length_kb, index=expr.index), axis=0)
        total = rate.sum(axis=0).replace(0, np.nan)
        out = rate.div(total, axis=1) * scale
        n_bad = int(np.isnan(length_kb).sum())
        if n_bad:
            warn(f"遺伝子長が不正な {n_bad} 遺伝子は 0 として扱います")
        log(f"TPM 正規化 (遺伝子長で割ってから列和 -> {scale:g})")
        return out.fillna(0.0)
    die(f"未知の --normalize: {method}")


def quantile_discretize(values, n_bins):
    """1 遺伝子分のベクトルを分位点 (等頻度) で 0..n_bins-1 に離散化。

    タイ (同値) が多く分位点が縮退する場合、実際の水準数は n_bins 未満に
    なることがある (定数遺伝子は事前に除去済み)。
    """
    qs = np.linspace(0, 1, n_bins + 1)[1:-1]     # 内側の境界のみ
    edges = np.quantile(values, qs)
    return np.digitize(values, edges, right=False).astype(int)


def uniform_discretize(values, n_bins):
    """1 遺伝子分のベクトルを等幅 (min..max を n_bins 等分) で離散化。"""
    lo, hi = float(np.min(values)), float(np.max(values))
    if hi <= lo:
        return np.zeros_like(values, dtype=int)
    edges = np.linspace(lo, hi, n_bins + 1)[1:-1]
    return np.digitize(values, edges, right=False).astype(int)


def discretize_matrix(mat, n_bins, method, compact):
    """genes x samples の行列を遺伝子ごとに離散化し、(codes, used_levels) を返す。"""
    fn = quantile_discretize if method == "quantile" else uniform_discretize
    disc = np.empty_like(mat, dtype=int)
    used = np.empty(mat.shape[0], dtype=int)
    for i in range(mat.shape[0]):
        codes = fn(mat[i], n_bins)
        if compact:
            # 空の水準による「穴」を詰めて 0..k-1 の連番にする
            # (順序は保持。fast_bn のカーディナリティ = 実際の水準数になる)
            uniq = np.unique(codes)
            if uniq.size and (uniq.max() != uniq.size - 1):
                remap = {old: new for new, old in enumerate(uniq)}
                codes = np.array([remap[c] for c in codes], dtype=int)
        disc[i] = codes
        used[i] = len(np.unique(codes))
    return disc, used


# ---------------------------------------------------------------------------
# ホワイトリスト・列名
# ---------------------------------------------------------------------------
def load_whitelist(path, csv_list):
    """必ず残す遺伝子 (ホワイトリスト) の集合。1 行 1 遺伝子, '#' はコメント。"""
    targets = set()
    if path:
        if not os.path.exists(path):
            warn(f"ホワイトリスト {path} が見つかりません (無効化)")
        else:
            with open(path, encoding="utf-8") as fp:
                for line in fp:
                    tok = line.strip()
                    if tok and not tok.startswith("#"):
                        targets.add(tok)
            log(f"ホワイトリスト読込: {path} ({len(targets)} 件)")
    if csv_list:
        for tok in csv_list.split(","):
            tok = tok.strip()
            if tok:
                targets.add(tok)
    return targets


def make_unique_names(gene_name, gene_id):
    """列ヘッダ用の一意な変数名を作る (重複は gene id を付与して一意化)。"""
    names, seen = [], {}
    for nm, gid in zip(gene_name, gene_id):
        base = str(nm).strip()
        if base in ("", "nan", "None", "NaN"):
            base = str(gid).strip()
        base = "_".join(base.split())          # 空白・タブ・改行を除去
        if base in seen:
            seen[base] += 1
            base = f"{base}__{gid}"
            if base in seen:                   # gene id も重複している場合
                base = f"{base}_{seen[base]}"
        seen.setdefault(base, 0)
        names.append(base)
    return names


# ---------------------------------------------------------------------------
# メイン
# ---------------------------------------------------------------------------
def build_parser():
    ap = argparse.ArgumentParser(
        description="バルク RNA 発現量データを fast_bn 入力 TSV に変換する汎用前処理",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    g = ap.add_argument_group("入力")
    g.add_argument("--input", required=True, help="入力ファイル (TSV/CSV/XLSX)")
    g.add_argument("--format", default="auto", choices=["auto", "tsv", "csv", "excel"],
                   help="入力形式 (auto は拡張子から判定)")
    g.add_argument("--sheet", default=None, help="Excel のシート名 (既定=先頭シート)")
    g.add_argument("--header-row", type=int, default=0,
                   help="ヘッダ行の位置 (0 始まり)。前に余分な行がある場合に指定")
    g.add_argument("--orientation", default="genes-in-rows",
                   choices=["genes-in-rows", "samples-in-rows"],
                   help="行と列の向き")
    g.add_argument("--id-col", default=None,
                   help="遺伝子 ID の列 (列名 or 0 始まり位置, 既定=先頭列)")
    g.add_argument("--name-col", default=None,
                   help="遺伝子シンボルの列 (省略時は ID を名前として使う)")
    g.add_argument("--length-col", default=None,
                   help="遺伝子長の列 (--normalize tpm のときに必要)")
    g.add_argument("--drop-cols", default="",
                   help="無視する注釈列 (カンマ区切りの列名 or 位置)")
    g.add_argument("--sample-id-col", default=None,
                   help="--orientation samples-in-rows のときのサンプル ID 列")

    g = ap.add_argument_group("サンプルメタデータ (群ラベル)")
    g.add_argument("--sample-meta", default=None,
                   help="サンプル ID と群ラベルの表 (TSV/CSV)。指定すると"
                        "群ごとにサンプルを並べ替え、群別解析に使える")
    g.add_argument("--meta-sample-col", default=None,
                   help="メタデータのサンプル ID 列 (既定=先頭列)")
    g.add_argument("--meta-group-col", default=None,
                   help="メタデータの群ラベル列 (既定='group' 列 or 2 列目)")
    g.add_argument("--group-order", default=None,
                   help="群の並び順 (カンマ区切り)。既定はメタデータの登場順")

    g = ap.add_argument_group("正規化・変換")
    g.add_argument("--normalize", default="none", choices=["none", "cpm", "tpm"],
                   help="ライブラリサイズ正規化。生カウントなら cpm (or tpm)")
    g.add_argument("--scale", type=float, default=1e6,
                   help="正規化後の列和 (CPM/TPM の 1e6)")
    g.add_argument("--log2", dest="log2", action="store_true", default=True,
                   help="log2(x + pseudocount) 変換を行う (既定 on)")
    g.add_argument("--no-log2", dest="log2", action="store_false",
                   help="log 変換を行わない")
    g.add_argument("--pseudocount", type=float, default=1.0,
                   help="log2(x + pseudocount) の擬似カウント")

    g = ap.add_argument_group("遺伝子フィルタ")
    g.add_argument("--min-mean-log", type=float, default=0.0,
                   help="変換後の平均発現がこの値以下の遺伝子を除外 (0 で無効)")
    g.add_argument("--min-detect-frac", type=float, default=0.0,
                   help="発現量 > --detect-threshold のサンプル割合がこの値未満の"
                        "遺伝子を除外 (0 で無効)")
    g.add_argument("--detect-threshold", type=float, default=0.0,
                   help="検出とみなす閾値 (正規化後・log 変換前の値)")
    g.add_argument("--top-var-genes", type=int, default=500,
                   help="分散上位 N 遺伝子のみ残す (0 で無効)")
    g.add_argument("--var-quantile", type=float, default=None,
                   help="分散がこの分位点以上の遺伝子のみ残す (0-1)。"
                        "--top-var-genes より優先")
    g.add_argument("--keep-genes-file", default=None,
                   help="フィルタを免除して必ず残す遺伝子のリスト (1 行 1 遺伝子)")
    g.add_argument("--keep-genes", default=None,
                   help="同上をカンマ区切りで指定")

    g = ap.add_argument_group("離散化")
    g.add_argument("--n-bins", type=int, default=3, help="離散化の段階数")
    g.add_argument("--disc-method", default="quantile", choices=["quantile", "uniform"],
                   help="quantile=等頻度 (分位点), uniform=等幅")
    g.add_argument("--compact-levels", dest="compact", action="store_true", default=True,
                   help="空の水準を詰めて 0..k-1 の連番にする (既定 on)")
    g.add_argument("--no-compact-levels", dest="compact", action="store_false",
                   help="離散化コードをそのまま出力する")

    g = ap.add_argument_group("出力")
    g.add_argument("--out", required=True, help="fast_bn 入力 TSV")
    g.add_argument("--out-map", default=None, help="列インデックス <-> 遺伝子 対応表")
    g.add_argument("--out-samples", default=None,
                   help="行番号 <-> サンプル ID / 群ラベル 対応表 "
                        "(make_groups.py が使用)")
    ap.add_argument("--quiet", action="store_true", help="ログを抑制")
    return ap


def main():
    global QUIET
    args = build_parser().parse_args()
    QUIET = args.quiet
    args.drop_cols = [s.strip() for s in args.drop_cols.split(",") if s.strip()]

    if args.n_bins < 2:
        die("--n-bins は 2 以上を指定してください")

    # --- 1. 読み込み --------------------------------------------------------
    gene_id, gene_name, expr, gene_len = load_expression(args)

    # --- 2. サンプル整列 (群ラベル) ----------------------------------------
    groups = None
    if args.sample_meta:
        meta_rows = load_sample_meta(args.sample_meta, args.meta_sample_col,
                                     args.meta_group_col)
        expr, groups = align_samples(expr, meta_rows, args.group_order)

    n_samples = expr.shape[1]
    if n_samples < 3:
        warn(f"サンプル数が {n_samples} 件しかありません。"
             "推定されるネットワークの信頼性は極めて限定的です。")
    if n_samples < args.n_bins:
        warn(f"サンプル数 ({n_samples}) < --n-bins ({args.n_bins})。"
             "離散化の水準数を減らすことを検討してください。")

    # --- ホワイトリスト ----------------------------------------------------
    targets = load_whitelist(args.keep_genes_file, args.keep_genes)
    is_target = gene_name.isin(targets) | gene_id.isin(targets)
    if targets:
        found = set(gene_name[is_target]) | set(gene_id[is_target])
        missing = sorted(t for t in targets if t not in found)
        log(f"ホワイトリスト: 指定 {len(targets)} 件中 {int(is_target.sum())} 遺伝子が該当")
        if missing:
            log(f"  元データに無い指定: {', '.join(missing[:20])}"
                f"{' ...' if len(missing) > 20 else ''}")

    # --- 3. 正規化 ---------------------------------------------------------
    norm = normalize(expr, args.normalize, gene_len, args.scale)

    # 検出率フィルタは (正規化後・log 変換前の) 値で判定する
    detected_frac = (norm > args.detect_threshold).sum(axis=1) / float(n_samples)

    # --- 4. log 変換 -------------------------------------------------------
    if args.log2:
        values = np.log2(norm + args.pseudocount)
        log(f"log2(x + {args.pseudocount}) 変換完了")
    else:
        values = norm
        log("log 変換なし")

    # --- 5. フィルタ -------------------------------------------------------
    keep = pd.Series(True, index=values.index)
    if args.min_mean_log > 0:
        keep &= (values.mean(axis=1) > args.min_mean_log) | is_target.to_numpy()
        log(f"低発現フィルタ (mean > {args.min_mean_log}): {int(keep.sum())} 遺伝子が残存")
    if args.min_detect_frac > 0:
        keep &= (detected_frac >= args.min_detect_frac) | is_target.to_numpy()
        log(f"検出率フィルタ (> {args.detect_threshold} が {args.min_detect_frac:.0%} "
            f"以上のサンプル): {int(keep.sum())} 遺伝子が残存")

    # 定数 (分散 0) 遺伝子は離散化できないため必ず除外 (ホワイトリストでも例外不可)
    variance = values.var(axis=1, ddof=0)
    dropped_const = sorted(set(gene_name[is_target.to_numpy() & (variance == 0).to_numpy()]))
    if dropped_const:
        warn("ホワイトリスト指定だが全サンプルで定数のため除外: "
             + ", ".join(dropped_const[:20]))
    keep &= variance > 0

    values = values[keep.to_numpy()]
    variance = variance[keep.to_numpy()]
    gene_id = gene_id[keep.to_numpy()]
    gene_name = gene_name[keep.to_numpy()]
    is_target = is_target[keep.to_numpy()]
    detected_frac = detected_frac[keep.to_numpy()]
    log(f"定数遺伝子除去後: {values.shape[0]} 遺伝子")
    if values.shape[0] == 0:
        die("フィルタ後に遺伝子が残りませんでした。閾値を見直してください。")

    # 分散フィルタ
    if args.var_quantile is not None:
        thr = variance.quantile(args.var_quantile)
        sel = variance >= thr
        log(f"分散 >= {args.var_quantile} 分位点 ({thr:.4f}): {int(sel.sum())} 遺伝子を選択")
    elif args.top_var_genes and args.top_var_genes > 0:
        n = min(args.top_var_genes, len(variance))
        top_idx = variance.sort_values(ascending=False).index[:n]
        sel = pd.Series(variance.index.isin(top_idx), index=variance.index)
        log(f"分散上位 {n} 遺伝子を選択")
    else:
        sel = pd.Series(True, index=variance.index)
        log("分散フィルタなし (全遺伝子を使用)")

    n_added = int((is_target.to_numpy() & ~sel.to_numpy()).sum())
    sel = pd.Series(sel.to_numpy() | is_target.to_numpy(), index=variance.index)
    if n_added:
        log(f"ホワイトリストにより分散フィルタ外の {n_added} 遺伝子を追加保持")

    m = sel.to_numpy()
    values, variance = values[m], variance[m]
    gene_id, gene_name, is_target = gene_id[m], gene_name[m], is_target[m]
    detected_frac = detected_frac[m]
    log(f"最終的に使用する遺伝子数: {values.shape[0]}")

    # --- 6. 離散化 ---------------------------------------------------------
    mat = values.to_numpy()
    disc, used_levels = discretize_matrix(mat, args.n_bins, args.disc_method,
                                          args.compact)
    log(f"{args.n_bins} 段階離散化完了 ({args.disc_method}; 実水準数 "
        f"平均={used_levels.mean():.2f}, 最小={used_levels.min()}, "
        f"最大={used_levels.max()})")

    # 離散化後に 1 水準しか取らない遺伝子は BN に寄与しないため除外
    informative = used_levels >= 2
    n_drop = int((~informative).sum())
    is_target_arr = is_target.to_numpy()
    dropped_disc = sorted(set(np.asarray(gene_name)[is_target_arr & ~informative]))
    if dropped_disc:
        warn("ホワイトリスト指定だが離散化後に定数となり除外: "
             + ", ".join(dropped_disc[:20]))
    if n_drop:
        disc = disc[informative]
        gene_id, gene_name = gene_id[informative], gene_name[informative]
        variance = variance[informative]
        detected_frac = detected_frac[informative]
        used_levels = used_levels[informative]
        is_target_arr = is_target_arr[informative]
        log(f"離散化後に定数となった {n_drop} 遺伝子を除外 -> 残り {disc.shape[0]} 遺伝子")
    if disc.shape[0] == 0:
        die("離散化後に有効な遺伝子が残りませんでした。--n-bins や"
            "フィルタ条件を見直してください。")
    if targets:
        kept = sorted(set(np.asarray(gene_name)[is_target_arr]))
        log(f"ホワイトリスト最終結果: {len(kept)} 遺伝子を保持 "
            f"({', '.join(kept) if kept else 'なし'})")

    # --- 7. 出力 -----------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    col_names = make_unique_names(list(gene_name), list(gene_id))
    out_df = pd.DataFrame(disc.T, columns=col_names)      # 行=サンプル, 列=遺伝子
    out_df.to_csv(args.out, sep="\t", index=False)
    log(f"出力: {args.out} (行={out_df.shape[0]} サンプル, 列={out_df.shape[1]} 遺伝子)")

    if args.out_map:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_map)) or ".", exist_ok=True)
        pd.DataFrame({
            "index": range(len(col_names)),
            "column_name": col_names,
            "gene_id": list(gene_id),
            "gene_name": list(gene_name),
            "variance": np.asarray(variance),
            "detected_frac": np.asarray(detected_frac),
            "used_levels": used_levels,
            "whitelisted": is_target_arr.astype(int),
        }).to_csv(args.out_map, sep="\t", index=False)
        log(f"対応表を出力: {args.out_map}")

    if args.out_samples:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_samples)) or ".",
                    exist_ok=True)
        sdf = pd.DataFrame({
            "row_index": range(len(expr.columns)),
            "sample_id": [str(c) for c in expr.columns],
            "group": groups if groups is not None else ["ALL"] * len(expr.columns),
        })
        sdf.to_csv(args.out_samples, sep="\t", index=False)
        log(f"サンプル表を出力: {args.out_samples}")


if __name__ == "__main__":
    main()

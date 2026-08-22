#!/usr/bin/env python3
# =============================================================================
# inspect_matrix.py
#   受け取った発現量行列を点検し、FastBN パイプラインの設定値を推奨する。
#
#   判定するもの:
#     * 行と列の向き (行=遺伝子 / 行=サンプル)
#     * 注釈列 (遺伝子 ID / シンボル / 遺伝子長など) と数値ブロック
#     * 値の種類 (生カウント / 線形正規化済み / log 済み / 既に離散化済み)
#     * サンプル数 N・変数数 D・欠損・重複・定数列・ゼロ率
#
#   出力するもの:
#     * 上記の要約 (標準出力)
#     * 推奨パラメータ (NORMALIZE / LOG2 / N_BINS / MAX_PARENTS / TOP_VAR_GENES /
#       ITERS / ITERS_BS / BOOTSTRAP / SEEDS) と、その導出根拠
#     * --emit-config で解析ディレクトリ用の config.sh
#
#   使い方:
#     python3 inspect_matrix.py --input counts.tsv --meta sample_meta.tsv
#     python3 inspect_matrix.py --input expr.xlsx --sheet Sheet1 --emit-config my/config.sh
# =============================================================================
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd


# --- 読み込み ----------------------------------------------------------------

def detect_format(path, fmt):
    if fmt != "auto":
        return fmt
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        return "excel"
    if ext == ".csv":
        return "csv"
    if ext in (".tsv", ".txt", ".tab"):
        return "tsv"
    with open(path, "r", errors="replace") as fh:
        head = fh.readline()
    return "csv" if head.count(",") > head.count("\t") else "tsv"


def count_lines(path):
    n = 0
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(1 << 22)
            if not buf:
                break
            n += buf.count(b"\n")
    return n


def load(path, fmt, sheet, header_row, profile_rows):
    """先頭 profile_rows 行だけ読む (大きなファイルでも軽い)。総行数は別に数える。"""
    if fmt == "excel":
        df = pd.read_excel(path, sheet_name=sheet or 0, header=header_row)
        return df, len(df), False
    sep = "," if fmt == "csv" else "\t"
    df = pd.read_csv(path, sep=sep, header=header_row, nrows=profile_rows,
                     low_memory=False)
    total = count_lines(path) - (header_row + 1)
    return df, max(total, len(df)), total > len(df)


# --- 構造の判定 --------------------------------------------------------------

GENE_ID_HINTS = ("ens", "gene", "symbol", "id", "name", "probe", "feature", "protein")

# 数値だがサンプルではない列 (遺伝子長・座標など) の名前に含まれる語
NUMERIC_ANNOT_HINTS = ("length", "len", "width", "chr", "chrom", "start", "end",
                       "strand", "biotype", "exon", "gc_content", "tss", "entrez",
                       "locus", "position")
LENGTH_HINTS = ("length", "len", "width")


def looks_like_numeric_annotation(col):
    c = str(col).lower()
    return any(h in c for h in NUMERIC_ANNOT_HINTS)


def numeric_mask(df):
    """列ごとに「数値として読めるか」を判定する。"""
    out = []
    for c in df.columns:
        s = df[c]
        if pd.api.types.is_numeric_dtype(s):
            out.append(True)
            continue
        conv = pd.to_numeric(s, errors="coerce")
        out.append(conv.notna().mean() > 0.9)
    return np.array(out)


def guess_orientation(df, nmask, n_rows_total):
    n_num = int(nmask.sum())
    n_ann = int((~nmask).sum())
    evidence = []
    score_gir = 0  # genes-in-rows
    if n_ann >= 1:
        score_gir += 1
        evidence.append(f"数値でない列が {n_ann} 本ある (遺伝子 ID / 注釈列とみられる)")
    if n_rows_total > n_num:
        score_gir += 1
        evidence.append(f"行数 {n_rows_total} > 数値列数 {n_num} (遺伝子の方が多いのが普通)")
    first = str(df.columns[0]).lower()
    if any(h in first for h in GENE_ID_HINTS):
        score_gir += 1
        evidence.append(f"先頭列名 '{df.columns[0]}' が遺伝子 ID 列に見える")
    orient = "genes-in-rows" if score_gir >= 2 else "samples-in-rows"
    return orient, evidence


def classify_values(mat):
    finite = mat[np.isfinite(mat)]
    if finite.size == 0:
        return "unknown", {}
    vmin, vmax = float(finite.min()), float(finite.max())
    integer = bool(np.allclose(finite, np.rint(finite)))
    uniq = np.unique(finite[: min(finite.size, 2_000_000)])
    n_uniq = int(uniq.size)
    zero_frac = float((finite == 0).mean())
    stats = dict(min=vmin, max=vmax, integer=integer, n_unique=n_uniq,
                 zero_frac=zero_frac)
    if integer and n_uniq <= 12 and vmax <= 11 and vmin >= 0:
        kind = "discretized"
    elif integer and vmax > 50:
        kind = "raw_counts"
    elif vmax <= 30 and vmin >= 0:
        kind = "log_transformed"
    elif vmax > 100:
        kind = "linear_normalized"
    else:
        kind = "continuous_unknown"
    return kind, stats


def variance_profile(mat, orientation, kind):
    """変数ごとの分散を (必要なら log スケールで) 求め、四分位を返す。"""
    x = mat.astype(float, copy=True)
    if kind in ("raw_counts", "linear_normalized"):
        with np.errstate(invalid="ignore"):
            x = np.log2(np.clip(x, 0, None) + 1.0)
    axis = 0 if orientation == "samples-in-rows" else 1
    with np.errstate(invalid="ignore"):
        v = np.nanvar(x, axis=axis)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return None
    q = np.percentile(v, [0, 25, 50, 75, 100])
    return dict(q0=float(q[0]), q25=float(q[1]), med=float(q[2]),
                q75=float(q[3]), q100=float(q[4]),
                n=int(v.size),
                # 中位数の半分を下回る「低分散のかたまり」を数える
                n_low=int(np.count_nonzero(v < 0.5 * q[2])),
                n_keep=int(np.count_nonzero(v >= 0.5 * q[2])))


# --- 推奨値の導出 ------------------------------------------------------------

def p_eff(n_samples, r, max_parents=3):
    if n_samples < 10 * r:
        return 0
    return min(max_parents, int(math.floor(math.log(n_samples / 10.0, r))))


def choose_bins(n_samples):
    """P_eff >= 2 を確保できる最大の段階数を選ぶ (詳細は parameter-sizing.md)。"""
    if p_eff(n_samples, 3) >= 2:
        return 3, "N が十分あり 3 値でも親を 2 個保持できる"
    if p_eff(n_samples, 2) >= 2:
        return 2, "3 値だと親 1 個しか保持できないため 2 値にする"
    if p_eff(n_samples, 2) >= 1:
        return 2, "サンプルが少ないので 2 値。親は 1 個までしか保持されない"
    return 2, "サンプルが不足 (N < 20)。この規模では構造学習は困難"


def round_up(x, step=100, floor=500):
    return max(floor, int(math.ceil(x / step) * step))


def suggest(n_samples, n_vars_total, kind, stats, meta_groups, cpu, var_prof=None):
    s = {}
    notes = []

    if kind == "raw_counts":
        s["NORMALIZE"], s["LOG2"] = "cpm", 1
        notes.append("生リードカウントとみて CPM 正規化 + log2(x+1)")
    elif kind == "linear_normalized":
        s["NORMALIZE"], s["LOG2"] = "none", 1
        notes.append("既に線形スケールで正規化済み (TPM/FPKM/CPM) とみて log2 のみ")
    elif kind == "log_transformed":
        s["NORMALIZE"], s["LOG2"] = "none", 0
        notes.append("既に log 変換済みとみて再変換しない (二重 log を避ける)")
    elif kind == "discretized":
        s["NORMALIZE"], s["LOG2"] = "none", 0
        notes.append("既に離散化済み。preprocess.sh は使わず expr_disc.tsv として直接使う")
    else:
        s["NORMALIZE"], s["LOG2"] = "none", 1
        notes.append("値の種類を自動判定できなかった。正規化はユーザに確認する")

    bins, why_bins = choose_bins(n_samples)
    if kind == "discretized":
        bins = int(stats.get("max", 1)) + 1
        why_bins = "データに含まれる水準数から決定"
    s["N_BINS"] = bins
    pe = p_eff(n_samples, bins) if bins >= 2 else 0
    s["MAX_PARENTS"] = max(1, min(3, pe)) if pe > 0 else 1
    notes.append(f"N_BINS={bins} ({why_bins}) → P_eff = floor(log_{bins}(N/10)) = {pe}")
    if pe == 0:
        notes.append("警告: P_eff = 0。この N と離散化では罰則付きスコアが親を保持できない。"
                     "サンプルを増やす / N_BINS を下げる / 変数を絞る が必要")

    # 解析対象の変数数
    if kind == "discretized":
        d = n_vars_total
        s["TOP_VAR_GENES"] = 0
    elif n_vars_total > 500:
        d = 500
        s["TOP_VAR_GENES"] = 500
        notes.append(f"変数が {n_vars_total} 本あるので分散上位 500 に絞る "
                     "(全部使うなら TOP_VAR_GENES=0。計算時間は D に比例)")
    elif var_prof and var_prof["n_low"] > 0.15 * var_prof["n"]:
        # 低分散 (ほぼ一定) の変数が多い。残すと偽のエッジを増やすので落とす。
        keep = max(20, var_prof["n_keep"])
        d = min(n_vars_total, keep)
        s["TOP_VAR_GENES"] = d
        frac = var_prof["n_low"] / var_prof["n"]
        notes.append(f"分散が中位数の半分未満の変数が {frac:.0%} ある "
                     f"(低分散のかたまり。残すと偽のエッジを増やす) ので分散上位 {d} に絞る")
    else:
        d = n_vars_total
        s["TOP_VAR_GENES"] = 0
    s["_D"] = d
    s["_FIRST_PASS_TOP_VAR"] = min(d, 150)

    pe_use = max(1, pe)
    s["ITERS"] = round_up(3 * d * pe_use)
    s["ITERS_BS"] = round_up(d * pe_use)
    notes.append(f"ITERS = 3 x D x P_eff = 3 x {d} x {pe_use}、ITERS_BS = D x P_eff "
                 "(切り上げ)")

    seeds = max(1, min(cpu, 10))
    target_b = 200
    s["SEEDS"] = seeds
    s["BOOTSTRAP"] = int(math.ceil(target_b / seeds))
    s["MAX_JOBS"] = seeds
    b = s["SEEDS"] * s["BOOTSTRAP"]
    notes.append(f"総リサンプル数 B = {b} (SE <= {0.5 / math.sqrt(b):.3f})。"
                 "まず BOOTSTRAP=3 SEEDS=2 で 1 回計測してから本番の回数を決める")

    zero_frac = stats.get("zero_frac", 0.0)
    if kind in ("raw_counts", "linear_normalized", "log_transformed"):
        s["MIN_DETECT_FRAC"] = 0.2 if zero_frac > 0.5 else 0.5
        if zero_frac > 0.5:
            notes.append(f"ゼロが {zero_frac:.0%} を占める (単一細胞のような疎なデータ)。"
                         "検出率フィルタは緩める")
    else:
        s["MIN_DETECT_FRAC"] = 0.0

    s["SCORE"] = "bdeu"
    s["ESS"] = 10 if n_samples >= 100 else 1
    s["SCORE_IMP"] = "bic"
    s["THRESHOLD_PROB"] = 0.3
    s["_GROUPS"] = meta_groups
    return s, notes


# --- config.sh の生成 --------------------------------------------------------

def emit_config(path, args, s, kind, abs_expr, abs_meta, orientation, ann):
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    disc = kind == "discretized"
    L = []
    a = L.append
    a("#!/usr/bin/env bash")
    a("# =============================================================================")
    a("# config.sh — inspect_matrix.py が生成 (推奨値。実データに合わせて必ず確認する)")
    a("#   使い方: source ./config.sh してから ${BN_SCRIPTS}/*.sh を呼ぶ")
    a("#           変数の一覧は ${FASTBN_HOME}/script/README.md")
    a("# =============================================================================")
    a("")
    a("# --- 出力先 (RUNDIR を変えると設定違いの実験を並べて比較できる) --------------")
    a('export RUNDIR="${RUNDIR:-.}"')
    a('export DATADIR="${DATADIR:-${RUNDIR}/data}"')
    a('export OUTDIR="${OUTDIR:-${RUNDIR}/out}"')
    a('export BSDIR="${BSDIR:-${RUNDIR}/bs}"')
    a('export GROUPDIR="${GROUPDIR:-${RUNDIR}/groups}"')
    a('export FIGDIR="${FIGDIR:-${RUNDIR}/figures}"')
    a('export FIGDIR_BS="${FIGDIR_BS:-${RUNDIR}/figures_bs}"')
    a('export REPORT_HTML="${REPORT_HTML:-${RUNDIR}/report.html}"')
    a('export TARGET_FILE="${TARGET_FILE:-${RUNDIR}/target_genes.txt}"  # 注目遺伝子 (1 行 1 名前)')
    a("")
    a("# --- fast_bn 入力 (前処理の出力) ---------------------------------------------")
    a('export INPUT="${INPUT:-${DATADIR}/expr_disc.tsv}"')
    a('export VARMAP="${VARMAP:-${DATADIR}/var_map.tsv}"')
    a('export SAMPLES="${SAMPLES:-${DATADIR}/samples.tsv}"')
    a("")
    if disc:
        a("# --- 入力データ ---------------------------------------------------------------")
        a("# このデータは既に離散化済み。preprocess.sh は使わない。")
        a(f"#   元ファイル: {abs_expr}")
        a("# ${DATADIR}/expr_disc.tsv (行=サンプル, 列=変数, 0 始まりの整数コード) と")
        a("# var_map.tsv / samples.tsv を用意すること")
        a("# (実例: ${FASTBN_HOME}/example_sc/prepare_data.py, 01prepare.sh)。")
        a(f'export SRC_MATRIX="{abs_expr}"')
        if abs_meta:
            a(f'export SAMPLE_META="{abs_meta}"')
    else:
        a("# --- 入力データ ---------------------------------------------------------------")
        a(f'export EXPR_INPUT="${{EXPR_INPUT:-{abs_expr}}}"')
        a(f'export ORIENTATION="${{ORIENTATION:-{orientation}}}"')
        for key, col in ann.items():
            if key == "_LENGTH_COL":
                a(f'#export LENGTH_COL="{col}"   # NORMALIZE=tpm にするなら有効化し DROP_COLS から外す')
            elif col is None:
                a(f'#export {key}="..."     # 該当する列があれば指定')
            else:
                a(f'export {key}="${{{key}:-{col}}}"')
        if abs_meta:
            a(f'export SAMPLE_META="${{SAMPLE_META:-{abs_meta}}}"   # sample_id + group の表')
        else:
            a('#export SAMPLE_META="/path/to/sample_meta.tsv"   # 群別解析にはこれが必要')
        a("")
        a("# --- 前処理 -------------------------------------------------------------------")
        a(f'export NORMALIZE="${{NORMALIZE:-{s["NORMALIZE"]}}}"        # none | cpm | tpm')
        a(f'export LOG2="${{LOG2:-{s["LOG2"]}}}"                  # 1 で log2(x+1)')
        a(f'export MIN_DETECT_FRAC="${{MIN_DETECT_FRAC:-{s["MIN_DETECT_FRAC"]}}}"')
        a(f'export TOP_VAR_GENES="${{TOP_VAR_GENES:-{s["TOP_VAR_GENES"]}}}"    # 0 = 全変数')
        a(f'export N_BINS="${{N_BINS:-{s["N_BINS"]}}}"')
        a('export DISC_METHOD="${DISC_METHOD:-quantile}"   # quantile | uniform')
    a("")
    a("# --- 構造学習 -----------------------------------------------------------------")
    a(f'export SCORE="${{SCORE:-{s["SCORE"]}}}"                # bic | k2 | bdeu')
    a(f'export ESS="${{ESS:-{s["ESS"]}}}"                     # BDeu の等価サンプルサイズ')
    a(f'export MAX_PARENTS="${{MAX_PARENTS:-{s["MAX_PARENTS"]}}}"')
    a(f'export ITERS="${{ITERS:-{s["ITERS"]}}}"')
    a('export TABU="${TABU:-30}"')
    a('export TOPK="${TOPK:-20}"                    # 候補親の上位 K')
    a('export REACH="${REACH:-lazy}"                # lazy | dense')
    a("")
    a("# --- エッジ重要度 -------------------------------------------------------------")
    a(f'export SCORE_IMP="${{SCORE_IMP:-{s["SCORE_IMP"]}}}"')
    a('export ALPHA="${ALPHA:-1.0}"')
    a("")
    a("# --- ブートストラップ (総リサンプル数 = BOOTSTRAP x SEEDS。最も重い) ----------")
    a(f'export BOOTSTRAP="${{BOOTSTRAP:-{s["BOOTSTRAP"]}}}"')
    a(f'export SEEDS="${{SEEDS:-{s["SEEDS"]}}}"')
    a(f'export MAX_JOBS="${{MAX_JOBS:-{s["MAX_JOBS"]}}}"           # 同時実行数 (CPU コア数)')
    a(f'export ITERS_BS="${{ITERS_BS:-{s["ITERS_BS"]}}}"')
    a(f'export THRESHOLD_PROB="${{THRESHOLD_PROB:-{s["THRESHOLD_PROB"]}}}"')
    a('export THRESHOLD_COUNT="${THRESHOLD_COUNT:-2}"')
    a("")
    a("# --- 可視化・レポート ---------------------------------------------------------")
    a('# viz*.sh はこの 2 つを読まないので --metrics / --top-n として明示的に渡すこと')
    a('export VIZ_METRICS="${VIZ_METRICS:-dlogL,dBIC}"   # dK2,dBDeu も指定可')
    a('export VIZ_TOP_N="${VIZ_TOP_N:-40}"')
    a('export REPORT_TITLE="${REPORT_TITLE:-ベイジアンネットワーク解析レポート}"')
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    os.chmod(path, 0o644)
    return path


# --- 本体 --------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="発現量行列を点検して FastBN の推奨設定を出す")
    ap.add_argument("--input", required=True, help="発現量ファイル (TSV/CSV/XLSX)")
    ap.add_argument("--format", default="auto", choices=["auto", "tsv", "csv", "excel"])
    ap.add_argument("--sheet", default=None, help="Excel のシート名")
    ap.add_argument("--header-row", type=int, default=0, help="ヘッダ行の位置 (0 始まり)")
    ap.add_argument("--orientation", default="auto",
                    choices=["auto", "genes-in-rows", "samples-in-rows"])
    ap.add_argument("--meta", default=None, help="サンプル ID + 群ラベルの表")
    ap.add_argument("--meta-sample-col", default=None)
    ap.add_argument("--meta-group-col", default="group")
    ap.add_argument("--profile-rows", type=int, default=5000,
                    help="値の profiling に読む行数 (既定 5000)")
    ap.add_argument("--emit-config", default=None, help="config.sh の出力先")
    args = ap.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"エラー: ファイルがありません: {args.input}")

    fmt = detect_format(args.input, args.format)
    df, n_rows_total, truncated = load(args.input, fmt, args.sheet,
                                       args.header_row, args.profile_rows)
    nmask = numeric_mask(df)
    ann_cols = [c for c, m in zip(df.columns, nmask) if not m]
    num_cols = [c for c, m in zip(df.columns, nmask) if m]

    if args.orientation == "auto":
        orientation, evidence = guess_orientation(df, nmask, n_rows_total)
    else:
        orientation, evidence = args.orientation, ["--orientation で明示指定"]

    # 数値だがサンプルではない列 (遺伝子長・座標など) を分離する
    num_annot_cols = []
    if orientation == "genes-in-rows":
        num_annot_cols = [c for c in num_cols if looks_like_numeric_annotation(c)]
        num_cols = [c for c in num_cols if c not in num_annot_cols]

    num = df[num_cols].apply(pd.to_numeric, errors="coerce")
    mat = num.to_numpy(dtype=float)          # 読んだ範囲 (行=ファイルの行)
    if orientation == "genes-in-rows":
        n_samples, n_vars = len(num_cols), n_rows_total
        feature_names = df[ann_cols[0]].astype(str).tolist() if ann_cols else []
        sample_names = [str(c) for c in num_cols]
    else:
        n_samples, n_vars = n_rows_total, len(num_cols)
        feature_names = [str(c) for c in num_cols]
        sample_names = df[ann_cols[0]].astype(str).tolist() if ann_cols else []

    kind, stats = classify_values(mat)
    n_missing = int(np.count_nonzero(~np.isfinite(mat)))
    n_const = int(np.count_nonzero(np.nanstd(mat, axis=0 if orientation == "samples-in-rows" else 1) == 0))
    dup = 0
    if feature_names:
        dup = len(feature_names) - len(set(feature_names))

    # サンプルメタデータ
    meta_groups, meta_info = None, []
    if args.meta:
        if not os.path.isfile(args.meta):
            sys.exit(f"エラー: --meta がありません: {args.meta}")
        msep = "," if detect_format(args.meta, "auto") == "csv" else "\t"
        mdf = pd.read_csv(args.meta, sep=msep)
        gcol = args.meta_group_col if args.meta_group_col in mdf.columns else None
        if gcol is None:
            cand = [c for c in mdf.columns if str(c).lower() in
                    ("group", "condition", "treatment", "tissue", "celltype", "cell_type")]
            gcol = cand[0] if cand else None
        scol = args.meta_sample_col or str(mdf.columns[0])
        meta_info.append(f"列: {list(mdf.columns)}")
        meta_info.append(f"サンプル ID 列 = {scol} / 群ラベル列 = {gcol}")
        if gcol:
            vc = mdf[gcol].value_counts()
            meta_groups = vc.to_dict()
            meta_info.append("群ごとのサンプル数: " +
                             ", ".join(f"{k}={v}" for k, v in vc.items()))
            small = [k for k, v in vc.items() if v < 2]
            if small:
                meta_info.append(f"警告: サンプル数 1 の群があり群別解析で落ちる: {small}")
        ids = set(mdf[scol].astype(str))
        overlap = len(ids & set(sample_names)) if sample_names else 0
        meta_info.append(f"行列のサンプル名との一致: {overlap} / {len(sample_names)}")
        if sample_names and overlap == 0:
            meta_info.append("警告: サンプル名が 1 つも一致しない。列の向きか ID 列を確認する")

    var_prof = variance_profile(mat, orientation, kind)
    cpu = os.cpu_count() or 4
    s, notes = suggest(n_samples, n_vars, kind, stats, meta_groups, cpu, var_prof)

    KIND_JA = {
        "raw_counts": "生リードカウント (整数, 大きな値)",
        "linear_normalized": "線形スケールの正規化済み (TPM/FPKM/CPM など)",
        "log_transformed": "log 変換済み",
        "discretized": "既に離散化済み (小さな整数コード)",
        "continuous_unknown": "連続値 (種類を特定できず)",
        "unknown": "判定不能",
    }

    out = print
    out("=" * 70)
    out(f" 入力: {args.input}  (形式 {fmt})")
    out("=" * 70)
    out(f" 行 x 列            : {n_rows_total} x {len(df.columns)}"
        + ("  (値の profiling は先頭 %d 行のみ)" % len(df) if truncated else ""))
    out(f" 向きの推定          : {orientation}")
    for e in evidence:
        out(f"   - {e}")
    out(f" 注釈列 (非数値)     : {ann_cols if ann_cols else 'なし'}")
    if num_annot_cols:
        out(f" 注釈列 (数値だが列名から注釈と判断): {num_annot_cols}")
        out("   -> サンプル列から除外した。誤りなら --orientation / 列指定で修正する")
    out(f" => サンプル数 N     : {n_samples}")
    out(f" => 変数数 D         : {n_vars}")
    out(f" 値の種類            : {KIND_JA.get(kind, kind)}")
    out(f"   範囲 {stats.get('min')} 〜 {stats.get('max')} / 整数 {stats.get('integer')} "
        f"/ 異なる値 {stats.get('n_unique')} / ゼロ率 {stats.get('zero_frac', 0):.1%}")
    if n_missing:
        out(f" 警告: 欠損・非数値が {n_missing} セル (前処理前に埋めるか除く)")
    if n_const:
        out(f" 注意: 定数 (分散 0) の変数が {n_const} 本 (離散化で水準 1 になる)")
    if var_prof:
        out(f" 変数ごとの分散      : 最小 {var_prof['q0']:.3g} / 25% {var_prof['q25']:.3g}"
            f" / 中位 {var_prof['med']:.3g} / 75% {var_prof['q75']:.3g}"
            f" / 最大 {var_prof['q100']:.3g}"
            + ("  (log2 スケール)" if kind in ("raw_counts", "linear_normalized") else ""))
        if var_prof["n_low"]:
            out(f"   うち中位数の半分未満: {var_prof['n_low']} 本"
                f" / {var_prof['n']} 本 ({var_prof['n_low'] / var_prof['n']:.0%})")
    if dup:
        out(f" 注意: 重複する変数名が {dup} 件")
    if meta_info:
        out("-" * 70)
        out(f" サンプルメタデータ: {args.meta}")
        for m in meta_info:
            out(f"   {m}")
    out("-" * 70)
    out(" 推奨設定 (根拠は下の注記と references/parameter-sizing.md)")
    keys = ["NORMALIZE", "LOG2", "MIN_DETECT_FRAC", "TOP_VAR_GENES", "N_BINS",
            "SCORE", "ESS", "MAX_PARENTS", "ITERS",
            "BOOTSTRAP", "SEEDS", "MAX_JOBS", "ITERS_BS", "THRESHOLD_PROB", "SCORE_IMP"]
    for k in keys:
        if k in s:
            out(f"   {k:<16} = {s[k]}")
    out("-" * 70)
    for n in notes:
        out(f" * {n}")
    out(f" * 初回は TOP_VAR_GENES={s['_FIRST_PASS_TOP_VAR']} / BOOTSTRAP=3 SEEDS=2 で"
        " 1 周させて所要時間を測る")
    if kind == "discretized":
        out(" * 既に離散化済みなので preprocess.sh は使わない "
            "(expr_disc.tsv / var_map.tsv / samples.tsv を直接用意する)")
    if not args.meta:
        out(" * 群 (条件) 別の比較をするなら sample_id + group の表を用意して "
            "SAMPLE_META に指定する")
    out("=" * 70)

    if args.emit_config:
        ann = {}
        if orientation == "genes-in-rows":
            if ann_cols:
                ann["ID_COL"] = ann_cols[0]
                ann["NAME_COL"] = ann_cols[1] if len(ann_cols) > 1 else None
            drop = [str(c) for c in ann_cols[2:]] + [str(c) for c in num_annot_cols]
            length_col = next((c for c in num_annot_cols
                               if any(h in str(c).lower() for h in LENGTH_HINTS)), None)
            if drop:
                ann["DROP_COLS"] = ",".join(drop)
            if length_col:
                ann["_LENGTH_COL"] = str(length_col)
        elif ann_cols:
            ann["SAMPLE_ID_COL"] = ann_cols[0]
        p = emit_config(args.emit_config, args, s, kind,
                        os.path.abspath(args.input),
                        os.path.abspath(args.meta) if args.meta else None,
                        orientation, ann)
        print(f"config.sh を書きました: {p}")
        print("  中身を必ず確認すること (特に ID_COL / NAME_COL / DROP_COLS / NORMALIZE)")


if __name__ == "__main__":
    main()

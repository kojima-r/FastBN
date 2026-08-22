#!/usr/bin/env bash
# =============================================================================
# config.sh — example_sc の全設定
#   各ステップスクリプト (00〜08) は先頭でこのファイルを source する。
#   ここで export した変数が ../script/*.sh の既定値を上書きする。
#   変数の一覧は ../script/README.md を参照。
#
#   すべての設定は環境変数で上書きできる (このファイルを編集しなくてよい)。
#     DATASET=ss ./run_all.sh
#     DISC=tri NVARS=1000 ./02learn.sh
# =============================================================================

# --- どのデータを使うか (この 3 つが例題の切り替えスイッチ) -----------------
# DATASET : bbknn = 複数プロトコル統合 (BBKNN 補正後; 2488 遺伝子 x 240 サンプル)
#           ss    = Smart-seq2 のみ            (58246 遺伝子 x 450 サンプル)
export DATASET="${DATASET:-bbknn}"

# DISC    : bin = 2 値離散化 / tri = 3 値離散化
export DISC="${DISC:-bin}"

# NVARS   : 使用する遺伝子数 (先頭 N 列)。10 | 100 | 1000 | all
#           既定は all (そのデータセットの全遺伝子)。
#           all は変数が数千〜数万になり計算時間が大きく伸びる:
#             bbknn (2488 遺伝子)  -> run_all.sh 全体で 2 時間程度
#             ss    (58246 遺伝子) -> さらに桁違いに長い
#           まず流れだけ確認したいときは NVARS=100 (数分) を指定すること。
export NVARS="${NVARS:-all}"

# --- 上の 3 つから導出されるパス (通常は編集不要) ---------------------------
case "${DATASET}" in
  bbknn) export DATA_ROOT="./data_bbknn_r_tissue_disc" ;;
  ss)    export DATA_ROOT="./data_ss_r_tissue_disc" ;;
  *) echo "[config] エラー: DATASET は bbknn | ss のいずれかです (指定: ${DATASET})" >&2
     return 1 2>/dev/null || exit 1 ;;
esac

case "${DISC}" in
  bin) _disc_sfx="";     _tissue_sub="tissue" ;;
  tri) _disc_sfx="_tri"; _tissue_sub="tissue_tri" ;;
  *) echo "[config] エラー: DISC は bin | tri のいずれかです (指定: ${DISC})" >&2
     return 1 2>/dev/null || exit 1 ;;
esac

case "${NVARS}" in
  all) _nvars_sfx="" ;;
  ''|*[!0-9]*) echo "[config] エラー: NVARS は 10 | 100 | 1000 | all です (指定: ${NVARS})" >&2
     return 1 2>/dev/null || exit 1 ;;
  *) _nvars_sfx="${NVARS}" ;;
esac

# 解析対象の行列と、群 (組織) ラベルの元になる組織別ファイル
export SRC_MATRIX="${DATA_ROOT}/all_disc${_disc_sfx}${_nvars_sfx}.tsv"
export SRC_TISSUE_DIR="${DATA_ROOT}/${_tissue_sub}"

# --- 出力先 -------------------------------------------------------------------
# 設定を切り替えても結果が混ざらないよう、組み合わせごとに RUNDIR を分ける。
export RUNDIR="${RUNDIR:-./run_${DATASET}_${DISC}_${NVARS}}"
export DATADIR="${RUNDIR}/data"
export OUTDIR="${RUNDIR}/out"
export BSDIR="${RUNDIR}/bs"
export GROUPDIR="${RUNDIR}/groups"
export FIGDIR="${RUNDIR}/figures"
export FIGDIR_BS="${RUNDIR}/figures_bs"
export REPORT_HTML="${RUNDIR}/report.html"
export INPUT="${DATADIR}/expr_disc.tsv"
export VARMAP="${DATADIR}/var_map.tsv"
export SAMPLES="${DATADIR}/samples.tsv"
export TARGET_FILE="${RUNDIR}/target_genes.txt"

# --- 注目遺伝子 (任意) --------------------------------------------------------
# 図で強調したい遺伝子をカンマ区切りで指定する。データに無い名前は無視される。
# bbknn は遺伝子シンボル (Sox17 など)、ss は Ensembl ID (ENSG... ) が列名。
export TARGET_GENES="${TARGET_GENES:-}"

# --- 構造学習 (02learn.sh) ----------------------------------------------------
export SCORE="${SCORE:-bdeu}"          # bic | k2 | bdeu
export ESS="${ESS:-10}"                # BDeu の等価サンプルサイズ
export MAX_PARENTS="${MAX_PARENTS:-3}" # 各ノードの最大親数
export ITERS="${ITERS:-5000}"          # Hill-Climb の最大反復数
export TABU="${TABU:-30}"              # Tabu サーチの禁制期間
export TOPK="${TOPK:-20}"              # 候補親の上位 K 制限
export JINDEX_CACHE="${JINDEX_CACHE:-1024}"
export REACH="${REACH:-lazy}"          # dense | lazy
# ログは既定で詳細に出る (VERBOSE=1 で --verbose を明示指定できる)。

# --- エッジ重要度 (03importance.sh / 05importance_groups.sh) ------------------
export SCORE_IMP="${SCORE_IMP:-bic}"   # 重要度評価に使うスコア
export ALPHA="${ALPHA:-1.0}"           # スムージング係数

# --- ブートストラップ (04bootstrap.sh) ---------------------------------------
export BOOTSTRAP="${BOOTSTRAP:-10}"    # 1 シードあたりのリサンプリング回数
export SEEDS="${SEEDS:-5}"             # シード数 (総リサンプル数 = 積 = 50)
export MAX_JOBS="${MAX_JOBS:-5}"       # 同時実行プロセス数 (CPU コア数に合わせる)
export ITERS_BS="${ITERS_BS:-2000}"    # リサンプルごとの反復数 (学習より少なめ)
export THRESHOLD_PROB="${THRESHOLD_PROB:-0.3}"   # コンセンサス採用の確率閾値
export THRESHOLD_COUNT="${THRESHOLD_COUNT:-2}"   # コンセンサス採用の回数閾値
export WARM_START="${WARM_START:-1}"   # 1 で out/edges.tsv を初期構造に使う

# --- 群 (組織) 別解析 (05importance_groups.sh) -------------------------------
# 群ラベルは 01prepare.sh が作る data/samples.tsv の group 列 (= 組織名)。
export MIN_GROUP_SAMPLES="${MIN_GROUP_SAMPLES:-2}"

# --- 可視化・レポート (06visualize.sh / 07report.sh) -------------------------
export VIZ_METRICS="${VIZ_METRICS:-dlogL,dBIC}"  # 全て: dlogL,dBIC,dK2,dBDeu
export VIZ_TOP_N="${VIZ_TOP_N:-40}"              # 強調するエッジ本数
export REPORT_TITLE="${REPORT_TITLE:-単一細胞発現データ BN 解析レポート (${DATASET}/${DISC}/${NVARS})}"

# --- データ取得 (00download.sh) ----------------------------------------------
# アーカイブは <DATA_BASE_URL>/<ディレクトリ名>.tar.gz として取得する。
export DATA_BASE_URL="${DATA_BASE_URL:-https://github.com/kojima-r/FastBN/releases/download/v0.1}"

# 汎用スクリプトの場所 (このディレクトリからの相対)
export BN_SCRIPTS="${BN_SCRIPTS:-../script}"

unset _disc_sfx _tissue_sub _nvars_sfx

#!/usr/bin/env bash
# =============================================================================
# config.sh — example_bulk の全設定
#   各ステップスクリプト (00〜08) は先頭でこのファイルを source する。
#   自分のデータで解析する場合は、このファイルをコピーして書き換えるだけで
#   ../script/ 以下の汎用スクリプトをそのまま使える。
#
#   ここで export した変数が ../script/*.sh の既定値を上書きする。
#   一覧は ../script/README.md を参照。
# =============================================================================

# --- 入力データ --------------------------------------------------------------
# 発現量ファイル (このサンプルでは 00make_data.sh が生成するダミーデータ)
export EXPR_INPUT="./data/counts.tsv"
export ID_COL="gene_id"          # 遺伝子 ID の列
export NAME_COL="gene_name"      # 遺伝子シンボルの列 (図のラベルに使う)
export DROP_COLS="gene_length"   # 解析に使わない注釈列
#export LENGTH_COL="gene_length" # NORMALIZE=tpm にする場合はこちらを有効化
export SAMPLE_META="./data/sample_meta.tsv"   # サンプル ID と群ラベルの表
export TARGET_FILE="./target_genes.txt"       # 注目遺伝子 (必ず残す)

# --- 前処理 ------------------------------------------------------------------
export NORMALIZE="cpm"           # 生カウントなので CPM 正規化 (TPM 済みなら none)
export LOG2=1                    # log2(x + 1) 変換
export MIN_DETECT_FRAC=0.5       # 半数以上のサンプルで検出される遺伝子のみ
export TOP_VAR_GENES=60          # 分散上位 60 遺伝子を解析対象に
export N_BINS=3                  # 3 段階 (低/中/高) に離散化
export DISC_METHOD="quantile"    # 等頻度 (分位点) 離散化

# --- 構造学習 ----------------------------------------------------------------
export SCORE="bdeu"              # bic | k2 | bdeu
export ESS=10                    # BDeu の等価サンプルサイズ
export MAX_PARENTS=2             # 各ノードの最大親数 (サンプルが少ないほど小さく)
export ITERS=5000                # Hill-Climb の最大反復数
export TABU=30                   # Tabu サーチの禁制期間
export TOPK=20                   # 候補親の上位 K 制限

# --- エッジ重要度 ------------------------------------------------------------
export SCORE_IMP="bic"           # 重要度評価に使うスコア
export ALPHA=1.0                 # スムージング係数

# --- ブートストラップ (エッジ安定性) ----------------------------------------
export BOOTSTRAP=10              # 1 シードあたりのリサンプリング回数
export SEEDS=5                   # シード数 (総リサンプル数 = BOOTSTRAP x SEEDS)
export MAX_JOBS=5                # 同時実行プロセス数 (CPU コア数に合わせる)
export ITERS_BS=1500             # リサンプルごとの反復数 (学習より少なめでよい)
export THRESHOLD_PROB=0.3        # コンセンサス採用のブートストラップ確率閾値
export THRESHOLD_COUNT=2         # コンセンサス採用の出現回数閾値

# --- 可視化・レポート --------------------------------------------------------
export VIZ_METRICS="dlogL,dBIC"  # 図を作るメトリクス (全て: dlogL,dBIC,dK2,dBDeu)
export VIZ_TOP_N=40              # 強調するエッジ本数
export REPORT_TITLE="バルク RNA 発現データ BN 解析レポート (ダミーデータ)"

# --- ダミーデータ生成 (00make_data.sh 用) ------------------------------------
export DUMMY_GENES=60            # 真の DAG に参加する遺伝子数
export DUMMY_NOISE=30            # 無相関・低分散のノイズ遺伝子数 (フィルタで落ちる)
export DUMMY_GROUPS="Control,TreatA,TreatB,Combo"
export DUMMY_REPLICATES=30       # 各群のサンプル数 (合計 120 サンプル)
export DUMMY_MAX_PARENTS=2       # 真の DAG の最大親数
export DUMMY_EDGE_PROB=0.06      # 真の DAG の密度 (概ね 1 ノードあたり 1 エッジ)
export DUMMY_SIGNAL_FRAC=0.65    # 変動のうち親で説明される割合 (親子相関 ≈ 0.8)
export DUMMY_SEED=7              # 乱数シード (再現性)

# --- 出力先 (既定のままで良い場合は変更不要) --------------------------------
export DATADIR="./data"
export OUTDIR="./out"
export BSDIR="./bs"
export GROUPDIR="./groups"
export INPUT="${DATADIR}/expr_disc.tsv"
export VARMAP="${DATADIR}/var_map.tsv"
export SAMPLES="${DATADIR}/samples.tsv"

# 汎用スクリプトの場所 (このディレクトリからの相対)
export BN_SCRIPTS="../script"

#include <string>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <exception>

#include "fast_bn_lib.hpp"
#include "fast_bn_dataset.hpp"
#include "fast_bn_dag.hpp"
#include "fast_bn_score.hpp"

#ifdef __NVCOMPILER
#include <openacc.h>
#endif

static void print_help() {
    const char* env_lang = std::getenv("LANG");
    std::string lang = (env_lang ? std::string(env_lang) : "ja");

    if (lang.substr(0,2) == "en") {
        std::cout <<
R"(Bayesian Network Structure Learning Tool (Hill-Climb + Tabu + Bootstrap)
---------------------------------------------------------------------
Usage:
  fast_bn --input <data.tsv> --score <bic|k2|bdeu> [options]

General options:
  -h, --help
      Show this help message and exit.
  --input <file>
      Input dataset file (CSV/TSV with header, separator auto-detected).
  --score <bic|k2|bdeu>
      Scoring function (default=bic).
  --init <file>
      Initial network structure (tab-separated u v).
  --iters <n>
      Max iterations for Hill-Climb (default=1000).
  --tabu <len>
      Tabu search tenure. 0 disables tabu (default=0).
  --max-parents <k>
      Maximum number of parents per node (default=5).
  --max-children <k>
      Maximum number of children per node (default=unlimited).
  --ess <val>
      Equivalent sample size for BDeu (default=1.0).
  --alpha <a>
      Smoothing coefficient for likelihood (default=1.0).
  --verbose
      Enable verbose output.

Candidate parent selection:
  --cand-metric mi|chi2
      Measure for candidate parent selection: mi=Mutual Information, chi2=Chi-square (default=mi).
  --topk <K>
      Keep top-K strongest candidates (0=no limit, default=50).
  --mi-threshold <t>
      Keep variables with MI ≥ t (in nats, default=0.0).
  --chi2-p-threshold <p>
      Keep variables with p ≤ threshold (default=0.05).
  --mi-sample <n>
      Number of samples for MI/Chi2 calculation (0=use all, default=0).
  --mi-budget <n>
      Candidate variable limit per node (0=all, default=0).
  --reach <dense|lazy>
      Reachability check mode (default=lazy).
  --jindex-cache <cap>
      Parent configuration cache size (default=32).

Bootstrap mode:
  --bootstrap <B>
      Run bootstrap sampling B times (default=0).
  --seed <s>
      Random seed (default=2025).
  --save-bootstrap-counts <path>
      Output base filename; "_seed####.tsv" will be appended automatically.
  --bootstrap-include-zero
      Include edges that never appeared (default=off).

Edge importance mode:
  --edge-importance
      Evaluate score changes by removing each edge.
  --score-dataset <file>
      New dataset for evaluation.
  --init <file>
      Trained network structure (tab-separated u v).
  --counts <file>
      Trained all_counts.tsv file.
  --save-edge-importance <path>
      Output file for edge importance (TSV).

Output & runtime:
  --verbose
      Enable detailed logging.
---------------------------------------------------------------------
Example:
./fast_bn --input data.tsv --score bic --iters 2000 --tabu 20 --max-parent 3
./fast_bn --input gene.tsv --score bdeu --ess 10 --iters 3000

Example: bootstrap
./fast_bn --input data.tsv --score bic \
  --bootstrap 100 \
  --save-bootstrap-counts results/boot_edges.tsv

 
Example: edge importance
./fast_bn --input test.tsv --score bic \
  --edge-importance \
  --init init_edges.tsv \
  --counts all_counts.tsv \
  --save-edge-importance edge_imp.tsv
---------------------------------------------------------------------
)";
    }else{
        std::cout <<
R"(Bayesian Network Structure Learning Tool (Hill-Climb + Tabu + Bootstrap)
---------------------------------------------------------------------
Usage:
  fast_bn --input <data.tsv> --score <score> [options]


General options:
  -h, --help
      このヘルプを表示して終了します。
  --input <file>
      入力データファイル（CSV/TSV, ヘッダ付き, 区切り自動判定）。
  --score <bic|k2|bdeu>
      スコア関数を指定します（default=bic）。
  --init <file>
      初期構造ファイル (タブ区切り u v)。
  --iters <n>
      Hill-Climb の最大反復数（default=1000）。
  --tabu <len>
      Tabu サーチの禁制期間。0で無効（default=0）。
  --max-parents <k>
      各ノードの最大親数（default=5）。
  --max-children <k>
      各ノードの最大子数（default=無制限）。
  --ess <val>
      BDeu スコアの等価サンプルサイズ（default=1.0）。
  --alpha <a>
      スムージング係数（default=1.0）。
  --verbose
      詳細ログを出力。

Candidate parent selection:
  --cand-metric mi|chi2
      候補親の関連度指標。mi=相互情報量, chi2=カイ二乗p値（default=mi）。
  --topk <K>
      候補親を上位K個に制限。0で制限なし（default=50）。
  --mi-threshold <t>
      MI >= t (nats単位) のみ採用（default=0.0）。
  --chi2-p-threshold <p>
      p <= p_threshold のみ採用（default=0.05）。
  --mi-sample <n>
      MI/chi2 計算時のサンプル数（0=全行, default=0）。
  --mi-budget <n>
      MI/chi2 計算時の候補数上限（0=全変数, default=0）。
  --reach <dense|lazy>
      到達可能性チェックモード（default=lazy）。
  --jindex-cache <cap>
      親配置キャッシュサイズ（default=32）。

Bootstrap mode:
  --bootstrap <B>
      ブートストラップ試行回数。B>0で実行モード（default=0）。
  --seed <s>
      乱数シード値（default=2025）。
  --save-bootstrap-counts <path>
      出力ファイルベース名。自動で "_seed####.tsv" が付与されます。
  --bootstrap-include-zero
      出現しなかったエッジも出力（default=off）。

Edge importance mode:
  --edge-importance
      エッジ除去によるスコア変化を評価。
  --score-dataset <file>
      評価用データセット（新しいCSV/TSV）。
  --init <file>
      学習済み構造（タブ区切り u v）。
  --counts <file>
      学習済み all_counts.tsv。
  --save-edge-importance <path>
      エッジ重要度を出力するTSVファイル。

Output & runtime:
  --save-bootstrap-counts <path>
      ブートストラップの出力ファイル（ベース名）。
  --seed <n>
      ブートストラップの乱数シード（default=2025）。
  --bootstrap-include-zero
      出現しないエッジも出力。
  --verbose
      詳細ログ出力。

---------------------------------------------------------------------
例:
./fast_bn --input data.tsv --score bic --iters 2000 --tabu 20 --max-parent 3
./fast_bn --input gene.tsv --score bdeu --ess 10 --iters 3000

例: bootstrap
./fast_bn --input data.tsv --score bic \
  --bootstrap 100 \
  --save-bootstrap-counts results/boot_edges.tsv

 
例: エッジ貢献
./fast_bn --input test.tsv --score bic \
  --edge-importance \
  --init init_edges.tsv \
  --counts all_counts.tsv \
  --save-edge-importance edge_imp.tsv
---------------------------------------------------------------------
)";
    }
}

int main(int argc, char** argv){
#ifdef __NVCOMPILER
    acc_init(acc_device_nvidia);
    std::cout << "Device Type: " << acc_get_device_type() << std::endl;
#endif
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc == 1){
        std::cerr << "Usage:\n"
             << "  "<<argv[0]<<" --input <csv_path> --score <score(bic|k2|bdeu)>\n"
             << "    [--ess E] [--init init_edges.txt]\n"
             << "    [--tabu T] [--iters N] [--quiet]\n"
             << "    [--max-parents M] [--max-children M]\n"
             << "    [--topk K] [--mi-sample S] [--mi-budget B]\n"
             << "    [--reach dense|lazy] [--jindex-cache C]\n";
        return 1;
    }

    std::string input_path;
    std::string score_name;
    ScoreType sc = ScoreType::BIC;

    // --- help オプション ---
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help();
            return 0;
        }
    }

    // ===== 既定値（D=10k〜30k を想定した省メモリ寄り） =====
    double ess = 1.0;
    std::string init_path;
    int tabu_tenure = 20;       // 0でタブー無効（純HC）
    int iters = 5000;
    bool verbose = true;
    int max_parents = 3, max_children = 16;
    int K = 50;                 // 候補親K
    int mi_sample = 8000;       // MI 行サンプル数（0で全行）
    int mi_budget = 1500;       // 1ノードあたり評価する相手変数数（0で全変数）
    std::string reach_mode_str = "lazy"; // 省メモリ既定
    int jindex_cache_cap = 1024;  // j_index LRU 容量（ノード数）
    std::string save_init_path;        // インデックスで保存
    std::string save_init_names_path;  // 変数名で保存（任意）
    std::string save_counts_path;   // カウントを保存する単一ファイル
    std::string score_dataset_path;  // 尤度を計算する新データセット
    std::string counts_in_path;      // all_counts.tsv の入力
    double alpha_ij = 0.0;      // 平滑化パラメータ（既定0=MLE）
    bool edge_importance_mode=false;
    std::string save_edge_importance_path;
    //
    int bootstrap_B = 0;                 // >0 でブートストラップ実行
    uint64_t seed = 2025;  // 既定シード
    std::string save_bootstrap_counts;   // 出力TSV
//    bool bootstrap_include_zero = false; // 0カウントも出すか
    std::string cand_metric_str = "mi";   // "mi" | "chi2"
    double mi_threshold = 0.0;           // MIの閾値：単位は ナチュラル対数（nats）
    double chi2_p_threshold = 1.0;       // 採用は p <= しきい値

    //
    //
    // 追加オプションのパース
    for (int i=1;i<argc;++i){
        std::string a=argv[i];
        auto need=[&](const char* f){return a==f && i+1<argc;};
        if (need("--input")) input_path = argv[++i];
        else if (need("--score")) score_name = argv[++i];
        else if (need("--ess")) ess = std::stod(argv[++i]);
        else if (need("--init")) init_path = argv[++i];
        else if (need("--tabu")) tabu_tenure = std::stoi(argv[++i]);
        else if (need("--iters")) iters = std::stoi(argv[++i]);
        else if (need("--max-parents")) max_parents = std::stoi(argv[++i]);
        else if (need("--max-children")) max_children = std::stoi(argv[++i]);
        else if (need("--topk")) K = std::stoi(argv[++i]);
        else if (need("--mi-sample")) mi_sample = std::stoi(argv[++i]);
        else if (need("--mi-budget")) mi_budget = std::stoi(argv[++i]);
        else if (need("--mi-threshold")) mi_threshold = std::stod(argv[++i]);
        else if (need("--cand-metric"))        cand_metric_str = argv[++i];  // "mi" or "chi2"
        else if (need("--chi2-p-threshold"))   chi2_p_threshold = std::stod(argv[++i]);

        else if (need("--reach")) reach_mode_str = argv[++i];
        else if (need("--jindex-cache")) jindex_cache_cap = std::stoi(argv[++i]);
        else if (need("--save"))        save_init_path = argv[++i];
        else if (need("--save-names"))  save_init_names_path = argv[++i];
        else if (need("--save-counts")) save_counts_path = argv[++i];
        else if (need("--score-dataset"))  score_dataset_path = argv[++i];
        else if (need("--counts"))         counts_in_path = argv[++i];
        else if (need("--alpha"))          alpha_ij = std::stod(argv[++i]);
        else if (need("--edge-importance")) edge_importance_mode = true;
        else if (need("--save-edge-importance")) save_edge_importance_path = argv[++i];
        else if (need("--bootstrap"))             bootstrap_B = std::stoi(argv[++i]);
        else if (need("--seed"))        seed = std::stoull(argv[++i]);
        else if (need("--save-bootstrap-counts")) save_bootstrap_counts = argv[++i];
        else if (a=="--quiet") verbose=false;
        else { std::cerr<<"Unknown option: "<<a<<"\n"; return 1; }
    }
    if(!score_name.empty()){
        sc=parseScore(score_name);
    }

    if (sc!=ScoreType::BDeu && ess!=1.0) std::cerr << "[warn] --ess is ignored for non-BDeu.\n";
    using clock = std::chrono::steady_clock;
    auto t_start = clock::now();
    
    // metric 解釈
    CandMetric cand_metric = CandMetric::MI;
    if (cand_metric_str == "chi2" || cand_metric_str == "CHI2") cand_metric = CandMetric::CHI2;
    // ===== ブートストラップ・モード（学習は行わず、B回の再学習＋集計のみ） =====
    if (bootstrap_B > 0) {
        if (save_bootstrap_counts.empty()) {
            std::cerr << "Error: --save-bootstrap-counts <path> is required with --bootstrap\n";
            return 2;
        }
        // 入力データ（CSV/TSV・ヘッダあり）
        Dataset ds = Dataset::fromCSV(input_path);

        const int D = ds.D;
        const int N = ds.N;
        const int* ds_x_ptr = ds.X_flat.data();
        const int* ds_r_ptr = ds.r.data();
        #pragma acc enter data copyin(ds_x_ptr[0:N*D],ds_r_ptr[0:D])

        // 既存のオプション（score, ess, init, tabu, iters, max-parents/children, topK, mi_sample, mi_budget, reach_mode, jindex_cache_cap）
        runBootstrapStructureCounts(
          ds, sc, ess, init_path,
          tabu_tenure, iters,
          max_parents, max_children,
          cand_metric,
          K, mi_sample, mi_budget, mi_threshold, chi2_p_threshold,
          (reach_mode_str=="dense"? Reachability::DENSE : Reachability::LAZY),
          jindex_cache_cap,
          bootstrap_B, seed,
          save_bootstrap_counts
        );
        #pragma acc exit data delete(ds_x_ptr[0:N*D],ds_r_ptr[0:D])
        // --- 終了時間計測 ---
        auto t_end = clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;
        double sec = elapsed.count();
        std::cerr << "[info] total runtime: " << sec << " seconds\n";

        return 0;
    }else if (edge_importance_mode) {
        try {
            Dataset ds_new = Dataset::fromCSV(score_dataset_path);

            const int D = ds_new.D;
            const int N = ds_new.N;
            const int* ds_new_x_ptr = ds_new.X_flat.data();
            const int* ds_new_r_ptr = ds_new.r.data();
            #pragma acc enter data copyin(ds_new_x_ptr[0:N*D],ds_new_r_ptr[0:D])

            DAG g_base = loadInitEdges(ds_new.D, init_path);
            std::vector<Counts> C_loaded = loadAllCountsTSV(counts_in_path, ds_new.D);

            computeEdgeImportanceScores(ds_new, g_base, C_loaded, alpha_ij, ess, save_edge_importance_path);

            #pragma acc exit data delete(ds_new_x_ptr[0:N*D],ds_new_r_ptr[0:D])
            // --- 終了時間計測 ---
            auto t_end = clock::now();
            std::chrono::duration<double> elapsed = t_end - t_start;
            double sec = elapsed.count();
            std::cerr << "[info] total runtime: " << sec << " seconds\n";
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "Error (edge-importance): " << e.what() << "\n";
            return 2;
        }
    }
    // 予測モード
    if (!score_dataset_path.empty()) {
        try {
            // 1) 新データセットのロード（CSV/TSV自動判定・ヘッダあり）
            Dataset ds_new = Dataset::fromCSV(score_dataset_path);

            const int D = ds_new.D;
            const int N = ds_new.N;
            const int* ds_new_x_ptr = ds_new.X_flat.data();
            const int* ds_new_r_ptr = ds_new.r.data();
            #pragma acc enter data copyin(ds_new_x_ptr[0:N*D],ds_new_r_ptr[0:D])

            // 2) 構造（DAG）のロード（TSV/スペース; 既存の loadInitEdges を利用）
            if (init_path.empty())
                throw std::runtime_error("--init <init_edges.tsv> is required for scoring mode");
            DAG g_new = loadInitEdges(ds_new.D, init_path);

            // 3) all_counts.tsv のロード（Counts 配列を復元）
            if (counts_in_path.empty())
                throw std::runtime_error("--counts <all_counts.tsv> is required for scoring mode");
            std::vector<int> child_r_from_counts;
            std::vector<Counts> C_loaded = loadAllCountsTSV(counts_in_path, ds_new.D, &child_r_from_counts);

            // 子の基数の簡易チェック（新データと counts の r_i が食い違っていないか）
            for (int v=0; v<ds_new.D; ++v){
                if (child_r_from_counts[v] > 0 && child_r_from_counts[v] != ds_new.r[v]) {
                    std::cerr << "[warn] child cardinality differs at node " << v
                         << " (counts r_i=" << child_r_from_counts[v]
                         << ", new-data r_i=" << ds_new.r[v] << ")\n";
                }
            }

            // 4) 尤度計算
            double avgN=0.0, avgND=0.0; long long zero_hits=0;
            double LL = computeLogLikelihoodOnDataset(ds_new, g_new, C_loaded, alpha_ij, &avgN, &avgND, &zero_hits);

            std::cout << std::fixed << std::setprecision(6);
            std::cout << "# scoring_mode\n";
            std::cout << "log_likelihood_total\t" << LL << "\n";
            std::cout << "log_likelihood_per_sample\t" << avgN << "\n";
            std::cout << "log_likelihood_per_variable\t" << avgND << "\n";
            std::cout << "alpha_ij\t" << alpha_ij << "\n";

            #pragma acc exit data delete(ds_new_x_ptr[0:N*D],ds_new_r_ptr[0:D])
            // --- 終了時間計測 ---
            auto t_end = clock::now();
            std::chrono::duration<double> elapsed = t_end - t_start;
            double sec = elapsed.count();
            std::cerr << "[info] total runtime: " << sec << " seconds\n";
            return 0;
        } catch (const std::exception& e){
            std::cerr << "Error (scoring): " << e.what() << "\n";
            return 2;
        }
    }
    // 単純な探索モード
    try{
        Dataset ds = Dataset::fromCSV(input_path);

        const int D = ds.D;
        const int N = ds.N;
        const int* ds_x_ptr = ds.X_flat.data();
        const int* ds_r_ptr = ds.r.data();
        #pragma acc enter data copyin(ds_x_ptr[0:N*D],ds_r_ptr[0:D])

        DAG init = loadInitEdges(ds.D, init_path);


        // --- 候補親Kの前処理（必要なときのみ） ---
        std::vector<std::vector<int>> topKlist;
        auto t_start_klist = clock::now();
        if (K>0 || mi_threshold>0.0 || chi2_p_threshold<1.0){
            std::mt19937_64 rng(seed);
            // 行サンプルを作る（mi_sample==0 なら全行）
            std::vector<int> rows;
            if (mi_sample>0 && mi_sample < ds.N){
                rows.resize(ds.N); iota(rows.begin(), rows.end(), 0);
                shuffle(rows.begin(), rows.end(), rng);
                rows.resize(mi_sample); sort(rows.begin(), rows.end());
            } else {
                rows.resize(ds.N); iota(rows.begin(), rows.end(), 0);
            }
            int budget = (mi_budget>0? mi_budget : (ds.D-1));
            //topKlist = MICandidates::compute(ds, K, budget, rows, rng, mi_threshold);
            topKlist = AssocCandidates::compute(
                ds, K, budget, rows, rng, cand_metric, mi_threshold, chi2_p_threshold
            );
            if (verbose) std::cerr << "[mi] topK="<<K<<" sample="<<rows.size()<<" budget="<<budget<<"\n";
        }
        {
            for(int i=0;i<topKlist.size();++i){
                std::cout << i << " ";
                for(int j=0;j<topKlist[i].size();++j){
                    std::cout << topKlist[i][j] << " ";
                }
            std::cout<<std::endl;
            }
        }

        // --- 終了時間計測 ---
        auto t_end_klist = clock::now();
        std::chrono::duration<double> elapsed_klist = t_end_klist - t_start_klist;
        double sec_klist = elapsed_klist.count();
        std::cerr << "[info] top-K list runtime: " << sec_klist << " seconds\n";

        // 到達性モード確定
        Reachability::Mode rmode = (reach_mode_str=="dense"? Reachability::DENSE : Reachability::LAZY);

        // 探索器のセットアップ
        auto t_start_iter = clock::now();
        HillClimber hc(ds, sc, ess, init, rmode, jindex_cache_cap);
        hc.max_iter = iters;
        hc.verbose = verbose;
        hc.tabu_tenure = tabu_tenure;
        hc.max_parents = max_parents;
        hc.max_children = max_children;
        if (K>0){ hc.candParents = move(topKlist); hc.topK = K; }

        auto [g, score, it] = hc.run(/*use_tabu=*/(tabu_tenure>0));

        // --- イテレーションあたり時間計測 ---
        auto t_end_iter = clock::now();
        std::chrono::duration<double> elapsed_iter = t_end_iter - t_start_iter;
        double sec_iter = elapsed_iter.count();
        std::cerr << "[info] time / iter: " << sec_iter/it << " seconds\n";
        
        // --- トータル終了時間計測 ---
        auto t_end = clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;
        double sec = elapsed.count();
        std::cerr << "[info] total time: " << sec << " seconds\n";


        // 結果出力
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "# learned_score="<<score<<"\n";
        // グラフ出力
        /*
        cout << "# edges (u v as u->v):\n";
        for (auto &e : g.edges()) {
            if (!ds.var_names.empty()) {
                cout << ds.var_names[e.first] << "\t" << ds.var_names[e.second] << "\n";
            } else {
                cout << e.first << "\t" << e.second << "\n";
            }
        }
        */
        
        // ===== 次回初期化用にファイル保存（インデックス版） =====
        if (!save_init_path.empty()) {
            std::ofstream fout(save_init_path);
            if (!fout) {
                std::cerr << "Error: failed to open --save path: " << save_init_path << "\n";
                return 2;
            }
            for (auto &e : g.edges()) {
                fout << e.first << "\t" << e.second << "\n";  // 0始まり u v 形式
            }
            std::cerr << "[info] wrote init edges (indices) to " << save_init_path << "\n";
        }

        // =====  参考用に変数名でも保存（再利用不可・人間可読） =====
        if (!save_init_names_path.empty()) {
            std::ofstream foutn(save_init_names_path);
            if (!foutn) {
                std::cerr << "Error: failed to open --save-names path: " << save_init_names_path << "\n";
                return 2;
            }
            for (auto &e : g.edges()) {
                if (!ds.var_names.empty()) {
                    foutn << ds.var_names[e.first] << "\t" << ds.var_names[e.second] << "\n";
                } else {
                    // ヘッダが無い場合はインデックスでフォールバック
                    foutn << e.first << "\t" << e.second << "\n";
                }
            }
            std::cerr << "[info] wrote init edges (names) to " << save_init_names_path << "\n";
        }
        // 全ノード分のカウントを 1 ファイル（TSV）にまとめて出力
        if (!save_counts_path.empty()) {
            try {
                //saveAllCountsTSV(save_counts_path, ds, g, hc.nodeCounts);
                saveAllCountsTSV(save_counts_path, ds, g);
                std::cerr << "[info] wrote all CPT counts (TSV) to " << save_counts_path << "\n";
            } catch (const std::exception& e) {
                std::cerr << "Error: cannot write --save-counts: " << e.what() << "\n";
                return 2;
            }
        }
        #pragma acc exit data delete(ds_x_ptr[0:N*D],ds_r_ptr[0:D])

    } catch (const std::exception& e){
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }
    return 0;
}


#pragma once

#include <cmath>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <filesystem>

#ifdef __NVCOMPILER
#include <openacc.h>
#endif

#include "fast_bn_counts.hpp"

struct Scorer {
    const Dataset& ds;
    ScoreType type;
    double ess; // BDeu の等価事例数（他スコアでは未使用）

    Scorer(const Dataset& ds, ScoreType t, double ess=1.0): ds(ds), type(t), ess(ess) {}

    // BIC: 対数尤度 - (d/2)*log(N), d=(r_i-1)*q_i
    double nodeScoreBIC(const Counts& c) const {
        double ll=0.0;
        const long long* _nij  = c.n_ij.data();
        const long long* _nijk = c.n_ijk.data();
        int q_i = c.q_i;
        int r_i = c.r_i;
//        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
//            #pragma acc parallel loop reduction(+:ll) present(_nij[0:q_i],_nijk[0:q_i*r_i])
            for (int j=0;j<q_i;++j){
                double nij=(double)_nij[j];
                if (nij>0){
                    for (int k=0;k<r_i;++k){
                        long long nijk=_nijk[(size_t)j*r_i+k];
                        if (nijk>0){
                            ll += nijk * (log((double)nijk) - log(nij));
                        }
                    }
                }
            }
        }
        int d=(c.r_i-1)*c.q_i;
        double pen = 0.5 * d * log((double)std::max(1, ds.N));
        return ll - pen;
    }

    // K2: Dirichlet(1) 事前
    double nodeScoreK2(const Counts& c) const {
        double s=0.0;
        const long long* _nij  = c.n_ij.data();
        const long long* _nijk = c.n_ijk.data();
        int q_i = c.q_i;
        int r_i = c.r_i;
//        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
//            #pragma acc parallel loop reduction(+:s) present(_nij[0:q_i],_nijk[0:q_i*r_i])
            for (int j=0;j<q_i;++j){
                double nij=(double)_nij[j];
                s += lgamma((double)r_i) - lgamma(nij + (double)r_i);
                for (int k=0;k<r_i;++k){
                    double nijk=(double)_nijk[(size_t)j*r_i + k];
                    s += lgamma(nijk + 1.0); // - lgamma(1) = 0
                }
            }
        }
        return s;
    }

    // BDeu: 一様ハイパーパラメータ（等価事例数 ess を q_i, r_i に均等割）
    double nodeScoreBDeu(const Counts& c) const {
        double s=0.0;
        double alpha_ij = ess / (double)std::max(1,c.q_i);
        double alpha_ijk = alpha_ij / (double)c.r_i;
        const long long* _nij  = c.n_ij.data();
        const long long* _nijk = c.n_ijk.data();
        int q_i = c.q_i;
        int r_i = c.r_i;
//        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
//            #pragma acc parallel loop reduction(+:s) present(_nij[0:q_i],_nijk[0:q_i*r_i])
            for (int j=0;j<q_i;++j){
                double nij=(double)_nij[j];
                s += lgamma(alpha_ij) - lgamma(nij + alpha_ij);
                for (int k=0;k<r_i;++k){
                    double nijk=(double)_nijk[(size_t)j*r_i + k];
                    s += lgamma(nijk + alpha_ijk) - lgamma(alpha_ijk);
                }
            }
        }
        return s;
    }

    double nodeScore(const Counts& c) const {
//        c.check_gpu("nodeScore");
        switch(type){
            case ScoreType::BIC: return nodeScoreBIC(c);
            case ScoreType::K2:  return nodeScoreK2(c);
            case ScoreType::BDeu:return nodeScoreBDeu(c);
        }
        return -INFINITY;
    }
};

//==================== j_index LRU キャッシュ ====================

/*
  j_index[v][n] = 「現在の親集合 Pa(v) に基づくサンプル n の親配置インデックス j」
  - ADD 試行時に j' = j * r[u] + x_u[n] で O(1) 更新できるため、O(N) で新カウントを分割生成可能。
  - ただし全ノード v について N 要素の配列を常駐させるとメモリが大きいので、
    LRU で保持ノード数を制限（--jindex-cache）。
*/
struct JIndexCache {
    struct Entry { int v; std::vector<int> j; std::list<int>::iterator it; };
    int cap; // 最大保持ノード数（0 ならキャッシュ無効）
    const Dataset* ds;
    const DAG* g;
    std::unordered_map<int, std::unique_ptr<Entry>> mp;
    std::list<int> lru; // front が最新

    JIndexCache(int cap, const Dataset* ds, const DAG* g): cap(cap), ds(ds), g(g) {}

    void touch_(int v){
        if (!mp.count(v)) return;
        lru.erase(mp[v]->it);
        lru.push_front(v);
        mp[v]->it = lru.begin();
    }
    void evictIfNeeded_(){
        while (cap>0 && (int)mp.size()>cap){
            int victim = lru.back(); lru.pop_back();
            mp.erase(victim);
        }
    }
    void invalidate(int v){
        if (cap==0) return;
        if (!mp.count(v)) return;
        lru.erase(mp[v]->it);
        mp.erase(v);
    }
    void onParentsChanged(const std::vector<int>& vs){ if (cap==0) return; for (int v: vs) invalidate(v); }

    // 親集合 Pa(v) に合わせて j_index を構築（O(N·|Pa|)）
    void build(int v, std::vector<int>& out) const {
        auto pa = g->parents(v);
        const int N = ds->N;
        out.assign(N, 0);
        if (pa.empty()) return;
        std::vector<int> radix;
        build_mixed_radix(pa, *ds, radix);
        const int P = (int)pa.size();
        const int D = ds->D;
        const int* pa_ptr = pa.data();
        const int* ds_x_ptr = ds->X_flat.data();
        int* rdx_ptr = radix.data();
        int* out_ptr = out.data();
        #pragma acc data present(ds_x_ptr[0:N*D])
        {
            #pragma acc parallel loop present(ds_x_ptr[0:N*D]) copy(pa_ptr[0:P],rdx_ptr[0:P])
            for (int n=0;n<N;++n){
                const int j = mixed_radix_index_row(D, ds_x_ptr, n, P, pa_ptr, rdx_ptr);
                out_ptr[n]=j;
            }
        }
    }

    // 取得：キャッシュにあれば返す。無ければ構築してキャッシュ。
    const std::vector<int>& get(int v){
        if (cap==0){
            // 省メモリモード：都度 build して一時オブジェクトを返す（実用では cap>0 を推奨）
            static std::vector<int> tmp;
            build(v, tmp);
            // 注意: 本簡易実装では cap==0 の場合にメモリリークを許容（サンプル用途）
            return *(new std::vector<int>(tmp));
        }
        if (mp.count(v)) { touch_(v); return mp[v]->j; }
        auto e = std::make_unique<Entry>();
        e->v = v;
        build(v, e->j);
        lru.push_front(v);
        e->it = lru.begin();
        mp[v] = std::move(e);
        evictIfNeeded_();
        return mp[v]->j;
    }
};

//==================== 候補親Kの前処理（MI: 相互情報量） ====================

/*
  目的:
    - 各ノード v に対して、相互情報量 MI(u; v) が大きい上位 K 個の u を候補親に限定
    - 検索空間を大幅に縮小

  計算量削減の工夫:
    - --mi-sample S: 行サンプル数 S のみで MI を近似（0 なら全行）
    - --mi-budget B: 1ノードあたり B 変数のみを評価（0 なら全変数）
*/
enum class CandMetric { MI, CHI2 };

// rows: サンプリングした行インデックス（全件なら 0..N-1）
struct AssocCandidates {
    // MI（nats）
    static double mutual_info_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows);

    // カイ二乗の p 値（独立性検定、上側確率）と統計量
    static std::pair<double,double> chi2_p_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows);

    // metric に応じて候補親を返す
    //  - metric=MI    : mi >= mi_threshold_nats を残し、K>0なら MI 降順で上位K
    //  - metric=CHI2  : p <= chi2_p_threshold を残し、K>0なら -log(p) 降順で上位K
    static std::vector<std::vector<int>> compute(const Dataset& ds,
                                                 int K,
                                                 int budget,
                                                 const std::vector<int>& rows,
                                                 std::mt19937_64& rng,
                                                 CandMetric metric,
                                                 double mi_threshold_nats,
                                                 double chi2_p_threshold);
};

struct MICandidates {
    static double mutual_info_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows);

    static std::vector<std::vector<int>> compute(const Dataset& ds, int K, int budget, const std::vector<int>& rows, std::mt19937_64& rng, double mi_threshold_nats);
};

//==================== 探索（HC/Tabu） ====================

struct HillClimber {
    const Dataset& ds;
    Scorer scorer;
    DAG g;
    Reachability reach;
    int max_iter = 100000;
    bool verbose = true;

    // 構造制約（次数上限）：q を小さく保つため重要（メモリ・計算量の上限に効く）
    int max_parents = 3;   // 既定: 大規模向け
    int max_children = 16;

    // タブー
    int tabu_tenure = 0;

    // 候補親（MI の上位 K）
    std::vector<std::vector<int>> candParents; // ソート済ベクタ
    int topK = 50;

    // j_index LRU キャッシュ
    JIndexCache jcache;

    // 現在のローカルカウント & スコア（各ノードで保持）
    // - ADD: 分割で新カウントを生成 → 更新
    // - REMOVE: **完全インクリメンタル**にマージ → 更新
    std::vector<Counts> nodeCounts;
    std::vector<double> nodeScoreNow;
    double totalNow = 0.0;
    Counts tmpCounts_a = {
      std::vector<long long>(4096),
      std::vector<long long>(128),
      128,
      32
    };
    Counts tmpCounts_r = {
      std::vector<long long>(4096),
      std::vector<long long>(128),
      128,
      32
    };


    HillClimber(const Dataset& ds, ScoreType t, double ess,
                const DAG& init, Reachability::Mode rmode,
                int jcache_cap)
        : ds(ds), scorer(ds,t,ess), g(init), reach(init, rmode), jcache(jcache_cap, &ds, &g)
    {
        nodeCounts.reserve(ds.D);
        nodeScoreNow.reserve(ds.D);

        // 初期カウントを各ノードで構築（O(N·|Pa|) を D 回）
        for (int v=0; v<ds.D; ++v){
            auto Pa = g.parents(v);
            std::vector<int> radix;
            build_mixed_radix(Pa, ds, radix);
            nodeCounts.emplace_back(computeCountsForNode_full(v, Pa, ds, radix));
            nodeScoreNow.emplace_back(scorer.nodeScore(nodeCounts[v]));
            totalNow += nodeScoreNow[v];
            nodeCounts[v].check_gpu();
        }
    }

    struct Move { enum Type { ADD, REMOVE, REVERSE, NONE } type=NONE; int u=-1, v=-1; double delta=0; };

    // 追加時の次数上限チェック
    bool addAllowedDegree(int u,int v) const {
        if (g.parentCount(v) >= max_parents) return false;
        if (g.childCount(u)  >= max_children) return false;
        return true;
    }

    bool addAllowed(const int u, const int v) const {
        if (!addAllowedDegree(u,v)) return false;
        if (!reach.canAddAcyclic(g,u,v)) return false;
        if (topK>0){
            const auto& C = candParents[v];
            if (!binary_search(C.begin(), C.end(), u)) return false;
        }
        return true;
    }

    bool reverseAllowed(int u,int v) const {
        if (!reach.canReverseAcyclic(g,u,v)) return false;
        if (g.parentCount(u) + 1 > max_parents) return false; // v->u の追加で u の親+1
        if (g.childCount(v)  + 1 > max_children) return false; // v の子+1
        if (topK>0){
            const auto& C = candParents[u];
            if (!binary_search(C.begin(), C.end(), v)) return false; // 逆向きも候補に含まれるか
        }
        return true;
    }

    //==================== Δスコア計算（インクリメンタル） ====================

    // ADD(u->v): 既存の j_index[v]（Pa(v) 用）を使い、q' = q * r_u の配列に「分割」して再構成（O(N)）
    // 返値: Δスコア（after - before）。生成した newC は呼び出し側で適用に使う。
    double deltaAdd_andBuildNewCounts(int u, int v, Counts& newC);

    // REMOVE(u->v): **完全インクリメンタル**
    // - 現在の counts（Pa(v)）から、削除対象 u の「桁」を落として合算するだけ（O(q·r_i)）
    // - データの再走査無し
    double deltaRemove_andBuildNewCounts(int u, int v, const std::vector<int> &Pa, Counts& newC);

    // REVERSE(u->v): v 側は REMOVE（マージ）、u 側は ADD（分割）
//    double deltaReverse_buildNewCounts(int u, int v, const std::vector<int> &Pa, Counts& newC_v, Counts& newC_u, double& d_v, double& d_u){
//        d_v = deltaRemove_andBuildNewCounts(u, v, Pa, newC_v); // v: 親 u を外す（マージ）
//        d_u = deltaAdd_andBuildNewCounts(v, u, newC_u);    // u: 親に v を加える（分割）
//        return d_v + d_u;
//    }

    //==================== 実行（HC / Tabu） ====================

    std::tuple<DAG,double,int> run(bool use_tabu){
        double cur = totalNow;
        DAG bestG = g; double bestScore = cur;

        if (verbose) std::cerr << std::fixed << std::setprecision(6)
                          << "[start] score="<<cur<<" edges="<<g.edges().size()
                          << " mode="<<(use_tabu?"tabu":"greedy")<<"\n";

        tmpCounts_a.check_gpu("run a");
        tmpCounts_r.check_gpu("run r");

        const int D = ds.D;
        // タブー配列（属性タブー：直近操作の巻き戻しを禁止）
        std::vector<std::vector<int>> tabu_add_until(D, std::vector<int>(D,-1));
        std::vector<std::vector<int>> tabu_remove_until(D, std::vector<int>(D,-1));
        std::vector<std::vector<int>> tabu_reverse_until(D, std::vector<int>(D,-1));

        auto isTabu = [&](const Move& mv, int it)->bool{
            if (!use_tabu) return false;
            if (mv.type==Move::ADD)    return tabu_add_until[mv.u][mv.v] > it;
            if (mv.type==Move::REMOVE) return tabu_remove_until[mv.u][mv.v] > it;
            if (mv.type==Move::REVERSE)return tabu_reverse_until[mv.u][mv.v] > it;
            return false;
        };
        auto setTabuAfter = [&](const Move& mv, int it){
            if (!use_tabu) return;
            int until = it + tabu_tenure;
            if (mv.type==Move::ADD){
                // add の直後：remove と reverse（戻し）を禁止
                tabu_remove_until[mv.u][mv.v] = until;
                tabu_reverse_until[mv.u][mv.v] = until;
            } else if (mv.type==Move::REMOVE){
                // remove の直後：add（元に戻す）を禁止
                tabu_add_until[mv.u][mv.v] = until;
            } else if (mv.type==Move::REVERSE){
                // reverse の直後：逆向き reverse（巻き戻し）を禁止
                tabu_reverse_until[mv.v][mv.u] = until;
            }
        };

        int it=0;
        for (; it<max_iter; ++it){
            Move best; best.type=Move::NONE; best.delta=use_tabu?-1e300:0.0;
            Move bestNonTabu; bestNonTabu.type=Move::NONE; bestNonTabu.delta=-1e300;
            std::cout<<"D="<<D<<std::endl;
            // 近傍列挙（ADD は候補親 K のみ、REMOVE/REVERSE は現辺）
            // Counts _newC,_newCu,_newCv;
            for (int v=0; v<D; ++v){
                // --- ADD 候補 ---
                if (topK>0){
                    for (int u: candParents[v]){
                        if (!addAllowed(u,v)) continue;
//                        Counts newC;
                        double d = deltaAdd_andBuildNewCounts(u, v, tmpCounts_a);
                        Move mv{Move::ADD,u,v,d};
                        bool tabu = isTabu(mv, it);
                        bool asp = (cur + d > bestScore + 1e-12); // アスピレーション
                        if (!tabu) { if (d > bestNonTabu.delta) bestNonTabu = mv; }
                        if (!tabu || asp) { if (d > best.delta) best = mv; }
                    }
                } else {
                    for (int u=0; u<D; ++u) if (u!=v){
                        if (!addAllowed(u,v)) continue;
//                        Counts newC;
                        double d = deltaAdd_andBuildNewCounts(u, v, tmpCounts_a);
                        Move mv{Move::ADD,u,v,d};
                        bool tabu = isTabu(mv, it);
                        bool asp = (cur + d > bestScore + 1e-12);
                        if (!tabu) { if (d > bestNonTabu.delta) bestNonTabu = mv; }
                        if (!tabu || asp) { if (d > best.delta) best = mv; }
                    }
                }

                // --- REMOVE 候補（現にある辺）---
                std::vector<int> Pa = g.parents(v);
                const int* _pa = Pa.data();
                int pa_size = Pa.size();
                #pragma acc enter data copyin(_pa[0:pa_size])
                for (int u=0; u<D; ++u) if (g.adj[u][v]){
//                    Counts newC;
                    double d = deltaRemove_andBuildNewCounts(u, v, Pa, tmpCounts_r);
                    Move mv{Move::REMOVE,u,v,d};
                    bool tabu = isTabu(mv, it);
                    bool asp = (cur + d > bestScore + 1e-12);
                    if (!tabu) { if (d > bestNonTabu.delta) bestNonTabu = mv; }
                    if (!tabu || asp) { if (d > best.delta) best = mv; }
                }

                // --- REVERSE 候補（現にある辺で逆向き追加が可能なもの）---

                for (int u=0; u<D; ++u) if (g.adj[u][v]){
                    if (!reverseAllowed(u,v)) continue;
//                    Counts newCv, newCu;
                    double dv = deltaRemove_andBuildNewCounts(u, v, Pa, tmpCounts_r); // v: 親 u を外す（マージ）
                    double du = deltaAdd_andBuildNewCounts(v, u, tmpCounts_a);    // u: 親に v を加える（分割）
                    double d = du + dv;
                    Move mv{Move::REVERSE,u,v,d};
                    bool tabu = isTabu(mv, it);
                    bool asp = (cur + d > bestScore + 1e-12);
                    if (!tabu) { if (d > bestNonTabu.delta) bestNonTabu = mv; }
                    if (!tabu || asp) { if (d > best.delta) best = mv; }
                }
                #pragma acc exit data delete(_pa[0:pa_size])
            }

            // タブーによる最良不可→非タブー最良を選ぶ
            Move chosen = best;
            if (chosen.type==Move::NONE && use_tabu && bestNonTabu.type!=Move::NONE) chosen = bestNonTabu;
            if (chosen.type==Move::NONE || (!use_tabu && chosen.delta<=1e-12)){
                if (verbose) std::cerr << "[stop] no improving move.\n";
                break;
            }

            // ===== 実適用（DAG / reach / j_index / counts / scores を整合させる） =====
            if (chosen.type==Move::ADD){
                double d = deltaAdd_andBuildNewCounts(chosen.u, chosen.v, tmpCounts_a);
                tmpCounts_a.check_gpu("run ADD");
                g.addEdge(chosen.u, chosen.v);
                reach.onAddEdge(g, chosen.u, chosen.v);
                jcache.onParentsChanged({chosen.v});          // v の親が変わったので j_index[v] を無効化
                std::cout << "ADD (q,r) " <<
                    "(" << nodeCounts[chosen.v].q_i << "," << nodeCounts[chosen.v].r_i << ") => " <<
                    "(" << tmpCounts_a.q_i << "," << tmpCounts_a.r_i << ")" << std::endl;
                nodeCounts[chosen.v] = tmpCounts_a;            // v の counts を差し替え
                nodeCounts[chosen.v].check_gpu("run ADD commit");
                nodeScoreNow[chosen.v] += d;                  // v のスコアだけ増分更新
                totalNow += d;                                 // 全体スコア更新

            } else if (chosen.type==Move::REMOVE){
                std::vector<int> Pa = g.parents(chosen.v);
                const int* _pa = Pa.data();
                int pa_size = Pa.size();
                #pragma acc enter data copyin(_pa[0:pa_size])
                double d = deltaRemove_andBuildNewCounts(chosen.u, chosen.v, Pa, tmpCounts_r);
                tmpCounts_r.check_gpu("run REMOVE");
                g.removeEdge(chosen.u, chosen.v);
                reach.onRemoveEdge(g, chosen.u, chosen.v);
                jcache.onParentsChanged({chosen.v});
                std::cout << "REMOVE (q,r) " <<
                    "(" << nodeCounts[chosen.v].q_i << "," << nodeCounts[chosen.v].r_i << ") => " <<
                    "(" << tmpCounts_r.q_i << "," << tmpCounts_r.r_i << ")" << std::endl;
                nodeCounts[chosen.v] = tmpCounts_r;
                nodeCounts[chosen.v].check_gpu("run REMOVE commit");
                nodeScoreNow[chosen.v] += d;
                totalNow += d;
                #pragma acc exit data delete(_pa[0:pa_size])
            } else if (chosen.type==Move::REVERSE){
                std::vector<int> Pa = g.parents(chosen.v);
                const int* _pa = Pa.data();
                int pa_size = Pa.size();
                #pragma acc enter data copyin(_pa[0:pa_size])
                double dv = deltaRemove_andBuildNewCounts(chosen.u, chosen.v, Pa, tmpCounts_r); // v: 親 u を外す（マージ）
                double du = deltaAdd_andBuildNewCounts(chosen.v, chosen.u, tmpCounts_a);    // u: 親に v を加える（分割）
                tmpCounts_a.check_gpu("run REVERSE");
                tmpCounts_r.check_gpu("run REVERSE");
                double d = du + dv;
                g.reverseEdge(chosen.u, chosen.v);
                reach.onReverseEdge(g, chosen.u, chosen.v);
                jcache.onParentsChanged({chosen.u, chosen.v});
                std::cout << "REVERSE (q,r) " <<
                    "(" << nodeCounts[chosen.v].q_i << "," << nodeCounts[chosen.v].r_i << ") => " <<
                    "(" << tmpCounts_r.q_i << "," << tmpCounts_r.r_i << ")" << std::endl;
                std::cout << "REVERSE (q,r) " <<
                    "(" << nodeCounts[chosen.u].q_i << "," << nodeCounts[chosen.u].r_i << ") => " <<
                    "(" << tmpCounts_a.q_i << "," << tmpCounts_a.r_i << ")" << std::endl;
                nodeCounts[chosen.v] = tmpCounts_r;
                nodeCounts[chosen.u] = tmpCounts_a;
                nodeCounts[chosen.v].check_gpu("run REVERSE commit");
                nodeCounts[chosen.u].check_gpu("run REVERSE commit");
                nodeScoreNow[chosen.v] += dv;
                nodeScoreNow[chosen.u] += du;
                totalNow += d;
                #pragma acc exit data delete(_pa[0:pa_size])
            }

            setTabuAfter(chosen, it);

            if (totalNow > bestScore + 1e-12){ bestScore=totalNow; bestG=g; if (verbose) std::cerr<<"[*] new best "<<bestScore<<"\n"; }
            if (verbose){
                std::string t = (chosen.type==Move::ADD?"ADD": chosen.type==Move::REMOVE?"REM":"REV");
                std::cerr << "[it "<<it+1<<"] "<<t<<" "<<chosen.u<<"->"<<chosen.v
                     <<"  delta="<<chosen.delta<<"  cur="<<totalNow<<std::endl;
            }
        }
        return {bestG, bestScore, it};
    }
};

//==================== 引数処理 / 初期構造ロード / MI 前処理 ====================

static ScoreType parseScore(const std::string& s){
    std::string t=s; for (auto& c:t) c=(char)tolower(c);
    if (t=="bic") return ScoreType::BIC;
    if (t=="k2")  return ScoreType::K2;
    if (t=="bdeu")return ScoreType::BDeu;
    throw std::runtime_error("Unknown score: "+s+" (use: bic|k2|bdeu)");
}

// 初期構造（エッジリスト）を読み込む
// ・タブ区切り（TSV）およびスペース区切りに対応
// ・コメント行 (#で始まる) は無視
static DAG loadInitEdges(int D, const std::string& path) {
    DAG g(D);
    if (path.empty()) return g;

    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Failed to open init edge list: " + path);

    std::string line;
    int line_no = 0;
    while (getline(fin, line)) {
        ++line_no;
        if (line.empty()) continue;
        if (line[0] == '#') continue; // コメント行スキップ

        // 区切り文字自動判定（タブ or スペース）
        char delim = (line.find('\t') != std::string::npos) ? '\t' : ' ';

        std::stringstream ss(line);
        std::string u_str, v_str;
        if (!getline(ss, u_str, delim)) continue;
        if (!getline(ss, v_str, delim)) continue;

        // 前後空白除去
        auto trim = [](std::string& s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char c){return !isspace(c);} ));
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c){return !isspace(c);} ).base(), s.end());
        };
        trim(u_str);
        trim(v_str);

        if (u_str.empty() || v_str.empty()) continue;

        int u = std::stoi(u_str);
        int v = std::stoi(v_str);
        if (u < 0 || u >= D || v < 0 || v >= D || u == v)
            throw std::runtime_error("Invalid edge at line " + std::to_string(line_no) + ": " + u_str + " " + v_str);

        if (!g.adj[u][v]) g.addEdge(u, v);
    }

    std::cerr << "[info] loaded init edges from " << path << " (" << g.edges().size() << " edges)\n";
    return g;
}

// ================================================
// 最終構造に基づく全ノード分のカウントを 1 つの TSV に出力
// - outfile: 出力パス（TSV）
// - ds     : データセット（基数・変数名に使用）
// - g      : 最終 DAG（親集合取得に使用）
// - counts : 各ノードの最終カウント（nodeCounts[v]）
// 形式：タブ区切り（TSV）
//   ・各ノードの前にメタ情報を # で出力（親の並び・基数・変数名など）
//   ・本体の行は 2 種類：
//        v <tab> j <tab> k <tab> n_ijk
//        v <tab> j <tab> * <tab> n_ij     （* は合計を表す）
//   ・n=0 の行は省略（ファイルサイズ節約）
// ================================================
inline void saveAllCountsTSV(const std::string& path, const Dataset& ds, const DAG& g) {
    std::ofstream fout(path);
    if (!fout) throw std::runtime_error("open failed: " + path);

    for (int v = 0; v < ds.D; ++v) {
        auto pa = g.parents(v);

        // メタ情報（親が空なら空行を出す）
        fout << "# --- node " << v << " ---\n";
        fout << "# node_name\t" << (v < (int)ds.var_names.size() ? ds.var_names[v] : ("X"+std::to_string(v))) << "\n";
        fout << "# parents_indices\t";
        for (size_t t=0; t<pa.size(); ++t) fout << (t? ",":"") << pa[t];
        fout << "\n# parents_names\t";
        for (size_t t=0; t<pa.size(); ++t) {
            int p = pa[t];
            fout << (t? ",":"") << (p < (int)ds.var_names.size() ? ds.var_names[p] : ("X"+std::to_string(p)));
        }
        fout << "\n# parents_cardinalities\t";
        for (size_t t=0; t<pa.size(); ++t) fout << (t? ",":"") << ds.r[pa[t]];
        fout << "\n# child_cardinality\t" << ds.r[v] << "\n";

        // 再度カウント（q_i は親配置数）
        std::vector<int> radix;
        build_mixed_radix(pa, ds, radix);
        Counts C = computeCountsForNode_full(v, pa, ds, radix);
        const int r_i = C.r_i;         // 子の取りうる値の数
        int q_i = C.q_i;               // 親配置数
        if (pa.empty()) q_i = 1;       // 念のため明示

        // j は 0..q_i-1 で回すのが正しい
        for (int j = 0; j < q_i; ++j) {
            long long nij = C.n_ij[j];

            // 各 k 行（ゼロは省略
            for (int k = 0; k < r_i; ++k) {
                long long nijk = C.n_ijk[(size_t)j*r_i + k];
                if (nijk == 0) continue;  // 省略したくない場合は消す
                fout << v << "\t" << j << "\t" << k << "\t" << nijk << "\n";
            }
            // 合計行は必ず出す
            fout << v << "\t" << j << "\t*\t" << nij << "\n";
        }
    }
}

// ============ all_counts.tsv ロード ============
// all_counts.tsv（単一ファイル）から、各ノード v の Counts を復元する。
// フォーマット：データ行 "v \t j \t k \t n" （k='*' のとき n_ij）
static std::vector<Counts> loadAllCountsTSV(const std::string& path, int expected_D, std::vector<int>* out_child_r = nullptr) {
    std::cout << "loadAllCountsTSV " << path << std::endl;
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Failed to open all-counts file: " + path);

    // 一時的に「ノードごと」に (j,k,n_ijk) と (j, n_ij) を保持
    struct KEntry { int j; int k; long long n; };
    struct JEntry { int j; long long n; };
    std::vector<std::vector<KEntry>> tmpK; // tmpK[v] に (j,k,n)
    std::vector<std::vector<JEntry>> tmpJ; // tmpJ[v] に (j,n)
    tmpK.resize(expected_D);
    tmpJ.resize(expected_D);

    // サイズ推定用
    std::vector<int> maxJ(expected_D, -1), maxK(expected_D, -1);

    std::string line;
    long long data_lines = 0;
    while (getline(fin, line)) {
        if (line.empty() || line[0]=='#') continue;

        // v \t j \t k \t n
        // k は整数か '*'（合計）を取る
        std::string vstr, jstr, kstr, nstr;
        {
            std::stringstream ss(line);
            if (!getline(ss, vstr, '\t')) continue;
            if (!getline(ss, jstr, '\t')) continue;
            if (!getline(ss, kstr, '\t')) continue;
            if (!getline(ss, nstr, '\t')) continue;
        }
        // trim 簡易
        auto trim = [](std::string& s){
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char c){return !isspace(c);} ));
            s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c){return !isspace(c);} ).base(), s.end());
        };
        trim(vstr); trim(jstr); trim(kstr); trim(nstr);

        int v = std::stoi(vstr);
        int j = std::stoi(jstr);
        if (v<0 || v>=expected_D) throw std::runtime_error("counts file: node index out of range: "+std::to_string(v));
        long long n = std::stoll(nstr);

        if (kstr == "*" || kstr == "'*'") {
            // n_ij
            tmpJ[v].push_back({j, n});
            maxJ[v] = std::max(maxJ[v], j);
        } else {
            int k = stoi(kstr);
            tmpK[v].push_back({j, k, n});
            maxJ[v] = std::max(maxJ[v], j);
            maxK[v] = std::max(maxK[v], k);
        }
        ++data_lines;
    }
    if (data_lines==0) throw std::runtime_error("counts file has no data lines: " + path);

    // 復元（配列確保→詰め込み）
    std::vector<Counts> C(expected_D);
    std::vector<int> child_r_detected(expected_D, 0);

    for (int v=0; v<expected_D; ++v){
        int q = (maxJ[v] >= 0 ? (maxJ[v]+1) : 1); // 観測なしなら q=1 とみなす
        int r = (maxK[v] >= 0 ? (maxK[v]+1) : 1); // 観測なしなら r=1 とみなす

        C[v].q_i = q;
        C[v].r_i = r;
        C[v].n_ij.assign(q, 0);
        C[v].n_ijk.assign((size_t)q * r, 0);

        for (auto &e : tmpJ[v]) {
            if (e.j < 0 || e.j >= q) throw std::runtime_error("n_ij j out of range at node "+std::to_string(v));
            C[v].n_ij[e.j] += e.n;
        }
        for (auto &e : tmpK[v]) {
            if (e.j < 0 || e.j >= q) throw std::runtime_error("n_ijk j out of range at node "+std::to_string(v));
            if (e.k < 0 || e.k >= r) throw std::runtime_error("n_ijk k out of range at node "+std::to_string(v));
            C[v].n_ijk[(size_t)e.j * r + e.k] += e.n;
        }
        child_r_detected[v] = r;
    }

    if (out_child_r) *out_child_r = std::move(child_r_detected);
    std::cout << "complete loadAllCountsTSV" << std::endl;
    return C;
}

// ============ 新データセットに対する対数尤度 ============
// alpha_ij: 平滑化（Dirichlet 事前）パラメータ a
//   P(x_i=k | pa=j) = (n_ijk + a/r_i) / (n_ij + a)
//   a=0 なら MLE（ゼロ割・確率0→ -inf に注意）
static double computeLogLikelihoodOnDataset(const Dataset& ds_new,
                                            const DAG& g,
                                            const std::vector<Counts>& C,
                                            double alpha_ij,
                                            double* out_avg_per_sample = nullptr,
                                            double* out_avg_per_var = nullptr,
                                            long long* out_zero_hits = nullptr)
{
    if (g.D != ds_new.D) throw std::runtime_error("D mismatch: graph vs dataset");
    if ((int)C.size() != g.D) throw std::runtime_error("Counts vector size != D");

    const int N = ds_new.N;
    const int D = ds_new.D;

    // 親集合は g.parents(v) の順序を使用（counts の j もこの順序で作っている前提）
    double LL = 0.0;
    long long zero_hits = 0;

    // 事前のための定数（r_i は C[v].r_i を使う）
    for (int n=0; n<N; ++n){
        for (int v=0; v<D; ++v){
            const auto parents = g.parents(v);
            int j = 0;
            if (!parents.empty()){
                // 混合基数（右端の親が最下位桁）で j を計算
                int mult = 1;
                for (int t=(int)parents.size()-1; t>=0; --t){
                    int p = parents[t];
                    j += ds_new.x(n,p) * mult;
                    mult *= ds_new.r[p];
                }
            }
            int k = ds_new.x(n,v);
            const Counts& Cv = C[v];
            if (j<0 || j>=Cv.q_i) {
                // counts に存在しない親配置（例えば子や親の基数が異なる等）
                // → smoothing がゼロなら -inf、a>0 なら一様事前で評価
                if (alpha_ij<=0.0) return -INFINITY;
                double prob = (1.0 / (double)Cv.r_i); // n_ij=0 前提で (0+a/r)/(0+a)
                LL += log(prob);
                continue;
            }

            long long nij = Cv.n_ij[j];
            long long nijk = (k>=0 && k<Cv.r_i) ? Cv.n_ijk[(size_t)j*Cv.r_i + k] : 0;

            double prob;
            if (alpha_ij > 0.0){
                prob = ( (double)nijk + alpha_ij / (double)Cv.r_i ) / ( (double)nij + alpha_ij );
            } else {
                if (nij == 0 || nijk == 0) {
                    ++zero_hits;
                    return -INFINITY; // MLE で確率0に当たった
                }
                prob = (double)nijk / (double)nij;
            }
            LL += log(prob);
        }
    }

    if (out_avg_per_sample) *out_avg_per_sample = LL / (double)N;
    if (out_avg_per_var)    *out_avg_per_var    = LL / (double)(N * D);
    if (out_zero_hits)      *out_zero_hits      = zero_hits;
    return LL;
};

struct DeltaStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double stdev = std::numeric_limits<double>::quiet_NaN();
    long long finite_count = 0;  // 有限値として集計できたサンプル数
};

// 親集合に対する j_index（各サンプルの親配置インデックス）を構築
static std::vector<int> buildJIndexForParents(const Dataset& ds,
                                         const std::vector<int>& parents)
{
    std::vector<int> jidx(ds.N, 0);
    if (parents.empty()) return jidx;
    // 右端の親が最下位桁になる混合基数
    std::vector<int> radix;
    build_mixed_radix(parents, ds, radix);
    const int P = (int)parents.size();
    const int D = ds.D;
    const int N = ds.N;
    const int* ds_x_ptr = ds.X_flat.data();
    const int* pa_ptr = parents.data();
    const int* rdx_ptr = radix.data();
    #pragma acc data present(ds_x_ptr[0:N*D])
    {
        #pragma acc parallel present(ds_x_ptr[0:N*D]) loop copy(pa_ptr[0:P],rdx_ptr[0:P])
        for (int n=0; n<N; ++n) {
            const int j = mixed_radix_index_row(D, ds_x_ptr, n, P, pa_ptr, rdx_ptr);
            jidx[n]=j;
        }
    }
    return jidx;
}


// 親集合 before/after のみが異なるノード v について、
// 新データ ds_new 上の「サンプルごとの logL 差分 (after - before)」の
// 平均・標準偏差を計算して返す。
// alpha_ij>0 なら (n_ijk + a/r_i)/(n_ij + a) を用いた平滑化付き確率、0ならMLE。
// ※ MLEで0確率に当たるサンプルは統計から除外（finite のみ集計）。
DeltaStats perSampleDeltaLogLStats(const Dataset& ds_new,
                                   int v,
                                   const std::vector<int>& parents_before,
                                   const std::vector<int>& parents_after,
                                   double alpha_ij);

// ============================================================
// 各エッジ (u->v) を1本ずつ削除したときのスコア変化を計算
// - データセット ds_new : 新しい評価データ
// - g_base : 学習済み構造（init_edges.tsv）
// - counts : all_counts.tsv から読み込んだ n_ijk / n_ij
// - alpha_ij : 平滑化 (MLE時=0.0)
// - ess : BDeu スコア用の等価サンプルサイズ
// - outfile : 結果TSVファイル (u, v, ΔlogL, ΔBIC, ΔK2, ΔBDeu)
//
// 計算の流れ：
//   1. 元スコアを計算（logL, BIC, K2, BDeu）
//   2. 各エッジ u->v を削除した DAG を作成
//   3. そのノード v のみ再スコア（親集合からuを除く）
//   4. 各スコア差 Δ = (after - before) を出力
// ============================================================
void computeEdgeImportanceScores(const Dataset& ds_new,
                                 const DAG& g_base,
                                 const std::vector<Counts>& counts,
                                 double alpha_ij,
                                 double ess,
                                 const std::string& outfile);

// ================================================================
// ブートストラップ学習を B 回まわして、エッジ出現回数を TSV で保存
//
// - 入力: 元データ ds（ヘッダ/基数付き）, スコア種別, ESS(BDeu用),
//         初期構造ファイルパス（空文字なら空DAGで開始OK）,
//         既存HCのパラメタ（tabu/iters/max-parents/max-children/topK 等）,
//         reach モード, jindex キャッシュ, ブートストラップ回数/seed,
//         MI の sample/budget 制御（各バッチで再計算）
//
// - 出力: save_path に TSV で保存（u v count prob）
//         include_zero==true のとき、全 (u!=v) の 0 含む一覧を出力（大規模Dでは注意）
//
// 注意: メモリ節約のためカウントは unordered_map で疎管理
// ================================================================
void runBootstrapStructureCounts(const Dataset& ds,
                                 ScoreType score_type,
                                 double ess_for_bdeu,
                                 const std::string& init_path,
                                 int tabu_tenure,
                                 int iters,
                                 int max_parents,
                                 int max_children,
                                 CandMetric cand_metric,
                                 int topK,
                                 int mi_sample,    // 0: 全行
                                 int mi_budget,    // 0: 全変数
                                 double mi_threshold,
                                 double chi2_p_threshold,
                                 Reachability::Mode reach_mode,
                                 int jindex_cache_cap,
                                 int B,            // ブートストラップ反復回数
                                 uint64_t seed,
                                 const std::string& save_path);


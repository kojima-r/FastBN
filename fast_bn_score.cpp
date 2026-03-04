#include "fast_bn_dataset.hpp"
#include "fast_bn_dag.hpp"
#include "fast_bn_lib.hpp"
#include "fast_bn_score.hpp"

double AssocCandidates::mutual_info_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows) {
    int ru = ds.r[u], rv = ds.r[v];
    std::vector<long long> cu(ru,0), cv(rv,0);
    std::unordered_map<long long, long long> cuv; cuv.reserve(rows.size()*2);
    auto key = [&](int a,int b)->long long { return (long long)a * (long long)rv + b; };

    for (int idx: rows){
        int a = ds.x(idx,u), b = ds.x(idx,v);
        ++cu[a]; ++cv[b]; ++cuv[key(a,b)];
    }
    double N = (double)rows.size();
    double mi=0.0;
    for (auto &kv : cuv) {
        long long ab = kv.second;
        int a = (int)(kv.first / rv);
        int b = (int)(kv.first % rv);
        if (ab==0) continue;
        double pab = ab / N;
        double pa  = cu[a] / N;
        double pb  = cv[b] / N;
        mi += pab * std::log( (pab + 1e-300) / (pa * pb + 1e-300) );
    }
    return mi; // nats
}

std::pair<double,double> AssocCandidates::chi2_p_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows) {
    int ru = ds.r[u], rv = ds.r[v];
    std::vector<long long> cu(ru,0), cv(rv,0);
    std::vector<long long> tab((size_t)ru * rv, 0);

    for (int idx: rows){
        int a = ds.x(idx,u), b = ds.x(idx,v);
        ++cu[a]; ++cv[b];
        ++tab[(size_t)a * rv + b];
    }

    double N = (double)rows.size();
    if (N <= 0.0) return {1.0, 0.0};

    // カイ二乗統計量
    double chi2 = 0.0;
    for (int a=0; a<ru; ++a) {
        for (int b=0; b<rv; ++b) {
            double obs = (double)tab[(size_t)a * rv + b];
            double expv = (double)cu[a] * (double)cv[b] / N;
            if (expv > 0.0) {
                double diff = obs - expv;
                chi2 += diff * diff / expv;
            }
        }
    }
    int df = (ru - 1) * (rv - 1);
    if (df <= 0) return {1.0, 0.0};
    double p = chisq_p_upper(chi2, df);
    return {p, chi2};
}

std::vector<std::vector<int>> AssocCandidates::compute(
    const Dataset& ds,
    int K,
    int budget,
    const std::vector<int>& rows,
    std::mt19937_64& rng,
    CandMetric metric,
    double mi_threshold_nats,
    double chi2_p_threshold)
{
    int D = ds.D;
    std::vector<std::vector<int>> out(D);

    for (int v=0; v<D; ++v) {
        // 候補母集団（budget で間引き）
        std::vector<int> cand;
        cand.reserve(budget>0? budget : (D-1));
        if (budget>0 && budget < D-1) {
            std::unordered_set<int> used; used.reserve(budget*2);
            while ((int)cand.size()<budget) {
                int u = (int)(rng()%D);
                if (u==v) continue;
                if (used.insert(u).second) cand.push_back(u);
            }
        } else {
            for (int u=0; u<D; ++u) if (u!=v) cand.push_back(u);
        }

        if (metric == CandMetric::MI) {
            std::vector<std::pair<double,int>> scored; // (mi,u)
            scored.reserve(cand.size());
            for (int u: cand) {
                double mi = mutual_info_pair(ds, u, v, rows);
                if (mi >= mi_threshold_nats) scored.emplace_back(mi, u);
            }
            if (K>0 && (int)scored.size()>K) {
                std::nth_element(scored.begin(), scored.begin()+K, scored.end(),
                    [](auto& a, auto& b){ return a.first > b.first; });
                scored.resize(K);
            }
            std::vector<int> keep; keep.reserve(scored.size());
            for (auto& p: scored) keep.push_back(p.second);
            std::sort(keep.begin(), keep.end());
            out[v] = std::move(keep);
        } else { // CHI2
            std::vector<std::pair<double,int>> scored; // (score,u), score=-log(p)
            scored.reserve(cand.size());
            for (int u: cand) {
                auto [p, chi2] = chi2_p_pair(ds, u, v, rows);
                if (p <= chi2_p_threshold) {
                    double score = -std::log(std::max(p, 1e-300)); // 小pほど大
                    scored.emplace_back(score, u);
                }
            }
            if (K>0 && (int)scored.size()>K) {
                std::nth_element(scored.begin(), scored.begin()+K, scored.end(),
                    [](auto& a, auto& b){ return a.first > b.first; });
                scored.resize(K);
            }
            std::vector<int> keep; keep.reserve(scored.size());
            for (auto& p: scored) keep.push_back(p.second);
            std::sort(keep.begin(), keep.end());
            out[v] = std::move(keep);
        }
    }
    return out;
}

double MICandidates::mutual_info_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows){
    int ru = ds.r[u], rv = ds.r[v];
    std::vector<long long> cu(ru,0), cv(rv,0);
    std::unordered_map<long long, long long> cuv; cuv.reserve(rows.size()*2);
    auto key = [&](int a,int b)->long long { return (long long)a * (long long)rv + b; };

    for (int idx: rows){
        int au = ds.x(idx,u), bv = ds.x(idx,v);
        ++cu[au]; ++cv[bv]; ++cuv[key(au,bv)];
    }

    double N = (double)rows.size();
    double mi=0.0;
    for (auto &kv : cuv){
        long long ab = kv.second;
        int a = (int)(kv.first / rv);
        int b = (int)(kv.first % rv);
        if (ab==0) continue;
        double pab = ab / N;
        double pa  = cu[a] / N;
        double pb  = cv[b] / N;
        // 数値安定用に 1e-300 を加える
        mi += pab * log( (pab + 1e-300) / (pa * pb + 1e-300) );
    }
    return mi;
}

std::vector<std::vector<int>> MICandidates::compute(const Dataset& ds, int K, int budget, const std::vector<int>& rows, std::mt19937_64& rng, double mi_threshold_nats){
    int D = ds.D;
    std::vector<std::vector<int>> topK(D);

    for (int v=0; v<D; ++v){
        // 変数サンプリング：自身以外から budget 個
        std::vector<int> cand;
        if (budget>0 && budget < D-1) {
            cand.reserve(budget);
            std::unordered_set<int> used; used.reserve(budget*2);
            while ((int)cand.size()<budget){
                int u = (int)(rng()%D);
                if (u==v) continue;
                if (used.insert(u).second) cand.push_back(u);
            }
        } else {
            cand.reserve(D-1);
            for (int u=0;u<D;++u) if (u!=v) cand.push_back(u);
        }

        // (u, MI(u;v)) を計算
        std::vector<std::pair<double,int>> scored;
        scored.reserve(cand.size());
        for (int u: cand){
            double mi = mutual_info_pair(ds, u, v, rows); // nats
            if (mi >= mi_threshold_nats) {                 // しきい値でまず足切り
                scored.emplace_back(mi, u);
            }
        }

        // K>0 のときは上位 K に丸める（MI 降順）
        if (K>0 && (int)scored.size() > K){
            nth_element(scored.begin(), scored.begin()+K, scored.end(),
                        [](const auto& a, const auto& b){ return a.first > b.first; });
            scored.resize(K);
        }

        // u のみ抽出し、binary_search 用にソート（昇順）
        std::vector<int> keep; keep.reserve(scored.size());
        for (auto &p: scored) keep.push_back(p.second);
        std::sort(keep.begin(), keep.end());
        topK[v] = move(keep);
    }
    return topK;
}

double HillClimber::deltaAdd_andBuildNewCounts(int u, int v, Counts& newC) {
    const Counts& curC = nodeCounts[v];
    const int N = ds.N;
    const int ru = ds.r[u];
    const int r_i = curC.r_i;

    // 親なしの場合は j_index を呼ばない
    if (curC.q_i == 0) {
        const int q = 1;
        const int qp = q * ru; // = ru

        newC.q_i = qp;
        newC.r_i = r_i;
        newC.n_ij.assign(qp, 0);
        newC.n_ijk.assign((size_t)qp * r_i, 0);

        auto * __restrict nij  = newC.n_ij.data();
        auto * __restrict nijk = newC.n_ijk.data();

        // 各サンプルについて j' と k を1回で決めてカウント
        for (int n = 0; n < N; ++n) {
            const int xu = ds.x(n, u);   // フラット配列アクセス
            const int k  = ds.x(n, v);
            const int j2 = xu;           // j=0 なので j2 = xu

            ++nij[j2];
            ++nijk[(size_t)j2 * r_i + k];
        }

        double after = scorer.nodeScore(newC);
        return after - nodeScoreNow[v];
    }

    // 親ありの場合
    const std::vector<int>& jindex = jcache.get(v); // Pa(v) に対応する j_index（キャッシュから取得）
    const int q   = curC.q_i;
    const int qp  = q * ru;

    newC.q_i = qp;
    newC.r_i = r_i;
    newC.n_ij.assign(qp, 0);
    newC.n_ijk.assign((size_t)qp * r_i, 0);

    auto * __restrict nij  = newC.n_ij.data();
    auto * __restrict nijk = newC.n_ijk.data();
    const int * __restrict jptr = jindex.data();

    // 各サンプルについて j' と k を1回で決めてカウント
    for (int n = 0; n < N; ++n) {
        const int j  = jptr[n];
        const int xu = ds.x(n, u);   // フラット配列アクセス
        const int k  = ds.x(n, v);
        const int j2 = j * ru + xu;

        ++nij[j2];
        ++nijk[(size_t)j2 * r_i + k];
    }

    double after = scorer.nodeScore(newC);
    return after - nodeScoreNow[v];
}

DeltaStats perSampleDeltaLogLStats(const Dataset& ds_new,
                                   int v,
                                   const std::vector<int>& parents_before,
                                   const std::vector<int>& parents_after,
                                   double alpha_ij)
{
    // before/after の counts と j_index を用意
    Counts C_before = computeCountsForNode_full(v, parents_before, ds_new);
    Counts C_after  = computeCountsForNode_full(v, parents_after , ds_new);

    auto buildJ = [&](const std::vector<int>& pa){
        return buildJIndexForParents(ds_new, pa);
    };
    std::vector<int> j_before = buildJ(parents_before);
    std::vector<int> j_after  = buildJ(parents_after);

    auto prob_from_counts = [&](const Counts& C, int j, int k)->double{
        if (j<0 || j>=C.q_i || k<0 || k>=C.r_i) return 0.0;
        double nij  = (double)C.n_ij[j];
        double nijk = (double)C.n_ijk[(size_t)j*C.r_i + k];
        if (alpha_ij > 0.0) {
            return (nijk + alpha_ij/(double)C.r_i) / (nij + alpha_ij);
        } else {
            if (nij<=0.0 || nijk<=0.0) return 0.0;
            return nijk / nij;
        }
    };

    long long cnt = 0;
    long double sum = 0.0L, sum2 = 0.0L;

    for (int n=0; n<ds_new.N; ++n) {
        int kx = ds_new.x(n,v);
        double p_b = prob_from_counts(C_before, j_before[n], kx);
        double p_a = prob_from_counts(C_after , j_after [n], kx);
        double lb = (p_b>0.0 ? std::log(p_b) : -INFINITY);
        double la = (p_a>0.0 ? std::log(p_a) : -INFINITY);
        double d  = la - lb; // after - before
        if (std::isfinite(d)) {
            ++cnt;
            long double dl = (long double)d;
            sum  += dl;
            sum2 += dl*dl;
        }
    }

    DeltaStats st;
    st.finite_count = cnt;
    if (cnt > 0) {
        long double mu = sum / (long double)cnt;
        st.mean = (double)mu;
        if (cnt > 1) {
            long double var = (sum2 - (long double)cnt*mu*mu) / (long double)(cnt-1);
            st.stdev = (var>0 ? std::sqrt((double)var) : 0.0);
        } else {
            st.stdev = 0.0;
        }
    }
    return st;
}

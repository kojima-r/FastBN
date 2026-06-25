#include "fast_bn_dataset.hpp"
#include "fast_bn_dag.hpp"
#include "fast_bn_lib.hpp"
#include "fast_bn_score.hpp"
#include "fast_bn_utils.hpp"

#include <iostream>

double AssocCandidates::mutual_info_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows) {
    int ru = ds.r[u], rv = ds.r[v];
    size_t row_size = rows.size();
    int ds_N = ds.N;
    int ds_D = ds.D;
    std::vector<long long> cu(ru,0), cv(rv,0);
    std::vector<long long> keys(row_size,0);
    std::vector<long long> unique(row_size,0);
    std::vector<long long> counts(row_size,0);
    const int* _rows = rows.data();
    const int* ds_x_ptr = ds.X_flat.data();
    long long* _cu = cu.data();
    long long* _cv = cv.data();
    long long* _keys = keys.data();
    long long* _unique = unique.data();
    long long* _counts = counts.data();
    double mi = 0.0;

    #pragma acc data copyin(_rows[0:row_size]) \
                     copy(_cu[0:ru], _cv[0:rv], mi) \
                     create(_keys[0:row_size], _unique[0:row_size], _counts[0:row_size]) \
                     present(ds_x_ptr[0:ds_N*ds_D])
    {
    #pragma acc parallel loop independent
    for (int i=0; i< row_size; i++){
        int idx = _rows[i];
        int a = ds_x_ptr[(size_t)idx * ds_D + u];
        int b = ds_x_ptr[(size_t)idx * ds_D + v];
        #pragma acc atomic update
        ++_cu[a];
        #pragma acc atomic update
        ++_cv[b];
        _keys[i] = (long long)a*rv + b;
    }

    long long* d_keys  = nullptr;
    long long* d_unique_indices  = nullptr;
    long long* d_counts_indices = nullptr;
    #pragma acc host_data use_device(_keys,_unique,_counts)
    {
        d_keys = _keys;
        d_unique_indices  = _unique;
        d_counts_indices = _counts;
    }
    size_t num_unique = count_frequencies(d_keys, row_size, d_unique_indices, d_counts_indices);
    double N = (double)row_size;
    #pragma acc parallel loop reduction(+:mi)
    for(size_t i = 0; i < num_unique; i++ ){
        long long ab = _counts[i];
        int a = (int)(_unique[i] / rv);
        int b = (int)(_unique[i] % rv);
        if (ab==0) continue;
        double pab = ab / N;
        double pa  = _cu[a] / N;
        double pb  = _cv[b] / N;
        mi += pab * std::log( (pab + 1e-300) / (pa * pb + 1e-300) );
    }
    }
    return mi; // nats
}

std::pair<double,double> AssocCandidates::chi2_p_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows) {
    int ru = ds.r[u], rv = ds.r[v];
    size_t row_size = rows.size();
    if (row_size <= 0) return {1.0, 0.0};
    int df = (ru - 1) * (rv - 1);
    if (df <= 0) return {1.0, 0.0};
    int ds_N = ds.N;
    int ds_D = ds.D;
    std::vector<long long> cu(ru,0), cv(rv,0);
    std::vector<long long> tab((size_t)ru * rv, 0);
    const int* _rows = rows.data();
    const int* ds_x_ptr = ds.X_flat.data();
    long long* _cu = cu.data();
    long long* _cv = cv.data();
    long long* _tab = tab.data();
    // カイ二乗統計量
    double chi2 = 0.0;

    #pragma acc data copyin(_rows[0:row_size]) \
                     copy(_cu[0:ru], _cv[0:rv], _tab[0:ru*rv], chi2) \
                     present(ds_x_ptr[0:ds_N*ds_D])
    {
    #pragma acc parallel loop independent
    for (int i=0; i< row_size; i++){
        int idx = _rows[i];
        int a = ds_x_ptr[(size_t)idx * ds_D + u];
        int b = ds_x_ptr[(size_t)idx * ds_D + v];
        #pragma acc atomic update
        ++_cu[a];
        #pragma acc atomic update
        ++_cv[b];
        #pragma acc atomic update
        ++_tab[(size_t)a * rv + b];
    }

    double N = (double)row_size;

    #pragma acc parallel loop collapse(2) reduction(+:chi2)
    for (int a=0; a<ru; ++a) {
        for (int b=0; b<rv; ++b) {
            double obs = (double)_tab[(size_t)a * rv + b];
            double expv = (double)_cu[a] * (double)_cv[b] / N;
            if (expv > 0.0) {
                double diff = obs - expv;
                chi2 += diff * diff / expv;
            }
        }
    }
    }
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
    size_t row_size = rows.size();
    int ds_N = ds.N;
    int ds_D = ds.D;
    std::vector<long long> cu(ru,0), cv(rv,0);
    std::vector<long long> keys(row_size,0);
    std::vector<long long> unique(row_size,0);
    std::vector<long long> counts(row_size,0);
    const int* _rows = rows.data();
    const int* ds_x_ptr = ds.X_flat.data();
    long long* _cu = cu.data();
    long long* _cv = cv.data();
    long long* _keys = keys.data();
    long long* _unique = unique.data();
    long long* _counts = counts.data();
    double mi = 0.0;

    #pragma acc data copyin(_rows[0:row_size]) \
                     copy(_cu[0:ru], _cv[0:rv], mi) \
                     create(_keys[0:row_size], _unique[0:row_size], _counts[0:row_size]) \
                     present(ds_x_ptr[0:ds_N*ds_D])
    {
    #pragma acc parallel loop independent
    for (int i=0; i< row_size; i++){
        int idx = _rows[i];
        int au = ds_x_ptr[(size_t)idx * ds_D + u];
        int bv = ds_x_ptr[(size_t)idx * ds_D + v];
        #pragma acc atomic update
        ++_cu[au];
        #pragma acc atomic update
        ++_cv[bv];
        _keys[i] = (long long)au*rv + bv;
    }
    long long* d_keys  = nullptr;
    long long* d_unique_indices  = nullptr;
    long long* d_counts_indices = nullptr;
    #pragma acc host_data use_device(_keys,_unique,_counts)
    {
        d_keys = _keys;
        d_unique_indices  = _unique;
        d_counts_indices = _counts;
    }
    size_t num_unique = count_frequencies(d_keys, row_size, d_unique_indices, d_counts_indices);
    double N = (double)row_size;
    #pragma acc parallel loop reduction(+:mi)
    for(size_t i = 0; i < num_unique; i++ ){
        long long ab = _counts[i];
        int a = (int)(_unique[i] / rv);
        int b = (int)(_unique[i] % rv);
        if (ab==0) continue;
        double pab = ab / N;
        double pa  = _cu[a] / N;
        double pb  = _cv[b] / N;
        mi += pab * std::log( (pab + 1e-300) / (pa * pb + 1e-300) );
    }
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
    const int D = ds.D;
    const int ru = ds.r[u];
    const int r_i = curC.r_i;
    const int q = (curC.q_i == 0) ? 1 : curC.q_i;
    const int qp = q * ru;

    newC.assign(qp,r_i);
//    newC.acc_update_host();

    const int* __restrict ds_x_ptr = ds.X_flat.data();

    long long* __restrict nij  = newC.n_ij.data();
    long long* __restrict nijk = newC.n_ijk.data();

    long long* __restrict keys_ij = this->work_keys_ij.data();
    long long* __restrict keys_ijk = this->work_keys_ijk.data();
    long long* __restrict _unique_indices = this->unique_indices.data();
    long long* __restrict _counts_indices = this->counts_indices.data();

    #pragma acc data present(keys_ij[0:N],keys_ijk[0:N],ds_x_ptr[0:N*D],nij[0:qp],nijk[0:qp*r_i])
    {
        if (curC.q_i == 0)
        {
            // 親なしの場合は j_index を呼ばない
            // 各サンプルについて j' と k を1回で決めてカウント
            #pragma acc parallel loop independent
            for (int n = 0; n < N; ++n) {
                const int xu = ds_x_ptr[n * D + u];
                const int k  = ds_x_ptr[n * D + v];
                const int j2 = xu;           // j=0 なので j2 = xu
                keys_ij[n] = (long long)j2;
                keys_ijk[n] = (long long)j2 * r_i + k;
            }
        }
        else{
            // 親ありの場合
            const int* __restrict jptr = jcache.get(v).data(); // Pa(v) に対応する j_index（キャッシュから取得）
            // 各サンプルについて j' と k を1回で決めてカウント
            #pragma acc parallel loop independent present(jptr[0:N])
            for (int n = 0; n < N; ++n) {
                const int j  = jptr[n];
                const int xu = ds_x_ptr[n * D + u];
                const int k  = ds_x_ptr[n * D + v];
                const int j2 = j * ru + xu;
                keys_ij[n] = (long long)j2;
                keys_ijk[n] = (long long)j2 * r_i + k;
            }
        }

        // OpenACCが管理しているkeys_ptrの「GPU側の生ポインタ」を取得して渡す
        long long* d_keys_ij  = nullptr;
        long long* d_keys_ijk = nullptr;
        long long* d_unique_indices  = nullptr;
        long long* d_counts_indices = nullptr;
        #pragma acc host_data use_device(keys_ij,keys_ijk,_unique_indices,_counts_indices)
        {
            d_keys_ij = keys_ij;
            d_keys_ijk = keys_ijk;
            d_unique_indices  = _unique_indices;
            d_counts_indices = _counts_indices;
        }
        // 頻度分布から nij nijk に代入する
        size_t num_unique_ij = count_frequencies(d_keys_ij, N, d_unique_indices, d_counts_indices);
        #pragma acc data present(nij[0:qp],_unique_indices[0:N],_counts_indices[0:N])
        {
            #pragma acc parallel loop independent
            for (int j=0; j<qp; ++j){
                nij[j] = 0;
            }
            #pragma acc parallel loop independent
            for (size_t i = 0; i < num_unique_ij; ++i) {
                nij[_unique_indices[i]] = _counts_indices[i];
            }
        }
        size_t num_unique_ijk = count_frequencies(d_keys_ijk, N, d_unique_indices, d_counts_indices);
        #pragma acc data present(nijk[0:qp*r_i],_unique_indices[0:N],_counts_indices[0:N])
        {
            #pragma acc parallel loop independent
            for (int j=0; j<qp*r_i; ++j){
                nijk[j] = 0;
            }
            #pragma acc parallel loop independent
            for (size_t i = 0; i < num_unique_ijk; ++i) {
                nijk[_unique_indices[i]] = _counts_indices[i];
            }
        }
    }
#ifndef RELAX
    newC.acc_update_host();
#endif
    return scorer.nodeScore(newC) - nodeScoreNow[v];
}

#pragma acc routine seq
int calc_right(const int* ds_r_ptr,int u,const int* _pa,int pa_size)
{
    int pos = -1;
    for (int t=0; t<pa_size; ++t){
        if(_pa[t]==u){
            pos=t;
        }
    }
    if(pos<0) return 0;
    int right = 1;
    for (int t=pos+1; t<pa_size; ++t){
        right *= ds_r_ptr[_pa[t]];
    }
    return right;
}

double HillClimber::deltaRemove_andBuildNewCounts(int u, int v, const std::vector<int> &Pa, Counts& newC)
{
    const Counts& curC = nodeCounts[v]; // Pa(v) に基づく現在の counts

    int right = 1;
    const int* ds_r_ptr = ds.r.data();
    int D = ds.D;
    const int* _pa = Pa.data();
    int pa_size = Pa.size();
    #pragma acc serial present(ds_r_ptr[0:D],_pa[0:pa_size]) copy(right)
    {
        right = calc_right(ds_r_ptr,u,_pa,pa_size);
    }
    if(right==0) throw std::runtime_error("deltaRemove: u is not a parent of v (logic error).");

    int r_u = ds.r[u];
    int r_i = curC.r_i;
    int q   = curC.q_i;
    int q2  = q / r_u; // u の桁を落とすと親配置数は r_u 倍減る

    int period = right * r_u; // 1つ上の繰り返し周期

    newC.assign(q2,r_i);

    long long int* __restrict new_nij  = newC.n_ij.data();
    long long int* __restrict new_nijk = newC.n_ijk.data();
    const long long int* __restrict cur_nij  = curC.n_ij.data();
    const long long int* __restrict cur_nijk = curC.n_ijk.data();

    // 旧インデックス j を新インデックス j' へ写像して合算
    // j' = floor(j / period) * right + (j % right)
    #pragma acc data present(new_nij[0:q2],new_nijk[0:q2*r_i],cur_nij[0:q],cur_nijk[0:q*r_i])
    {
//    #pragma acc kernels present(new_nij[0:q2],new_nijk[0:q2*r_i],cur_nij[0:q],cur_nijk[0:q*r_i])
    #pragma acc serial present(new_nij[0:q2],new_nijk[0:q2*r_i],cur_nij[0:q],cur_nijk[0:q*r_i])
    for (int j=0; j<q; ++j){
        int jp = (j / period) * right + (j % right);
        new_nij[jp] += cur_nij[j];
        size_t base_old = (size_t)j * r_i;
        size_t base_new = (size_t)jp * r_i;
        for (int k=0;k<r_i;++k){
            new_nijk[base_new + k] += cur_nijk[base_old + k];
        }
    }
    }
#ifndef RELAX
    newC.acc_update_host();
#endif
    return scorer.nodeScore(newC) - nodeScoreNow[v];
}

DeltaStats perSampleDeltaLogLStats(const Dataset& ds_new,
                                   int v,
                                   const std::vector<int>& parents_before,
                                   const std::vector<int>& parents_after,
                                   double alpha_ij)
{
    // before/after の counts と j_index を用意
    std::vector<int> radix_before;
    build_mixed_radix(parents_before, ds_new.r, radix_before);
    Counts C_before = computeCountsForNode_full(v, parents_before, ds_new, radix_before);

    std::vector<int> radix_after;
    build_mixed_radix(parents_after, ds_new.r, radix_after);
    Counts C_after  = computeCountsForNode_full(v, parents_after , ds_new, radix_after);

    std::vector<int> j_before = buildJIndexForParents(ds_new,parents_before);
    std::vector<int> j_after  = buildJIndexForParents(ds_new,parents_after);

    const int N = ds_new.N;
    const int D = ds_new.D;
    const int* ds_new_x_ptr = ds_new.X_flat.data();

    const long long* nb_ij = C_before.n_ij.data();
    const long long* nb_ijk = C_before.n_ijk.data();
    const long long* na_ij = C_after.n_ij.data();
    const long long* na_ijk = C_after.n_ijk.data();
    const int* jb_ptr = j_before.data();
    const int* ja_ptr = j_after.data();

    long long cnt = 0;
    long double sum = 0.0L, sum2 = 0.0L;

    #pragma acc parallel loop reduction(+:sum, sum2, cnt) \
        copyin(jb_ptr[0:N], ja_ptr[0:N]) \
        present(nb_ij[0:C_before.q_i], nb_ijk[0:(size_t)C_before.q_i * C_before.r_i], \
                na_ij[0:C_after.q_i],  na_ijk[0:(size_t)C_after.q_i * C_after.r_i], \
                ds_new_x_ptr[0:N*D])
    for (int n=0; n<N; ++n) {
        int kx = ds_new_x_ptr[(size_t)n * D + v];

        double p_b = 0.0;
        int jb = jb_ptr[n];
        if (jb>=0 && jb<C_before.q_i && kx<0 && kx>=C_before.r_i){
            double nij  = (double)nb_ij[jb];
            double nijk = (double)nb_ijk[(size_t)jb*C_before.r_i + kx];
            if (alpha_ij > 0.0) {
                p_b = (nijk + alpha_ij/(double)C_before.r_i) / (nij + alpha_ij);
            } else if (nij>0.0 && nijk>0.0){
                p_b = nijk / nij;
            }
        }

        double p_a = 0.0;
        int ja = ja_ptr[n];
        if (ja>=0 && ja<C_after.q_i && kx<0 && kx>=C_after.r_i){
            double nij  = (double)na_ij[ja];
            double nijk = (double)na_ijk[(size_t)ja*C_after.r_i + kx];
            if (alpha_ij > 0.0) {
                p_a = (nijk + alpha_ij/(double)C_after.r_i) / (nij + alpha_ij);
            } else if (nij>0.0 && nijk>0.0){
                p_a = nijk / nij;
            }
        }

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

void computeEdgeImportanceScores(const Dataset& ds_new,
                                 const DAG& g_base,
                                 const std::vector<Counts>& counts,
                                 double alpha_ij,
                                 double ess,
                                 const std::string& outfile)
{
    std::ofstream fout(outfile);
    if (!fout) throw std::runtime_error("Failed to open edge-importance file: " + outfile);
    fout << "u\tv\tΔlogL\tΔBIC\tΔK2\tΔBDeu\n";

    const int D = ds_new.D;
    const int N = ds_new.N;
    const double logN = std::log((double)std::max(1, N));

#ifdef RELAX
    // ========== スコア計算関数群 ==========
    auto nodeLogLikelihood = [&](const Counts &C) -> double {
        double ll = 0.0;
        const long long* _nij  = C.n_ij.data();
        const long long* _nijk = C.n_ijk.data();
        int q_i = C.q_i;
        int r_i = C.r_i;
        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
            #pragma acc parallel loop reduction(+:ll)
            for (int j = 0; j < q_i; ++j) {
                double local_ll = 0.0;
                for (int k = 0; k < r_i; ++k) {
                    double nij = (double)_nij[j];
                    if (nij > 0){
                        long long nijk = _nijk[j*r_i + k];
                        if (nijk > 0){
                            local_ll += nijk * (std::log((double)nijk) - std::log(nij));
                        }
                    }
                }
                ll += local_ll;
            }
        }
        return ll;
    };

    auto nodeBIC = [&](const Counts &C) -> double {
        double ll = 0.0;
        const long long* _nij  = C.n_ij.data();
        const long long* _nijk = C.n_ijk.data();
        int q_i = C.q_i;
        int r_i = C.r_i;
        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
            #pragma acc parallel loop reduction(+:ll)
            for (int j = 0; j < q_i; ++j) {
                double local_ll = 0.0;
                for (int k = 0; k < r_i; ++k) {
                    double nij = (double)_nij[j];
                    if (nij > 0){
                        long long nijk = _nijk[j*r_i + k];
                        if (nijk > 0){
                            local_ll += nijk * (std::log((double)nijk) - std::log(nij));
                        }
                    }
                }
                ll += local_ll;
            }
        }
        int d = (C.r_i - 1) * C.q_i;
        double pen = 0.5 * d * logN;
        return ll - pen;
    };

    auto nodeK2 = [&](const Counts &C) -> double {
        double s = 0.0;
        const long long* _nij  = C.n_ij.data();
        const long long* _nijk = C.n_ijk.data();
        int q_i = C.q_i;
        int r_i = C.r_i;
        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
            #pragma acc parallel loop reduction(+:s)
            for (int j = 0; j < q_i; ++j) {
                double nij = (double)_nij[j];
                s += std::lgamma((double)r_i) - std::lgamma(nij + (double)r_i);
                for (int k = 0; k < r_i; ++k) {
                    double nijk = (double)_nijk[j*r_i + k];
                    s += std::lgamma(nijk + 1.0);
                }
            }
        }
        return s;
    };

    auto nodeBDeu = [&](const Counts &C) -> double {
        double s = 0.0;
        if (C.q_i == 0) return -INFINITY;
        double alpha_ij_local = ess / (double)C.q_i;
        double alpha_ijk_base = alpha_ij_local / (double)C.r_i;
        const long long* _nij  = C.n_ij.data();
        const long long* _nijk = C.n_ijk.data();
        int q_i = C.q_i;
        int r_i = C.r_i;
        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
            #pragma acc parallel loop reduction(+:s) present(_nij[0:q_i],_nijk[0:q_i*r_i])
            for (int j = 0; j < q_i; ++j) {
                double nij = (double)_nij[j];
                s += std::lgamma(alpha_ij_local) - std::lgamma(nij + alpha_ij_local);
                double local_s = 0.0;
                for (int k = 0; k < r_i; ++k) {
                    double nijk = (double)_nijk[j*C.r_i + k];
                    local_s += std::lgamma(nijk + alpha_ijk_base) - std::lgamma(alpha_ijk_base);
                }
                s += local_s;
            }
        }
        return s;
    };
#else
    // ========== スコア計算関数群 ==========
    auto nodeLogLikelihood = [&](const Counts &C) -> double {
        double ll = 0.0;
        const long long* _nij  = C.n_ij.data();
        const long long* _nijk = C.n_ijk.data();
        int q_i = C.q_i;
        int r_i = C.r_i;
//        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
//            #pragma acc parallel loop collapse(2) reduction(+:ll)
            for (int j = 0; j < q_i; ++j) {
                for (int k = 0; k < r_i; ++k) {
                    double nij = (double)_nij[j];
                    if (nij > 0){
                        long long nijk = _nijk[j*r_i + k];
                        if (nijk > 0){
                            ll += nijk * (std::log((double)nijk) - std::log(nij));
                        }
                    }
                }
            }
        }
        return ll;
    };

    auto nodeBIC = [&](const Counts &C) -> double {
        double ll = 0.0;
        const long long* _nij  = C.n_ij.data();
        const long long* _nijk = C.n_ijk.data();
        int q_i = C.q_i;
        int r_i = C.r_i;
//        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
//            #pragma acc parallel loop collapse(2) reduction(+:ll)
            for (int j = 0; j < q_i; ++j) {
                for (int k = 0; k < r_i; ++k) {
                    double nij = (double)_nij[j];
                    if (nij > 0){
                        long long nijk = _nijk[j*r_i + k];
                        if (nijk > 0){
                            ll += nijk * (std::log((double)nijk) - std::log(nij));
                        }
                    }
                }
            }
        }
        int d = (C.r_i - 1) * C.q_i;
        double pen = 0.5 * d * logN;
        return ll - pen;
    };

    auto nodeK2 = [&](const Counts &C) -> double {
        double s = 0.0;
        const long long* _nij  = C.n_ij.data();
        const long long* _nijk = C.n_ijk.data();
        int q_i = C.q_i;
        int r_i = C.r_i;
//        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
//            #pragma acc parallel loop reduction(+:s)
            for (int j = 0; j < q_i; ++j) {
                double nij = (double)_nij[j];
                s += std::lgamma((double)r_i) - std::lgamma(nij + (double)r_i);
                for (int k = 0; k < r_i; ++k) {
                    double nijk = (double)_nijk[j*r_i + k];
                    s += std::lgamma(nijk + 1.0);
                }
            }
        }
        return s;
    };

    auto nodeBDeu = [&](const Counts &C) -> double {
        double s = 0.0;
        if (C.q_i == 0) return -INFINITY;
        double alpha_ij_local = ess / (double)C.q_i;
        double alpha_ijk_base = alpha_ij_local / (double)C.r_i;
        const long long* _nij  = C.n_ij.data();
        const long long* _nijk = C.n_ijk.data();
        int q_i = C.q_i;
        int r_i = C.r_i;
//        #pragma acc data present(_nij[0:q_i],_nijk[0:q_i*r_i])
        {
//            #pragma acc parallel loop reduction(+:s) present(_nij[0:q_i],_nijk[0:q_i*r_i])
            for (int j = 0; j < q_i; ++j) {
                double nij = (double)_nij[j];
                s += std::lgamma(alpha_ij_local) - std::lgamma(nij + alpha_ij_local);
                for (int k = 0; k < r_i; ++k) {
                    double nijk = (double)_nijk[j*C.r_i + k];
                    s += std::lgamma(nijk + alpha_ijk_base) - std::lgamma(alpha_ijk_base);
                }
            }
        }
        return s;
    };
#endif

    // ========== 元スコア計算 ==========
    double baseLL = 0.0, baseBIC = 0.0, baseK2 = 0.0, baseBDeu = 0.0;
    for (int v = 0; v < D; ++v) {
        const auto pa = g_base.parents(v);
        std::vector<int> radix;
        build_mixed_radix(pa, ds_new.r, radix);
        Counts C = computeCountsForNode_full(v, pa, ds_new, radix);

        baseLL  += nodeLogLikelihood(C);
        baseBIC += nodeBIC(C);
        baseK2  += nodeK2(C);
        baseBDeu+= nodeBDeu(C);
/*
        if (C.q_i == 0) baseBDeu = -INFINITY;
        else{
            double alpha_ij_local = ess / (double)C.q_i;
            double alpha_ijk_base = alpha_ij_local / (double)C.r_i;
            for (int j = 0; j < C.q_i; ++j) {
                double nij = (double)C.n_ij[j];
                baseK2 += std::lgamma((double)C.r_i) - std::lgamma(nij + (double)C.r_i);
                for (int k = 0; k < C.r_i; ++k) {
                    double nijk = (double)C.n_ijk[j*C.r_i + k];
                    baseK2 += std::lgamma(nijk + 1.0);
                }
                baseBDeu += std::lgamma(alpha_ij_local) - std::lgamma(nij + alpha_ij_local);
                for (int k = 0; k < C.r_i; ++k) {
                    double nijk = (double)C.n_ijk[j*C.r_i + k];
                    baseBDeu += std::lgamma(nijk + alpha_ijk_base) - std::lgamma(alpha_ijk_base);
                }
                if (nij <= 0) continue;
                for (int k = 0; k < C.r_i; ++k) {
                    long long nijk = C.n_ijk[j*C.r_i + k];
                    if (nijk == 0) continue;
                    baseLL += nijk * (std::log((double)nijk) - std::log(nij));
                    baseBIC += nijk * (std::log((double)nijk) - std::log(nij));
                }
            }
            int d = (C.r_i - 1) * C.q_i;
            double pen = 0.5 * d * logN;
            baseBIC -= pen;
        }
*/
    }

    // ========== 各エッジ削除時のスコア差 ==========
    for (int u = 0; u < D; ++u) {
        for (int v = 0; v < D; ++v) {
            if (!g_base.hasEdge(u, v)) continue;

            DAG g_mod = g_base;
            g_mod.removeEdge(u, v);

            // 対象ノード v のみ再スコア
            const auto pa_mod = g_mod.parents(v);
            std::vector<int> radix_mod;
            build_mixed_radix(pa_mod, ds_new.r, radix_mod);
            Counts C_mod = computeCountsForNode_full(v, pa_mod, ds_new, radix_mod);
            double ll_new   = nodeLogLikelihood(C_mod);
            double bic_new  = nodeBIC(C_mod);
            double k2_new   = nodeK2(C_mod);
            double bdeu_new = nodeBDeu(C_mod);

            // 元スコアとの差分（Δ）
            const auto pa_base = g_base.parents(v);
            std::vector<int> radix_base;
            build_mixed_radix(pa_base, ds_new.r, radix_base);
            Counts C_base = computeCountsForNode_full(v, pa_base, ds_new, radix_base);
            double deltaLL   = ll_new   - nodeLogLikelihood(C_base);
            double deltaBIC  = bic_new  - nodeBIC(C_base);
            double deltaK2   = k2_new   - nodeK2(C_base);
            double deltaBDeu = bdeu_new - nodeBDeu(C_base);

            auto pa_before = g_base.parents(v);
            std::vector<int> pa_after; pa_after.reserve(pa_before.size());
            for (int p : pa_before) if (p!=u) pa_after.push_back(p);
            DeltaStats st = perSampleDeltaLogLStats(ds_new, v, pa_before, pa_after, alpha_ij);

            fout << u << "\t" << v << "\t"
                 << deltaLL   << "\t"
                 << deltaBIC  << "\t"
                 << deltaK2   << "\t"
                 << deltaBDeu << "\t"
                 << std::setprecision(12) << st.mean << '\t'
                 << std::setprecision(12) << st.stdev << '\n';
        }
    }
    fout.close();
    std::cerr << "[info] wrote edge importance results to " << outfile << "\n";
}

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
                                 const std::string& save_path)
{
    if (B <= 0) throw std::runtime_error("bootstrap B must be > 0");
    std::mt19937_64 rng(seed);

    // 疎カウント：エッジ (u->v) の出現回数
    std::unordered_map<uint64_t, uint32_t> edge_counts;
    edge_counts.reserve((size_t)ds.D * 8); // 適当に初期予約（出現が疎な想定）

    for (int b = 0; b < B; ++b) {
        // 1) ブートストラップ標本の生成
        Dataset ds_b = bootstrapResampleDataset(ds, rng);

        const int D = ds_b.D;
        const int N = ds_b.N;
        const int* ds_b_x_ptr = ds_b.X_flat.data();
        const int* ds_b_r_ptr = ds_b.r.data();
        #pragma acc enter data copyin(ds_b_x_ptr[0:N*D],ds_b_r_ptr[0:D])


        // 2) 初期構造のロード
        DAG init = loadInitEdges(ds_b.D, init_path);

        // 3) 候補親Kの前処理（サンプル/バジェットは引数で制御）
        std::vector<std::vector<int>> topKlist;
        if (topK >0 || mi_threshold>0.0 || chi2_p_threshold<1.0){
            std::vector<int> rows;
            if (mi_sample > 0 && mi_sample < ds_b.N) {
                rows.resize(ds_b.N); iota(rows.begin(), rows.end(), 0);
                shuffle(rows.begin(), rows.end(), rng);
                rows.resize(mi_sample);
                std::sort(rows.begin(), rows.end());
            } else {
                rows.resize(ds_b.N); iota(rows.begin(), rows.end(), 0);
            }
            int budget = (mi_budget > 0 ? mi_budget : (ds_b.D - 1));
            //topKlist = MICandidates::compute(ds_b, topK, budget, rows, rng, mi_threshold);
            topKlist = AssocCandidates::compute(
                ds, topK, budget, rows, rng, cand_metric, mi_threshold, chi2_p_threshold
            );
        }
        
        // 4) 探索器の構築 & 実行（学習本体）
        HillClimber hc(ds_b, score_type, ess_for_bdeu, init, reach_mode, jindex_cache_cap);
        hc.max_iter = iters;
        hc.tabu_tenure = tabu_tenure;
        hc.max_parents = max_parents;
        hc.max_children = max_children;
        if (topK > 0) { hc.candParents = std::move(topKlist); hc.topK = topK; }
        hc.verbose = false; // ブートストラップ時は静かに

        auto [g_learned, score, it] = hc.run(/*use_tabu=*/(tabu_tenure > 0));

        // 5) 出現エッジを蓄積
        for (auto &e : g_learned.edges()) {
            uint64_t key = edgeKeyUV((uint32_t)e.first, (uint32_t)e.second);
            auto it = edge_counts.find(key);
            if (it == edge_counts.end()) edge_counts.emplace(key, 1u);
            else ++(it->second);
        }

        // 進捗表示（必要なら）
        if ( (b % 10) == 0 ) {
            std::cerr << "[bootstrap] " << (b+1) << "/" << B
                      << " edges_seen=" << edge_counts.size() << "\n";
        }

        #pragma acc exit data delete(ds_b_x_ptr[0:N*D],ds_b_r_ptr[0:D])
    }

    // 6) 出力：TSV（u, v, count, prob）
    //    include_zero_edges が true なら（u!=v）の全組み合わせで 0 も出力（巨大Dでは注意）
    namespace fs = std::filesystem;
    fs::path base_path(save_path);
    std::string stem = base_path.stem().string();   // 例: "boot_edges"
    std::string ext  = base_path.extension().string(); // 例: ".tsv"
    fs::path parent  = base_path.parent_path();

    std::ostringstream seed_str;
    seed_str << std::setfill('0') << std::setw(4) << (seed % 10000);
    // ファイル名: <stem>_seed<seed><ext>
    std::ostringstream oss;
    oss << stem << "_seed" << seed_str.str() << ext;
    fs::path out_path = parent / oss.str();

    std::ofstream fout(out_path);

    if (!fout) throw std::runtime_error("Failed to open --save-bootstrap-counts: " + save_path);

    fout << "u\tv\tcount\tprob\n";

    // 出現したエッジのみ
    // キーを (u,v) に分解して出力
    for (auto &kv : edge_counts) {
        uint64_t key = kv.first;
        uint32_t c = kv.second;
        uint32_t u = (uint32_t)(key >> 32);
        uint32_t v = (uint32_t)(key & 0xffffffffULL);
        double p = (double)c / (double)B;
        fout << u << "\t" << v << "\t" << c << "\t" << std::setprecision(12) << p << "\n";
    }

    fout.close();
    std::cerr << "[info] wrote bootstrap edge counts to " << save_path
              << " (B=" << B << ", unique_edges=" << edge_counts.size() << ")\n";
}

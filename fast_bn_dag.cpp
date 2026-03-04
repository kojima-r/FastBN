#include "fast_bn_dag.hpp"

// 親集合 Pa(i) に対する counts のフル再集計（O(N)）
// 本実装では、ADD のときは分割のため O(N) が不可避だが、REMOVE は完全インクリメンタル化する。
Counts computeCountsForNode_full(int i, const std::vector<int>& parents, const Dataset& ds){
    const int N   = ds.N;
    const int r_i = ds.r[i];
    int q_i = 1;
    for (int p: parents) q_i *= ds.r[p];
    
    std::vector<long long> nij(q_i,0), nijk((size_t)q_i*r_i,0);

    if (parents.empty()){
        // 親無し: j は常に 0
        for (int n=0;n<ds.N;++n){
            int k=ds.x(n, i);
            ++nijk[k];
        }
        nij.assign(1, (long long)N);
    } else {
        std::vector<int> radix;
        build_mixed_radix(parents, ds.r, radix);
        for (int n=0;n<N;++n){
            //int j = mixedRadixIndex(parents, ds.r, ds.X[n]);
            const int j = mixed_radix_index_row(ds, n, parents, radix);
            const int k = ds.x(n, i);
            ++nijk[(size_t)j*r_i + k];
            ++nij[j];
        }
    }
    return {move(nijk), move(nij), q_i, r_i};
}


#include <iostream>

#ifdef __NVCOMPILER
#include <openacc.h>
#endif

#include "fast_bn_counts.hpp"

// 親集合 Pa(i) に対する counts のフル再集計（O(N)）
// 本実装では、ADD のときは分割のため O(N) が不可避だが、REMOVE は完全インクリメンタル化する。
Counts computeCountsForNode_full(int i, const std::vector<int>& parents, const Dataset& ds, const std::vector<int>& radix){
    const int N   = ds.N;
    const int r_i = ds.r[i];
    int q_i = 1;
    for (int p: parents) q_i *= ds.r[p];
    
    std::vector<long long> nijk_vec((size_t)q_i * r_i, 0);
    std::vector<long long> nij_vec(q_i, 0);

    // ポインタを取り出す（OpenACCが配列範囲を認識しやすくするため）
    long long* nijk = nijk_vec.data();
    long long* nij = nij_vec.data();
    const int* ds_x_ptr = ds.X_flat.data();
    const int D = ds.D;

    const int* rdx_ptr = radix.data();
    const int* pa_ptr = parents.data();
    const int P = (int)parents.size();

    #pragma acc data present(ds_x_ptr[0:N*D])
    {
        if (parents.empty()) {
            // 親なしの場合の並列化
            #pragma acc serial present(ds_x_ptr[0:N*D]) copy(rdx_ptr[0:P], nijk[0:q_i*r_i], nij[0:q_i])
            {
                for (int n = 0; n < N; ++n) {
                    int k = ds_x_ptr[(size_t)n * D + i];
                    ++nijk[k];
                }
                nij[0] = (long long)N;
            }
        } else {
            // 親ありの場合の並列化
            #pragma acc serial present(ds_x_ptr[0:N*D]) copy(pa_ptr[0:P], rdx_ptr[0:P], nijk[0:q_i*r_i], nij[0:q_i])
            {
                for (int n = 0; n < N; ++n) {
                    int j = 0;
                    for (int t = 0; t < P; ++t) {
                        j += ds_x_ptr[(size_t)n * D + pa_ptr[t]] * rdx_ptr[t];
                    }
                    int k = ds_x_ptr[(size_t)n * D + i];
                    ++nijk[(size_t)j * r_i + k];
                    ++nij[j];
                }
            }
        }
    }
    return {move(nijk_vec), move(nij_vec), q_i, r_i};
}


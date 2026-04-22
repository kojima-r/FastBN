#include "fast_bn_dag.hpp"

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

    std::cout << "computeCountsForNode_full " <<
        " D " << D <<
        " N " << N <<
        " P " << P << std::endl;

    if (parents.empty()) {
        // 親なしの場合の並列化
        #pragma acc kernels present(ds_x_ptr[0:N*D]) copy(nijk[0:q_i*r_i], nij[0:q_i])
        {
            for (int n = 0; n < N; ++n) {
                int k = ds_x_ptr[(size_t)n * D + i];
                #pragma acc atomic
                ++nijk[k];
            }
            nij[0] = (long long)N;
        }
    } else {
        // 親ありの場合の並列化
        // data指示文で、データセットをGPUに送り、カウント配列をやり取りする範囲を指定
        #pragma acc data present(ds_x_ptr[0:N*D], rdx_ptr[0:P], pa_ptr[0:P]) copy(nijk[0:q_i*r_i], nij[0:q_i])
        {
            #pragma acc parallel loop
            for (int n = 0; n < N; ++n) {
                int j = 0;
                for (int t = 0; t < P; ++t) {
                    j += ds_x_ptr[(size_t)n * D + pa_ptr[t]] * rdx_ptr[t];
                }
                int k = ds_x_ptr[(size_t)n * D + i];

                // 競合を防ぐためにアトミック処理を指定
                #pragma acc atomic
                ++nijk[(size_t)j * r_i + k];
                #pragma acc atomic
                ++nij[j];
            }
        }
    }
    #pragma acc exit data delete(rdx_ptr[0:P], pa_ptr[0:P])
    return {move(nijk_vec), move(nij_vec), q_i, r_i};
}


#include "fast_bn_utils.hpp"

#ifdef __NVCOMPILER

#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/constant_iterator.h>
#include <openacc.h> // acc_malloc などのため

size_t count_frequencies(long long* d_input_keys, int N,
        long long* d_unique_keys, long long* d_counts
    ) {
    // 1. OpenACCの生ポインタをThrustのdevice_ptrにラップ
    thrust::device_ptr<long long> th_keys(d_input_keys);

    // 2. GPU上で高速ソート (アトミックなし、Radix Sort)
    thrust::sort(thrust::device, th_keys, th_keys + N);

    thrust::device_ptr<long long> th_o_keys(d_unique_keys);
    thrust::device_ptr<long long> th_o_counts(d_counts);

    // 4. 隣り合う同じキーを数え上げる (Run-Length Encoding的な集計)
    auto new_end = thrust::reduce_by_key(
        thrust::device,
        th_keys, th_keys + N,                     // 入力: ソート済みキー
        thrust::make_constant_iterator(1LL),      // 入力: すべての値に「1」を供給
        th_o_keys,                                // 出力: ユニークキー
        th_o_counts                               // 出力: 出現回数
    );

    // 一意になった要素数を計算
    return new_end.first - th_o_keys;
}
#endif

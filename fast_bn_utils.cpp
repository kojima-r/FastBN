#include "fast_bn_utils.hpp"

#ifdef __NVCOMPILER

#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/constant_iterator.h>
#include <openacc.h> // acc_malloc などのため

size_t count_frequencies(long long* d_input_keys, int N,
        long long* d_unique_keys, long long* d_counts) {

    thrust::device_ptr<long long> th_keys(d_input_keys);
    thrust::sort(thrust::device, th_keys, th_keys + N);

    thrust::device_ptr<long long> th_o_keys(d_unique_keys);
    thrust::device_ptr<long long> th_o_counts(d_counts);

    // ソート済みの th_keys すべてに 1 をあたえて、キーごとにまとめると出現回数が得られる
    auto new_end = thrust::reduce_by_key(
        thrust::device,
        th_keys, th_keys + N,
        thrust::make_constant_iterator(1LL),
        th_o_keys,
        th_o_counts
    );

    return new_end.first - th_o_keys;
}
#endif

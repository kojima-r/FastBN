#pragma once
#include <vector>

// 外部から呼ぶための関数宣言
// d_input_keys 大きさ N の配列
// d_unique_keys 配列に格納されている値
// d_counts 配列に格納されている値の個数
size_t count_frequencies(long long* d_input_keys, int N, long long* d_unique_keys, long long* d_counts);
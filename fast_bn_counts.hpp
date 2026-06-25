#pragma once

#include <vector>

#include "fast_bn_dataset.hpp"

#ifdef __NVCOMPILER
#include <cuda_runtime.h>
#endif

//==================== カウント（n_ijk/n_ij） & スコア ====================
static int count_max_q = 128;
static int count_max_r = 32;

struct Counts {
    // n_ijk: 親配置 j と 子状態 k の同時度数
    // n_ij : 親配置 j の度数合計
    // サイズ: n_ijk は q_i * r_i, n_ij は q_i
    std::vector<long long> n_ijk;
    std::vector<long long> n_ij;
    int q_i=0;
    int r_i=0;
    Counts(){
        n_ijk.resize(count_max_q*count_max_r);
        n_ij.resize(count_max_q);
        long long* nij  = n_ij.data();
        long long* nijk = n_ijk.data();
        int q = count_max_q;
        int r = count_max_r;
        #pragma acc enter data create(nijk[0:q*r],nij[0:q])
    }
    Counts(const Counts& other){
        n_ijk.resize(count_max_q*count_max_r);
        n_ij.resize(count_max_q);
        long long* nij  = n_ij.data();
        long long* nijk = n_ijk.data();
        int q = count_max_q;
        int r = count_max_r;
        #pragma acc enter data create(nijk[0:q*r],nij[0:q])
        *this = other;
        // 代入演算子で update device 済
    }
    Counts(
        std::vector<long long>&& _n_ijk,
        std::vector<long long>&& _n_ij,
        const int _q_i,
        const int _r_i
    )
    {
        n_ijk.resize(count_max_q*count_max_r);
        n_ij.resize(count_max_q);
        long long* nij  = n_ij.data();
        long long* nijk = n_ijk.data();
        int q = count_max_q;
        int r = count_max_r;
        #pragma acc enter data create(nijk[0:q*r],nij[0:q])
        n_ijk = _n_ijk;
        n_ij = _n_ij;
        q_i = _q_i;
        r_i = _r_i;
        acc_update_device();
    }
    void assign(int q,int r){
        this->q_i = q;
        this->r_i = r;
        n_ijk.resize((size_t)q_i * r_i);
        n_ij.resize(q_i);
        long long* nij  = n_ij.data();
        long long* nijk = n_ijk.data();
#if defined __NVCOMPILER
        #pragma acc host_data use_device(nij, nijk)
        {
            // OpenACCの async(1) に紐づく CUDA ストリームを取得
            cudaStream_t stream = (cudaStream_t)acc_get_cuda_stream(1);
            // バイト単位で指定するため、要素数 × sizeof(型) を指定
            cudaMemsetAsync(nij,  0, q_i * sizeof(long long), stream);
            cudaMemsetAsync(nijk, 0, q_i * r_i * sizeof(long long), stream);
        }
#else
        {
            for(int i=0;i<q_i;i++){
                nij[i] = 0;
            }
            for(int i=0;i<q_i*r_i;i++){
                nijk[i] = 0;
            }
        }
#endif
    }
    Counts& operator=(const Counts& other) {
        this->q_i = other.q_i;
        this->r_i = other.r_i;
        this->n_ij = other.n_ij;
        this->n_ijk = other.n_ijk;
        acc_update_device();
        return *this;
    }
    ~Counts(){
        acc_delete();
    }
    void acc_update_device(void){
        if( q_i>0 && r_i> 0 ){
            long long* nij  = n_ij.data();
            long long* nijk = n_ijk.data();
            #pragma acc update device(nij[0:q_i], nijk[0:q_i*r_i])
        }
    }
    void acc_delete(void){
        long long* __restrict nij  = n_ij.data();
        long long* __restrict nijk = n_ijk.data();
        int q = count_max_q;
        int r = count_max_r;
        #pragma acc exit data delete(nijk[0:q*r],nij[0:q])
    }
    void acc_update_host(void){
        if( q_i>0 && r_i> 0 ){
            long long* __restrict nij  = n_ij.data();
            long long* __restrict nijk = n_ijk.data();
            #pragma acc update host(nij[0:q_i], nijk[0:q_i*r_i])
        }
    }
    void check_gpu(const char* tag=""){
#ifdef __NVCOMPILER
        long long* nij  = n_ij.data();
        long long* nijk = n_ijk.data();
        if(acc_is_present(nij,q_i*8)){
            std::cout << tag << " GPU: [PRESENT] nij" << std::endl;
        }else{
            std::cout << tag << " GPU: [NOT FOUND] nij" << std::endl;
        }
        if(acc_is_present(nijk,q_i*r_i*8)){
            std::cout << tag << " GPU: [PRESENT] nijk" << std::endl;
        }else{
            std::cout << tag << " GPU: [NOT FOUND] nijk" << std::endl;
        }
#endif
    }
};

// 親集合 Pa(i) に対する counts のフル再集計（O(N)）
// 本実装では、ADD のときは分割のため O(N) が不可避だが、REMOVE は完全インクリメンタル化する。
Counts computeCountsForNode_full(int i, const std::vector<int>& parents, const Dataset& ds, const std::vector<int>& radix);

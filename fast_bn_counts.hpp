#pragma once

#include <vector>

#include "fast_bn_dataset.hpp"

//==================== カウント（n_ijk/n_ij） & スコア ====================

struct Counts {
    // n_ijk: 親配置 j と 子状態 k の同時度数
    // n_ij : 親配置 j の度数合計
    // サイズ: n_ijk は q_i * r_i, n_ij は q_i
    std::vector<long long> n_ijk;
    std::vector<long long> n_ij;
    int q_i=0;
    int r_i=0;
    Counts(){
        std::cout << "construct Counts " << q_i << " " << r_i << std::endl;
        n_ijk.reserve(4096);
        n_ij.reserve(128);
        long long* nij  = n_ij.data();
        long long* nijk = n_ijk.data();
        #pragma acc enter data create(nijk[0:4096],nij[0:128])
    }
    Counts(const Counts& other){
        std::cout << "copy construct Counts " << q_i << " " << r_i << std::endl;
        n_ijk.reserve(4096);
        n_ij.reserve(128);
        long long* nij  = n_ij.data();
        long long* nijk = n_ijk.data();
        #pragma acc enter data create(nijk[0:4096],nij[0:128])
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
        n_ijk.reserve(4096);
        n_ij.reserve(128);
//        std::cout << "n_ijk address " << n_ijk.data() << std::endl;
//        std::cout << "n_ij address " << n_ij.data() << std::endl;
        long long* nij  = n_ij.data();
        long long* nijk = n_ijk.data();
        #pragma acc enter data create(nijk[0:4096],nij[0:128])
        n_ijk = _n_ijk;
        n_ij = _n_ij;
        q_i = _q_i;
        r_i = _r_i;
        std::cout << "construct Counts " << q_i << " " << r_i << std::endl;
//        std::cout << "n_ijk address " << n_ijk.data() << std::endl;
//        std::cout << "n_ij address " << n_ij.data() << std::endl;
        acc_update_device();
    }
    void assign(int q,int r){
        this->q_i = q;
        this->r_i = r;
        n_ijk.resize((size_t)q_i * r_i);
        std::fill(n_ijk.begin(), n_ijk.end(), 0);
        n_ij.resize(q_i);
        std::fill(n_ij.begin(), n_ij.end(), 0);
        acc_update_device();
    }
    Counts& operator=(const Counts& other) {
//        std::cout << "before copy address " << n_ij.data() << std::endl;
        this->q_i = other.q_i;
        this->r_i = other.r_i;
        this->n_ij = other.n_ij;
        this->n_ijk = other.n_ijk;
        std::cout << "copy operator " << q_i << " " << r_i << std::endl;
//        std::cout << "after copy address " << n_ij.data() << std::endl;
        acc_update_device();
        return *this;
    }
    ~Counts(){
        std::cout << "destruct Counts " << q_i << " " << r_i << std::endl;
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
        #pragma acc exit data delete(nij[0:128], nijk[0:4096])
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

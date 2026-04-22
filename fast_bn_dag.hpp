#pragma once

#include <cmath>
#include <vector>
#include <stack>
#include <cstdint>

#ifdef __NVCOMPILER
#include <openacc.h>
#endif

#include "fast_bn_dataset.hpp"

//==================== DAG（隣接行列）と到達性（サイクル判定） ====================

struct DAG {
    // 隣接行列（疎な大規模でも扱いやすく、実装簡単。密だとメモリは D^2）
    // 変数が数万で密行列は厳しいが、学習中の DAG は通常疎なので許容されることが多い。
    int D;
    std::vector<std::vector<char>> adj;   // adj[u][v] = 1 (u->v)
    std::vector<int> child_deg, parent_deg;

    DAG(int D=0): D(D), adj(D, std::vector<char>(D,0)), child_deg(D,0), parent_deg(D,0) {}

    bool hasEdge(int u,int v) const { return adj[u][v]; }
    int parentCount(int v) const { return parent_deg[v]; }
    int childCount(int u)  const { return child_deg[u]; }

    std::vector<int> parents(int v) const {
        // v の親ノード番号一覧を返す
        std::vector<int> p; p.reserve(8);
        for (int u=0;u<D;++u) if (adj[u][v]) p.push_back(u);
        return p;
    }

    void addEdge(int u,int v){ adj[u][v]=1; child_deg[u]++; parent_deg[v]++; }
    void removeEdge(int u,int v){ adj[u][v]=0; child_deg[u]--; parent_deg[v]--; }
    void reverseEdge(int u,int v){ adj[u][v]=0; child_deg[u]--; parent_deg[v]--; adj[v][u]=1; child_deg[v]++; parent_deg[u]++; }

    std::vector<std::pair<int,int>> edges() const {
        std::vector<std::pair<int,int>> e;
        for (int u=0;u<D;++u) for (int v=0;v<D;++v) if (adj[u][v]) e.emplace_back(u,v);
        return e;
    }
};

// 到達性（サイクル検出）
// - dense: 各ノードからの到達先集合をビットセットで保持→追加辺で増分更新（高速）
// - lazy : 毎回 DFS で確認→メモリ極小（既定）
struct Reachability {
    enum Mode { DENSE, LAZY } mode = LAZY;
    int D=0, W=0;                       // W = 64bit ワード数（D/64 切上げ）
    std::vector<std::vector<uint64_t>> reach;     // reach[u][word] のビットが立っていれば到達可
    const DAG* g = nullptr;

    Reachability() {}
    Reachability(const DAG& dag, Mode m): mode(m), D(dag.D), g(&dag) {
        if (mode==DENSE){
            W = (D + 63) >> 6;
            reach.assign(D, std::vector<uint64_t>(W, 0));
            // 初期化：直接辺のみセット（推移閉包は追加操作で徐々に積み上げる）
            for (int u=0;u<D;++u)
                for (int v=0;v<D;++v)
                    if (dag.adj[u][v]) reach[u][v>>6] |= (1ULL<<(v&63));
        }
    }

    inline bool testBit(const std::vector<uint64_t>& B, int v) const {
        // v ビットが立っているか
        return (B[v>>6] >> (v&63)) & 1ULL;
    }

    // 辺追加時の増分更新（dense のみ）
    // u->v を追加したら reach[u] に {v} と reach[v] を取り込み、
    // さらに「u に到達可能な全ノード w」にも reach[u] を OR で波及
    void onAddEdge(const DAG& dag, int u, int v){
        if (mode==LAZY) return; // lazy は状態を持たない
        auto onehot = std::vector<uint64_t>(W,0); onehot[v>>6] |= (1ULL<<(v&63));
        std::vector<uint64_t> addv = reach[v];
        if (reach[u].empty()) reach[u].assign(W,0);
        bool changed=false;
        for (int w=0; w<W; ++w){
            uint64_t nv = reach[u][w] | addv[w] | onehot[w];
            if (nv != reach[u][w]) { reach[u][w]=nv; changed=true; }
        }
        if (!changed) return;
        // u に到達可能な全ノード w にも reach[u] を反映
        for (int w=0; w<dag.D; ++w){
            if (w==u) continue;
            if (testBit(reach[w], u)) {
                for (int k=0;k<W;++k) {
                    uint64_t nv = reach[w][k] | reach[u][k];
                    if (nv != reach[w][k]) reach[w][k] = nv;
                }
            }
        }
    }
    // 辺削除時（dense）に reach を厳密に縮めるのは高コスト。
    // ここでは保守的運用（縮めない）。false positive により一部 ADD/REVERSE が弾かれる可能性はあるが、
    // 計算量を重視して簡略化している（必要なら定期的な再構築オプションを追加可能）。
    void onRemoveEdge(const DAG&, int, int){ /* keep conservative */ }

    void onReverseEdge(const DAG& dag, int u, int v){
        if (mode==LAZY) return;
        onRemoveEdge(dag, u, v);
        onAddEdge(dag, v, u);
    }

    // u->v を追加してサイクルになるか？（true=追加OK）
    bool canAddAcyclic(const DAG& dag, int u, int v) const {
        if (u==v || dag.adj[u][v]) return false;
        if (mode==DENSE) return !testBit(reach[v], u);
        // LAZY: v から DFS して u に到達するか
        std::vector<char> vis(dag.D,0);
        std::stack<int> st; st.push(v);
        while(!st.empty()){
            int x=st.top(); st.pop();
            if (x==u) return false; // サイクル
            if (vis[x]) continue;
            vis[x]=1;
            for (int y=0;y<dag.D;++y) if (dag.adj[x][y]) st.push(y);
        }
        return true;
    }

    // u->v の反転（v->u 追加）でサイクルにならないか
    bool canReverseAcyclic(const DAG& dag, int u, int v) const {
        if (!dag.adj[u][v] || u==v) return false;
        if (mode==DENSE) return !testBit(reach[u], v); // 保守的判定
        // LAZY: 一時的に u->v を外して v->u 追加の可否を DFS で確認
        const_cast<DAG&>(dag).adj[u][v]=0;
        bool ok = canAddAcyclic(dag, v, u);
        const_cast<DAG&>(dag).adj[u][v]=1;
        return ok;
    }
};

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

// 親集合の混合基数インデックス（右側の親が下位桁）
// 例: parents = [p0,p1,p2], 基数 r[p0], r[p1], r[p2]
//     j = x[p0]*r[p1]*r[p2] + x[p1]*r[p2] + x[p2]
inline int mixedRadixIndex(const std::vector<int>& parents, const std::vector<int>& r, const std::vector<int>& row){
    int qidx=0, mult=1;
    for (int t=(int)parents.size()-1; t>=0; --t){
        int p=parents[t];
        qidx += row[p]*mult;
        mult *= r[p];
    }
    return qidx;
}

// 親集合 parents に対する右からの混合基数(radix)を前計算する。
//  parents: 親ノードのインデックス（例: [p0, p1, ..., p_{P-1}]）
//  ds.r   : 各変数の取りうる値の数 r[i]
//  radix  : 出力 (サイズ P)。j = Σ ds.x(n, parents[t]) * radix[t] で使用。
static inline void build_mixed_radix(const std::vector<int>& parents,
                                     const Dataset& ds,
                                     std::vector<int>& radix) noexcept
{
    const int P = (int)parents.size();
    if (P <= 0) {
        radix.clear();
        return;
    }
    // 既存容量を再利用しつつ、サイズだけ合わせる
    radix.assign(P, 1);
    const int* ds_r_ptr = ds.r.data();
    int* rdx_ptr = radix.data();
    const int* pa_ptr = parents.data();
    // 右から順に混合基数を構成
    for (int t = P - 2; t >= 0; --t) {
        rdx_ptr[t] = rdx_ptr[t + 1] * ds_r_ptr[pa_ptr[t + 1]];
    }
    #pragma acc enter data copyin(rdx_ptr[0:P], pa_ptr[0:P])
}

// 1サンプル n について、parents の値から混合基数インデックス j を計算。
//  事前に build_mixed_radix(parents, r, radix) 済みであることが前提。
#pragma acc routine seq
//static inline int mixed_radix_index_row(const Dataset& ds,
//                                        int n,
//                                        const std::vector<int>& parents,
//                                        const std::vector<int>& radix) noexcept
static inline int mixed_radix_index_row(
  int D,
  const int* ds_ptr,
  int n,
  int P,
  const int* pa,
  const int* rdx) noexcept
{
    int j = 0;
//    const int P = (int)parents.size();
//    const int D = ds.D;
    if( P<=0 ) return j;
//    const int* ds_ptr = ds.X_flat.data();
//    const int* rdx = radix.data();
//    const int* pa = parents.data();
    const size_t offset = (size_t)n * D;
//    #pragma acc parallel loop reduction(+:j) present(pa, ds_ptr, rdx)
    for (int t = 0; t < P; ++t) {
        const int p = pa[t];
        j += ds_ptr[offset + p] * rdx[t];
        // j += ds.x(n, p) * radix[t];   // ★フラット配列アクセス
    }
    return j;
}

// 親集合 Pa(i) に対する counts のフル再集計（O(N)）
// 本実装では、ADD のときは分割のため O(N) が不可避だが、REMOVE は完全インクリメンタル化する。
Counts computeCountsForNode_full(int i, const std::vector<int>& parents, const Dataset& ds, const std::vector<int>& radix);

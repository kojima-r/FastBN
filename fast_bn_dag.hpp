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

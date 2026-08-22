#include <bits/stdc++.h>
#include <random>
#include <chrono>
#include <tuple>

#include <cstdlib> // getenv()

using namespace std;

// ===== chi-square 上側確率 Q(k/2, x/2) を計算するための補助 =====
static inline double gammaln(double x) { return std::lgamma(x); }

// 正規化下側不完全ガンマ P(s,x)（級数展開）と
// 正規化上側不完全ガンマ Q(s,x)（連分数）を使い分け
static double regularized_gamma_P(double s, double x) {
    if (x <= 0.0) return 0.0;
    const int MAXIT = 1000;
    const double EPS = 1e-14;

    // 級数: P(s,x) = e^{-x} x^s / Γ(s) * sum_{n=0..∞} (x^n / Γ(s+n+1))
    double sum = 1.0 / s;
    double term = sum;
    for (int n=1; n<MAXIT; ++n) {
        term *= x / (s + n);
        sum += term;
        if (std::fabs(term) < std::fabs(sum) * EPS) break;
    }
    return std::exp(s*std::log(x) - x - gammaln(s)) * sum;
}

static double regularized_gamma_Q(double s, double x) {
    if (x <= 0.0) return 1.0;
    const int MAXIT = 1000;
    const double EPS = 1e-14;

    // 連分数: Q(s,x) を直接計算（Lentz 法）
    double C = 1.0 / 1e-300;
    double D = 0.0;
    double f = 0.0;

    // 初期化
    double a0 = 0.0;
    double b0 = x + 1.0 - s;
    D = 1.0 / std::max(1e-300, b0);
    C = std::max(1e-300, b0);
    f = D;

    for (int i=1; i<MAXIT; ++i) {
        double a = i * (s - i);
        double b = b0 + 2.0 * i;

        // D = 1 / (b + a*D)
        D = 1.0 / std::max(1e-300, b + a*D);
        C = std::max(1e-300, b + a/C);
        double delta = C * D;
        f *= delta;
        if (std::fabs(delta - 1.0) < EPS) break;
    }

    double pref = std::exp(s*std::log(x) - x - gammaln(s));
    return pref * f;
}

// P と Q の切替（数値安定化）
static double gamma_reg_P(double s, double x) {
    if (x < s + 1.0) return regularized_gamma_P(s, x);
    return 1.0 - regularized_gamma_Q(s, x);
}
static double gamma_reg_Q(double s, double x) {
    if (x < s + 1.0) return 1.0 - regularized_gamma_P(s, x);
    return regularized_gamma_Q(s, x);
}

// カイ二乗の上側確率（p 値）: df=k, statistic=chi2
static double chisq_p_upper(double chi2, int df) {
    if (chi2 < 0.0) return 1.0;
    double s = 0.5 * df;
    double x = 0.5 * chi2;
    return gamma_reg_Q(s, x); // p = Q(df/2, chi2/2)
}

/*
  ベイジアンネットワーク構造学習（スコアベース）
  - スコア: BIC / K2 / BDeu
  - 探索: Hill Climbing + Tabu Search（タブーテナー指定で有効化）
  - 高速化:
      (1) 候補親 K の前処理（MI: 相互情報量に基づく上位Kの候補親の絞り込み）
      (3) 到達可能性（サイクル判定）の効率化
          - dense: ビットセット推移閉包（増分更新・高速／中〜大メモリ）
          - lazy : DFS ベース（省メモリ／やや遅い）
      (4) 親配置インデックス j_index の LRU キャッシュ
          - ADD の試行時に O(N) で新カウントを「分割」再構成
  - さらに本版では **REMOVE を完全インクリメンタル**に実装
      - 混合基数インデックスの「桁落とし」により O(q·r_i) で counts をマージ（データ再走査不要）

  データ前提:
    - CSV 各セルは離散化済みの 0,1,2,... の整数
    - 欠損無しを想定（必要なら拡張）

  メモリ・速度の考え方（大規模 D=10k〜30k を想定）:
    - reach: 既定は lazy（DFS）で省メモリ。dense は D^2/64 ワードを要するため数万では重い場合あり。
    - j_index: N×int を全ノードに持つと巨大なので LRU（--jindex-cache）で保持ノード数を制限
    - topK: 候補親 K の前処理で探索分岐を大きく削減（--mi-sample と --mi-budget で計算量制御）
    - max-parents / max-children: 構造の疎性を強制することで q（親配置数）を抑制 → カウント配列が爆発しない
*/

//==================== 基本定義 ====================

enum class ScoreType { BIC, K2, BDeu };
struct Dataset {
    //vector<vector<int>> X;     // N×D データ本体
    vector<int> X_flat;
    
    int N=0, D=0;
    vector<int> r;             // 各列の基数
    vector<string> var_names;  // 変数名（CSV/TSV のヘッダ）
    // accessor
    inline int x(int n, int d) const noexcept {
        return X_flat[(size_t)n * D + d];
    }
    // CSV または TSV を読み込み（区切り文字を自動判定）
    static Dataset fromCSV(const string& path) {
        ifstream fin(path);
        if (!fin) throw runtime_error("Failed to open file: " + path);

        string line;
        vector<string> headers;
        vector<vector<int>> rows;
        bool first = true;
        char delim = ','; // デフォルトはカンマ

        while (getline(fin, line)) {
            if (line.empty()) continue;

            // 区切り文字の自動判定：最初の行の '\t' の有無で決める
            if (first) {
                if (line.find('\t') != string::npos)
                    delim = '\t';
            }

            stringstream ss(line);
            string tmp;

            if (first) {
                // --- 1行目はヘッダ ---
                while (getline(ss, tmp, delim)) {
                    // 前後の空白除去
                    tmp.erase(tmp.begin(), find_if(tmp.begin(), tmp.end(), [](unsigned char c){return !isspace(c);} ));
                    tmp.erase(find_if(tmp.rbegin(), tmp.rend(), [](unsigned char c){return !isspace(c);} ).base(), tmp.end());
                    headers.push_back(tmp);
                }
                first = false;
                continue;
            }

            // --- データ行 ---
            vector<int> row;
            while (getline(ss, tmp, delim)) {
                tmp.erase(tmp.begin(), find_if(tmp.begin(), tmp.end(), [](unsigned char c){return !isspace(c);} ));
                tmp.erase(find_if(tmp.rbegin(), tmp.rend(), [](unsigned char c){return !isspace(c);} ).base(), tmp.end());
                if (tmp.empty()) throw runtime_error("Empty cell in data row.");
                try {
                    row.push_back(stoi(tmp));
                } catch (...) {
                    throw runtime_error("Non-integer cell detected: \"" + tmp + "\"");
                }
            }
            if (!row.empty())
                rows.push_back(move(row));
        }

        if (rows.empty()) throw runtime_error("Empty dataset (no data rows).");
        int D = (int)rows[0].size();
        for (auto& r: rows)
            if ((int)r.size()!=D)
                throw runtime_error("Inconsistent column count in data rows.");

        if (!headers.empty() && (int)headers.size()!=D)
            throw runtime_error("Header count does not match column count.");

        // 各列の基数（最大値+1）
        vector<int> rcard(D,0);
        for (int j=0;j<D;++j){
            int mx=0;
            for (auto& row: rows) mx=max(mx, row[j]);
            rcard[j]=mx+1;
        }

        Dataset ds;
        //ds.X = move(rows);
        ds.N = (int)rows.size();
        ds.D = D;
        ds.r = move(rcard);
        ds.var_names = move(headers);
        ds.X_flat.resize((size_t)ds.N * ds.D);
        for (int n = 0; n < ds.N; ++n) {
            const auto &row = rows[n];
            std::copy(row.begin(), row.end(),
                      ds.X_flat.begin() + (size_t)n * ds.D);
        }

        cerr << "[info] Loaded " << ds.N << " samples, "
             << ds.D << " variables, delimiter='" << (delim=='\t' ? "\\t" : ",") << "'\n";
        return ds;
    }
};

static Dataset bootstrapResampleDataset(const Dataset& ds, std::mt19937_64& rng) {
    Dataset out;
    out.D = ds.D;
    out.N = ds.N;
    out.r = ds.r;
    out.var_names = ds.var_names;
    //out.X.resize(ds.N, std::vector<int>(ds.D, 0));
    out.X_flat.resize((size_t)out.N * out.D);

    std::uniform_int_distribution<int> uid(0, ds.N - 1);
    for (int n = 0; n < ds.N; ++n) {
        int src = uid(rng);          // 復元抽出
        //out.X[n] = ds.X[src];        // 行コピー
	std::copy(ds.X_flat.begin() + (size_t)src * ds.D,
                  ds.X_flat.begin() + (size_t)(src + 1) * ds.D,
                  out.X_flat.begin() + (size_t)n * ds.D);
    }
    return out;
}

// （疎）エッジ出現カウントのキー化ユーティリティ
static inline uint64_t edgeKeyUV(uint32_t u, uint32_t v) {
    return ( (uint64_t)u << 32 ) | (uint64_t)v;
}
//==================== DAG（隣接行列）と到達性（サイクル判定） ====================

struct DAG {
    // 隣接行列（疎な大規模でも扱いやすく、実装簡単。密だとメモリは D^2）
    // 変数が数万で密行列は厳しいが、学習中の DAG は通常疎なので許容されることが多い。
    int D;
    vector<vector<char>> adj;   // adj[u][v] = 1 (u->v)
    vector<int> child_deg, parent_deg;

    DAG(int D=0): D(D), adj(D, vector<char>(D,0)), child_deg(D,0), parent_deg(D,0) {}

    bool hasEdge(int u,int v) const { return adj[u][v]; }
    int parentCount(int v) const { return parent_deg[v]; }
    int childCount(int u)  const { return child_deg[u]; }

    vector<int> parents(int v) const {
        // v の親ノード番号一覧を返す
        vector<int> p; p.reserve(8);
        for (int u=0;u<D;++u) if (adj[u][v]) p.push_back(u);
        return p;
    }

    void addEdge(int u,int v){ adj[u][v]=1; child_deg[u]++; parent_deg[v]++; }
    void removeEdge(int u,int v){ adj[u][v]=0; child_deg[u]--; parent_deg[v]--; }
    void reverseEdge(int u,int v){ adj[u][v]=0; child_deg[u]--; parent_deg[v]--; adj[v][u]=1; child_deg[v]++; parent_deg[u]++; }

    vector<pair<int,int>> edges() const {
        vector<pair<int,int>> e;
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
    vector<vector<uint64_t>> reach;     // reach[u][word] のビットが立っていれば到達可
    const DAG* g = nullptr;

    Reachability() {}
    Reachability(const DAG& dag, Mode m): mode(m), D(dag.D), g(&dag) {
        if (mode==DENSE){
            W = (D + 63) >> 6;
            reach.assign(D, vector<uint64_t>(W, 0));
            // 初期化：直接辺のみセット（推移閉包は追加操作で徐々に積み上げる）
            for (int u=0;u<D;++u)
                for (int v=0;v<D;++v)
                    if (dag.adj[u][v]) reach[u][v>>6] |= (1ULL<<(v&63));
        }
    }

    inline bool testBit(const vector<uint64_t>& B, int v) const {
        // v ビットが立っているか
        return (B[v>>6] >> (v&63)) & 1ULL;
    }

    // 辺追加時の増分更新（dense のみ）
    // u->v を追加したら reach[u] に {v} と reach[v] を取り込み、
    // さらに「u に到達可能な全ノード w」にも reach[u] を OR で波及
    void onAddEdge(const DAG& dag, int u, int v){
        if (mode==LAZY) return; // lazy は状態を持たない
        auto onehot = vector<uint64_t>(W,0); onehot[v>>6] |= (1ULL<<(v&63));
        vector<uint64_t> addv = reach[v];
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
        vector<char> vis(dag.D,0);
        stack<int> st; st.push(v);
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
    vector<long long> n_ijk;
    vector<long long> n_ij;
    int q_i=0, r_i=0;
};

// 親集合の混合基数インデックス（右側の親が下位桁）
// 例: parents = [p0,p1,p2], 基数 r[p0], r[p1], r[p2]
//     j = x[p0]*r[p1]*r[p2] + x[p1]*r[p2] + x[p2]
inline int mixedRadixIndex(const vector<int>& parents, const vector<int>& r, const vector<int>& row){
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
//  r      : 各変数の取りうる値の数 r[i]
//  radix  : 出力 (サイズ P)。j = Σ ds.x(n, parents[t]) * radix[t] で使用。
static inline void build_mixed_radix(const std::vector<int>& parents,
                                     const std::vector<int>& r,
                                     std::vector<int>& radix) noexcept
{
    const int P = (int)parents.size();
    if (P <= 0) {
        radix.clear();
        return;
    }
    // 既存容量を再利用しつつ、サイズだけ合わせる
    radix.assign(P, 1);
    // 右から順に混合基数を構成
    for (int t = P - 2; t >= 0; --t) {
        radix[t] = radix[t + 1] * r[parents[t + 1]];
    }
}

// 1サンプル n について、parents の値から混合基数インデックス j を計算。
//  事前に build_mixed_radix(parents, r, radix) 済みであることが前提。
static inline int mixed_radix_index_row(const Dataset& ds,
                                        int n,
                                        const std::vector<int>& parents,
                                        const std::vector<int>& radix) noexcept
{
    const int P = (int)parents.size();
    int j = 0;
    for (int t = 0; t < P; ++t) {
        const int p = parents[t];
        j += ds.x(n, p) * radix[t];   // ★フラット配列アクセス
    }
    return j;
}

// 親集合 Pa(i) に対する counts のフル再集計（O(N)）
// 本実装では、ADD のときは分割のため O(N) が不可避だが、REMOVE は完全インクリメンタル化する。
Counts computeCountsForNode_full(int i, const vector<int>& parents, const Dataset& ds){
    const int N   = ds.N;
    const int r_i = ds.r[i];
    int q_i = 1;
    for (int p: parents) q_i *= ds.r[p];
    
    vector<long long> nij(q_i,0), nijk((size_t)q_i*r_i,0);

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

struct Scorer {
    const Dataset& ds;
    ScoreType type;
    double ess; // BDeu の等価事例数（他スコアでは未使用）

    Scorer(const Dataset& ds, ScoreType t, double ess=1.0): ds(ds), type(t), ess(ess) {}

    // BIC: 対数尤度 - (d/2)*log(N), d=(r_i-1)*q_i
    double nodeScoreBIC(const Counts& c) const {
        double ll=0.0;
        for (int j=0;j<c.q_i;++j){
            double nij=(double)c.n_ij[j];
            if (nij<=0) continue;
            for (int k=0;k<c.r_i;++k){
                long long nijk=c.n_ijk[(size_t)j*c.r_i+k];
                if (nijk==0) continue;
                ll += nijk * (log((double)nijk) - log(nij));
            }
        }
        int d=(c.r_i-1)*c.q_i;
        double pen = 0.5 * d * log((double)max(1, ds.N));
        return ll - pen;
    }

    // K2: Dirichlet(1) 事前
    double nodeScoreK2(const Counts& c) const {
        double s=0.0;
        for (int j=0;j<c.q_i;++j){
            double nij=(double)c.n_ij[j];
            s += lgamma((double)c.r_i) - lgamma(nij + (double)c.r_i);
            for (int k=0;k<c.r_i;++k){
                double nijk=(double)c.n_ijk[(size_t)j*c.r_i + k];
                s += lgamma(nijk + 1.0); // - lgamma(1) = 0
            }
        }
        return s;
    }

    // BDeu: 一様ハイパーパラメータ（等価事例数 ess を q_i, r_i に均等割）
    double nodeScoreBDeu(const Counts& c) const {
        double s=0.0;
        double alpha_ij = ess / (double)max(1,c.q_i);
        double alpha_ijk = alpha_ij / (double)c.r_i;
        for (int j=0;j<c.q_i;++j){
            double nij=(double)c.n_ij[j];
            s += lgamma(alpha_ij) - lgamma(nij + alpha_ij);
            for (int k=0;k<c.r_i;++k){
                double nijk=(double)c.n_ijk[(size_t)j*c.r_i + k];
                s += lgamma(nijk + alpha_ijk) - lgamma(alpha_ijk);
            }
        }
        return s;
    }

    double nodeScore(const Counts& c) const {
        switch(type){
            case ScoreType::BIC: return nodeScoreBIC(c);
            case ScoreType::K2:  return nodeScoreK2(c);
            case ScoreType::BDeu:return nodeScoreBDeu(c);
        }
        return -INFINITY;
    }
};

//==================== j_index LRU キャッシュ ====================

/*
  j_index[v][n] = 「現在の親集合 Pa(v) に基づくサンプル n の親配置インデックス j」
  - ADD 試行時に j' = j * r[u] + x_u[n] で O(1) 更新できるため、O(N) で新カウントを分割生成可能。
  - ただし全ノード v について N 要素の配列を常駐させるとメモリが大きいので、
    LRU で保持ノード数を制限（--jindex-cache）。
*/
struct JIndexCache {
    struct Entry { int v; vector<int> j; list<int>::iterator it; };
    int cap; // 最大保持ノード数（0 ならキャッシュ無効）
    const Dataset* ds;
    const DAG* g;
    unordered_map<int, unique_ptr<Entry>> mp;
    list<int> lru; // front が最新

    JIndexCache(int cap, const Dataset* ds, const DAG* g): cap(cap), ds(ds), g(g) {}

    void touch_(int v){
        if (!mp.count(v)) return;
        lru.erase(mp[v]->it);
        lru.push_front(v);
        mp[v]->it = lru.begin();
    }
    void evictIfNeeded_(){
        while (cap>0 && (int)mp.size()>cap){
            int victim = lru.back(); lru.pop_back();
            mp.erase(victim);
        }
    }
    void invalidate(int v){
        if (cap==0) return;
        if (!mp.count(v)) return;
        lru.erase(mp[v]->it);
        mp.erase(v);
    }
    void onParentsChanged(const vector<int>& vs){ if (cap==0) return; for (int v: vs) invalidate(v); }

    // 親集合 Pa(v) に合わせて j_index を構築（O(N·|Pa|)）
    void build(int v, vector<int>& out) const {
        auto pa = g->parents(v);
	const int N = ds->N;
        out.assign(N, 0);
        if (pa.empty()) return;
        std::vector<int> radix;
        build_mixed_radix(pa, ds->r, radix);
        for (int n=0;n<N;++n){
	    const int j = mixed_radix_index_row(*ds, n, pa, radix);
            out[n]=j;
        }
    }

    // 取得：キャッシュにあれば返す。無ければ構築してキャッシュ。
    const vector<int>& get(int v){
        if (cap==0){
            // 省メモリモード：都度 build して一時オブジェクトを返す（実用では cap>0 を推奨）
            static vector<int> tmp;
            build(v, tmp);
            // 注意: 本簡易実装では cap==0 の場合にメモリリークを許容（サンプル用途）
            return *(new vector<int>(tmp));
        }
        if (mp.count(v)) { touch_(v); return mp[v]->j; }
        auto e = make_unique<Entry>();
        e->v = v;
        build(v, e->j);
        lru.push_front(v);
        e->it = lru.begin();
        mp[v] = move(e);
        evictIfNeeded_();
        return mp[v]->j;
    }
};

//==================== 候補親Kの前処理（MI: 相互情報量） ====================

/*
  目的:
    - 各ノード v に対して、相互情報量 MI(u; v) が大きい上位 K 個の u を候補親に限定
    - 検索空間を大幅に縮小

  計算量削減の工夫:
    - --mi-sample S: 行サンプル数 S のみで MI を近似（0 なら全行）
    - --mi-budget B: 1ノードあたり B 変数のみを評価（0 なら全変数）
*/
enum class CandMetric { MI, CHI2 };

// rows: サンプリングした行インデックス（全件なら 0..N-1）
struct AssocCandidates {
    // MI（nats）
    static double mutual_info_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows) {
        int ru = ds.r[u], rv = ds.r[v];
        std::vector<long long> cu(ru,0), cv(rv,0);
        std::unordered_map<long long, long long> cuv; cuv.reserve(rows.size()*2);
        auto key = [&](int a,int b)->long long { return (long long)a * (long long)rv + b; };

        for (int idx: rows){
            int a = ds.x(idx,u), b = ds.x(idx,v);
            ++cu[a]; ++cv[b]; ++cuv[key(a,b)];
        }
        double N = (double)rows.size();
        double mi=0.0;
        for (auto &kv : cuv) {
            long long ab = kv.second;
            int a = (int)(kv.first / rv);
            int b = (int)(kv.first % rv);
            if (ab==0) continue;
            double pab = ab / N;
            double pa  = cu[a] / N;
            double pb  = cv[b] / N;
            mi += pab * std::log( (pab + 1e-300) / (pa * pb + 1e-300) );
        }
        return mi; // nats
    }

    // カイ二乗の p 値（独立性検定、上側確率）と統計量
    static std::pair<double,double> chi2_p_pair(const Dataset& ds, int u, int v, const std::vector<int>& rows) {
        int ru = ds.r[u], rv = ds.r[v];
        std::vector<long long> cu(ru,0), cv(rv,0);
        std::vector<long long> tab((size_t)ru * rv, 0);

        for (int idx: rows){
            int a = ds.x(idx,u), b = ds.x(idx,v);
            ++cu[a]; ++cv[b];
            ++tab[(size_t)a * rv + b];
        }

        double N = (double)rows.size();
        if (N <= 0.0) return {1.0, 0.0};

        // カイ二乗統計量
        double chi2 = 0.0;
        for (int a=0; a<ru; ++a) {
            for (int b=0; b<rv; ++b) {
                double obs = (double)tab[(size_t)a * rv + b];
                double expv = (double)cu[a] * (double)cv[b] / N;
                if (expv > 0.0) {
                    double diff = obs - expv;
                    chi2 += diff * diff / expv;
                }
            }
        }
        int df = (ru - 1) * (rv - 1);
        if (df <= 0) return {1.0, 0.0};
        double p = chisq_p_upper(chi2, df);
        return {p, chi2};
    }

    // metric に応じて候補親を返す
    //  - metric=MI    : mi >= mi_threshold_nats を残し、K>0なら MI 降順で上位K
    //  - metric=CHI2  : p <= chi2_p_threshold を残し、K>0なら -log(p) 降順で上位K
    static std::vector<std::vector<int>> compute(const Dataset& ds,
                                                 int K,
                                                 int budget,
                                                 const std::vector<int>& rows,
                                                 std::mt19937_64& rng,
                                                 CandMetric metric,
                                                 double mi_threshold_nats,
                                                 double chi2_p_threshold)
    {
        int D = ds.D;
        std::vector<std::vector<int>> out(D);

        for (int v=0; v<D; ++v) {
            // 候補母集団（budget で間引き）
            std::vector<int> cand;
            cand.reserve(budget>0? budget : (D-1));
            if (budget>0 && budget < D-1) {
                std::unordered_set<int> used; used.reserve(budget*2);
                while ((int)cand.size()<budget) {
                    int u = (int)(rng()%D);
                    if (u==v) continue;
                    if (used.insert(u).second) cand.push_back(u);
                }
            } else {
                for (int u=0; u<D; ++u) if (u!=v) cand.push_back(u);
            }

            if (metric == CandMetric::MI) {
                std::vector<std::pair<double,int>> scored; // (mi,u)
                scored.reserve(cand.size());
                for (int u: cand) {
                    double mi = mutual_info_pair(ds, u, v, rows);
                    if (mi >= mi_threshold_nats) scored.emplace_back(mi, u);
                }
                if (K>0 && (int)scored.size()>K) {
                    std::nth_element(scored.begin(), scored.begin()+K, scored.end(),
                        [](auto& a, auto& b){ return a.first > b.first; });
                    scored.resize(K);
                }
                std::vector<int> keep; keep.reserve(scored.size());
                for (auto& p: scored) keep.push_back(p.second);
                std::sort(keep.begin(), keep.end());
                out[v] = std::move(keep);
            } else { // CHI2
                std::vector<std::pair<double,int>> scored; // (score,u), score=-log(p)
                scored.reserve(cand.size());
                for (int u: cand) {
                    auto [p, chi2] = chi2_p_pair(ds, u, v, rows);
                    if (p <= chi2_p_threshold) {
                        double score = -std::log(std::max(p, 1e-300)); // 小pほど大
                        scored.emplace_back(score, u);
                    }
                }
                if (K>0 && (int)scored.size()>K) {
                    std::nth_element(scored.begin(), scored.begin()+K, scored.end(),
                        [](auto& a, auto& b){ return a.first > b.first; });
                    scored.resize(K);
                }
                std::vector<int> keep; keep.reserve(scored.size());
                for (auto& p: scored) keep.push_back(p.second);
                std::sort(keep.begin(), keep.end());
                out[v] = std::move(keep);
            }
        }
        return out;
    }
};

struct MICandidates {
    static double mutual_info_pair(const Dataset& ds, int u, int v, const vector<int>& rows){
        int ru = ds.r[u], rv = ds.r[v];
        vector<long long> cu(ru,0), cv(rv,0);
        unordered_map<long long, long long> cuv; cuv.reserve(rows.size()*2);
        auto key = [&](int a,int b)->long long { return (long long)a * (long long)rv + b; };

        for (int idx: rows){
            int au = ds.x(idx,u), bv = ds.x(idx,v);
            ++cu[au]; ++cv[bv]; ++cuv[key(au,bv)];
        }

        double N = (double)rows.size();
        double mi=0.0;
        for (auto &kv : cuv){
            long long ab = kv.second;
            int a = (int)(kv.first / rv);
            int b = (int)(kv.first % rv);
            if (ab==0) continue;
            double pab = ab / N;
            double pa  = cu[a] / N;
            double pb  = cv[b] / N;
            // 数値安定用に 1e-300 を加える
            mi += pab * log( (pab + 1e-300) / (pa * pb + 1e-300) );
        }
        return mi;
    }

    static vector<vector<int>> compute(const Dataset& ds, int K, int budget, const vector<int>& rows, mt19937_64& rng, double mi_threshold_nats){
        int D = ds.D;
        vector<vector<int>> topK(D);

        for (int v=0; v<D; ++v){
            // 変数サンプリング：自身以外から budget 個
            vector<int> cand;
            if (budget>0 && budget < D-1) {
                cand.reserve(budget);
                unordered_set<int> used; used.reserve(budget*2);
                while ((int)cand.size()<budget){
                    int u = (int)(rng()%D);
                    if (u==v) continue;
                    if (used.insert(u).second) cand.push_back(u);
                }
            } else {
                cand.reserve(D-1);
                for (int u=0;u<D;++u) if (u!=v) cand.push_back(u);
            }

            // (u, MI(u;v)) を計算
            vector<pair<double,int>> scored;
            scored.reserve(cand.size());
            for (int u: cand){
                double mi = mutual_info_pair(ds, u, v, rows); // nats
                if (mi >= mi_threshold_nats) {                 // しきい値でまず足切り
                    scored.emplace_back(mi, u);
                }
            }

            // K>0 のときは上位 K に丸める（MI 降順）
            if (K>0 && (int)scored.size() > K){
                nth_element(scored.begin(), scored.begin()+K, scored.end(),
                            [](const auto& a, const auto& b){ return a.first > b.first; });
                scored.resize(K);
            }

            // u のみ抽出し、binary_search 用にソート（昇順）
            vector<int> keep; keep.reserve(scored.size());
            for (auto &p: scored) keep.push_back(p.second);
            sort(keep.begin(), keep.end());
            topK[v] = move(keep);
        }
        return topK;
    }
};

//==================== 探索（HC/Tabu） ====================

struct HillClimber {
    const Dataset& ds;
    Scorer scorer;
    DAG g;
    Reachability reach;
    int max_iter = 100000;
    bool verbose = true;

    // 構造制約（次数上限）：q を小さく保つため重要（メモリ・計算量の上限に効く）
    int max_parents = 3;   // 既定: 大規模向け
    int max_children = 16;

    // タブー
    int tabu_tenure = 0;

    // 候補親（MI の上位 K）
    vector<vector<int>> candParents; // ソート済ベクタ
    int topK = 50;

    // j_index LRU キャッシュ
    JIndexCache jcache;

    // 現在のローカルカウント & スコア（各ノードで保持）
    // - ADD: 分割で新カウントを生成 → 更新
    // - REMOVE: **完全インクリメンタル**にマージ → 更新
    vector<Counts> nodeCounts;
    vector<double> nodeScoreNow;
    double totalNow = 0.0;

    HillClimber(const Dataset& ds, ScoreType t, double ess,
                const DAG& init, Reachability::Mode rmode,
                int jcache_cap)
        : ds(ds), scorer(ds,t,ess), g(init), reach(init, rmode), jcache(jcache_cap, &ds, &g)
    {
        nodeCounts.resize(ds.D);
        nodeScoreNow.resize(ds.D, 0.0);

        // 初期カウントを各ノードで構築（O(N·|Pa|) を D 回）
        for (int v=0; v<ds.D; ++v){
            auto Pa = g.parents(v);
            nodeCounts[v] = computeCountsForNode_full(v, Pa, ds);
            nodeScoreNow[v] = scorer.nodeScore(nodeCounts[v]);
            totalNow += nodeScoreNow[v];
        }
    }

    struct Move { enum Type { ADD, REMOVE, REVERSE, NONE } type=NONE; int u=-1, v=-1; double delta=0; };

    // 追加時の次数上限チェック
    bool addAllowedDegree(int u,int v) const {
        if (g.parentCount(v) >= max_parents) return false;
        if (g.childCount(u)  >= max_children) return false;
        return true;
    }

    bool addAllowed(const int u, const int v) const {
        if (!addAllowedDegree(u,v)) return false;
        if (!reach.canAddAcyclic(g,u,v)) return false;
        if (topK>0){
            const auto& C = candParents[v];
            if (!binary_search(C.begin(), C.end(), u)) return false;
        }
        return true;
    }

    bool reverseAllowed(int u,int v) const {
        if (!reach.canReverseAcyclic(g,u,v)) return false;
        if (g.parentCount(u) + 1 > max_parents) return false; // v->u の追加で u の親+1
        if (g.childCount(v)  + 1 > max_children) return false; // v の子+1
        if (topK>0){
            const auto& C = candParents[u];
            if (!binary_search(C.begin(), C.end(), v)) return false; // 逆向きも候補に含まれるか
        }
        return true;
    }

    //==================== Δスコア計算（インクリメンタル） ====================

    // ADD(u->v): 既存の j_index[v]（Pa(v) 用）を使い、q' = q * r_u の配列に「分割」して再構成（O(N)）
    // 返値: Δスコア（after - before）。生成した newC は呼び出し側で適用に使う。
    double deltaAdd_andBuildNewCounts(int u, int v, Counts& newC) {
        const Counts& curC = nodeCounts[v];
        const int N = ds.N;
        const int ru = ds.r[u];
        const int r_i = curC.r_i;

        // 親なしの場合は j_index を呼ばない
        if (curC.q_i == 0) {
            const int q = 1;
            const int qp = q * ru; // = ru

            newC.q_i = qp;
            newC.r_i = r_i;
            newC.n_ij.assign(qp, 0);
            newC.n_ijk.assign((size_t)qp * r_i, 0);

            auto * __restrict nij  = newC.n_ij.data();
            auto * __restrict nijk = newC.n_ijk.data();

            // 各サンプルについて j' と k を1回で決めてカウント
            for (int n = 0; n < N; ++n) {
                const int xu = ds.x(n, u);   // フラット配列アクセス
                const int k  = ds.x(n, v);
                const int j2 = xu;           // j=0 なので j2 = xu

                ++nij[j2];
                ++nijk[(size_t)j2 * r_i + k];
            }

            double after = scorer.nodeScore(newC);
            return after - nodeScoreNow[v];
        }

        // 親ありの場合
        const vector<int>& jindex = jcache.get(v); // Pa(v) に対応する j_index（キャッシュから取得）
        const int q   = curC.q_i;
        const int qp  = q * ru;

        newC.q_i = qp;
        newC.r_i = r_i;
        newC.n_ij.assign(qp, 0);
        newC.n_ijk.assign((size_t)qp * r_i, 0);

        auto * __restrict nij  = newC.n_ij.data();
        auto * __restrict nijk = newC.n_ijk.data();
        const int * __restrict jptr = jindex.data();

        // u の桁を「昇順の正しい位置」に挿入するための重みを求める。
        // 混合基数の桁順は g.parents(v)（= 昇順）に対応しており、
        // computeCountsForNode_full / JIndexCache / deltaRemove もこの順を前提に
        // している。ここで u を単に最下位桁として足すと（j2 = j*ru + xu）、
        // u より大きい添字の親が既にいる場合に桁順が正準順とずれ、以後の
        // REMOVE が別の親の桁をマージしてスコアが壊れる。
        //   right_u = u より添字が大きい既存の親の基数の総乗（= u の桁の重み）
        int right_u = 1;
        {
            const vector<int> Pa_v = g.parents(v);
            for (size_t t = 0; t < Pa_v.size(); ++t)
                if (Pa_v[t] > u) right_u *= ds.r[Pa_v[t]];
        }
        const int period_u = right_u * ru;

        // 各サンプルについて j' と k を1回で決めてカウント
        for (int n = 0; n < N; ++n) {
            const int j  = jptr[n];
            const int xu = ds.x(n, u);   // フラット配列アクセス
            const int k  = ds.x(n, v);
            // deltaRemove の写像 j -> (j/period)*right + (j%right) の逆写像
            const int j2 = (j / right_u) * period_u + xu * right_u + (j % right_u);

            ++nij[j2];
            ++nijk[(size_t)j2 * r_i + k];
        }

        double after = scorer.nodeScore(newC);
        return after - nodeScoreNow[v];
    }

    // REMOVE(u->v): **完全インクリメンタル**
    // - 現在の counts（Pa(v)）から、削除対象 u の「桁」を落として合算するだけ（O(q·r_i)）
    // - データの再走査無し
    double deltaRemove_andBuildNewCounts(int u, int v, Counts& newC) {
        const Counts& curC = nodeCounts[v]; // Pa(v) に基づく現在の counts

        // Pa(v) を取得して、u の位置（桁）と各基数を把握
        vector<int> Pa = g.parents(v);
        int pos = -1; vector<int> pa_r; pa_r.reserve(Pa.size());
        for (int t=0; t<(int)Pa.size(); ++t){
            pa_r.push_back(ds.r[Pa[t]]);
            if (Pa[t]==u) pos=t;
        }
        if (pos<0) throw runtime_error("deltaRemove: u is not a parent of v (logic error).");

        int r_u = ds.r[u];
        int r_i = curC.r_i;
        int q    = curC.q_i;
        int q2   = q / r_u; // u の桁を落とすと親配置数は r_u 倍減る

        // 右側（u より「下位桁」側）の基数の総乗
        int right = 1;
        for (int t=pos+1; t<(int)pa_r.size(); ++t) right *= pa_r[t];
        int period = right * r_u; // 1つ上の繰り返し周期

        newC.q_i = q2; newC.r_i = r_i;
        newC.n_ij.assign(q2, 0);
        newC.n_ijk.assign((size_t)q2*r_i, 0);

        // 旧インデックス j を新インデックス j' へ写像して合算
        // j' = floor(j / period) * right + (j % right)
        for (int j=0; j<q; ++j){
            int jp = (j / period) * right + (j % right);
            newC.n_ij[jp] += curC.n_ij[j];
            size_t base_old = (size_t)j * r_i;
            size_t base_new = (size_t)jp * r_i;
            for (int k=0;k<r_i;++k){
                newC.n_ijk[base_new + k] += curC.n_ijk[base_old + k];
            }
        }

        double after = scorer.nodeScore(newC);
        return after - nodeScoreNow[v];
    }

    // REVERSE(u->v): v 側は REMOVE（マージ）、u 側は ADD（分割）
    double deltaReverse_buildNewCounts(int u, int v, Counts& newC_v, Counts& newC_u, double& d_v, double& d_u){
        d_v = deltaRemove_andBuildNewCounts(u, v, newC_v); // v: 親 u を外す（マージ）
        d_u = deltaAdd_andBuildNewCounts(v, u, newC_u);    // u: 親に v を加える（分割）
        return d_v + d_u;
    }

    //==================== 実行（HC / Tabu） ====================

    tuple<DAG,double,int> run(bool use_tabu){
        double cur = totalNow;
        DAG bestG = g; double bestScore = cur;

        if (verbose) cerr << fixed << setprecision(6)
                          << "[start] score="<<cur<<" edges="<<g.edges().size()
                          << " mode="<<(use_tabu?"tabu":"greedy")<<"\n";

        const int D = ds.D;
        // タブー配列（属性タブー：直近操作の巻き戻しを禁止）
        vector<vector<int>> tabu_add_until(D, vector<int>(D,-1));
        vector<vector<int>> tabu_remove_until(D, vector<int>(D,-1));
        vector<vector<int>> tabu_reverse_until(D, vector<int>(D,-1));

        auto isTabu = [&](const Move& mv, int it)->bool{
            if (!use_tabu) return false;
            if (mv.type==Move::ADD)    return tabu_add_until[mv.u][mv.v] > it;
            if (mv.type==Move::REMOVE) return tabu_remove_until[mv.u][mv.v] > it;
            if (mv.type==Move::REVERSE)return tabu_reverse_until[mv.u][mv.v] > it;
            return false;
        };
        auto setTabuAfter = [&](const Move& mv, int it){
            if (!use_tabu) return;
            int until = it + tabu_tenure;
            if (mv.type==Move::ADD){
                // add の直後：remove と reverse（戻し）を禁止
                tabu_remove_until[mv.u][mv.v] = until;
                tabu_reverse_until[mv.u][mv.v] = until;
            } else if (mv.type==Move::REMOVE){
                // remove の直後：add（元に戻す）を禁止
                tabu_add_until[mv.u][mv.v] = until;
            } else if (mv.type==Move::REVERSE){
                // reverse の直後：逆向き reverse（巻き戻し）を禁止
                tabu_reverse_until[mv.v][mv.u] = until;
            }
        };

	int it=0;
        for (; it<max_iter; ++it){
            Move best; best.type=Move::NONE; best.delta=use_tabu?-1e300:0.0;
            Move bestNonTabu; bestNonTabu.type=Move::NONE; bestNonTabu.delta=-1e300;

            // 近傍列挙（ADD は候補親 K のみ、REMOVE/REVERSE は現辺）
            for (int v=0; v<D; ++v){
                // --- ADD 候補 ---
                if (topK>0){
                    for (int u: candParents[v]){
                        if (!addAllowed(u,v)) continue;
                        Counts newC;
                        double d = deltaAdd_andBuildNewCounts(u, v, newC);
                        Move mv{Move::ADD,u,v,d};
                        bool tabu = isTabu(mv, it);
                        bool asp = (cur + d > bestScore + 1e-12); // アスピレーション
                        if (!tabu) { if (d > bestNonTabu.delta) bestNonTabu = mv; }
                        if (!tabu || asp) { if (d > best.delta) best = mv; }
                    }
                } else {
                    for (int u=0; u<D; ++u) if (u!=v){
                        if (!addAllowed(u,v)) continue;
                        Counts newC;
                        double d = deltaAdd_andBuildNewCounts(u, v, newC);
                        Move mv{Move::ADD,u,v,d};
                        bool tabu = isTabu(mv, it);
                        bool asp = (cur + d > bestScore + 1e-12);
                        if (!tabu) { if (d > bestNonTabu.delta) bestNonTabu = mv; }
                        if (!tabu || asp) { if (d > best.delta) best = mv; }
                    }
                }

                // --- REMOVE 候補（現にある辺）---
                for (int u=0; u<D; ++u) if (g.adj[u][v]){
                    Counts newC;
                    double d = deltaRemove_andBuildNewCounts(u, v, newC);
                    Move mv{Move::REMOVE,u,v,d};
                    bool tabu = isTabu(mv, it);
                    bool asp = (cur + d > bestScore + 1e-12);
                    if (!tabu) { if (d > bestNonTabu.delta) bestNonTabu = mv; }
                    if (!tabu || asp) { if (d > best.delta) best = mv; }
                }

                // --- REVERSE 候補（現にある辺で逆向き追加が可能なもの）---
                for (int u=0; u<D; ++u) if (g.adj[u][v]){
                    if (!reverseAllowed(u,v)) continue;
                    Counts newCv, newCu; double dv=0.0, du=0.0;
                    double d = deltaReverse_buildNewCounts(u, v, newCv, newCu, dv, du);
                    Move mv{Move::REVERSE,u,v,d};
                    bool tabu = isTabu(mv, it);
                    bool asp = (cur + d > bestScore + 1e-12);
                    if (!tabu) { if (d > bestNonTabu.delta) bestNonTabu = mv; }
                    if (!tabu || asp) { if (d > best.delta) best = mv; }
                }
            }

            // タブーによる最良不可→非タブー最良を選ぶ
            Move chosen = best;
            if (chosen.type==Move::NONE && use_tabu && bestNonTabu.type!=Move::NONE) chosen = bestNonTabu;
            if (chosen.type==Move::NONE || (!use_tabu && chosen.delta<=1e-12)){
                if (verbose) cerr << "[stop] no improving move.\n";
                break;
            }

            // ===== 実適用（DAG / reach / j_index / counts / scores を整合させる） =====
            if (chosen.type==Move::ADD){
                Counts newC; double d = deltaAdd_andBuildNewCounts(chosen.u, chosen.v, newC);
                g.addEdge(chosen.u, chosen.v);
                reach.onAddEdge(g, chosen.u, chosen.v);
                jcache.onParentsChanged({chosen.v});          // v の親が変わったので j_index[v] を無効化
                nodeCounts[chosen.v] = move(newC);            // v の counts を差し替え
                nodeScoreNow[chosen.v] += d;                  // v のスコアだけ増分更新
                totalNow += d;                                 // 全体スコア更新

            } else if (chosen.type==Move::REMOVE){
                Counts newC; double d = deltaRemove_andBuildNewCounts(chosen.u, chosen.v, newC);
                g.removeEdge(chosen.u, chosen.v);
                reach.onRemoveEdge(g, chosen.u, chosen.v);
                jcache.onParentsChanged({chosen.v});
                nodeCounts[chosen.v] = move(newC);
                nodeScoreNow[chosen.v] += d;
                totalNow += d;

            } else if (chosen.type==Move::REVERSE){
                Counts newCv, newCu; double dv=0.0, du=0.0;
                double d = deltaReverse_buildNewCounts(chosen.u, chosen.v, newCv, newCu, dv, du);
                g.reverseEdge(chosen.u, chosen.v);
                reach.onReverseEdge(g, chosen.u, chosen.v);
                jcache.onParentsChanged({chosen.u, chosen.v});
                nodeCounts[chosen.v] = move(newCv);
                nodeCounts[chosen.u] = move(newCu);
                nodeScoreNow[chosen.v] += dv;
                nodeScoreNow[chosen.u] += du;
                totalNow += d;
            }

            setTabuAfter(chosen, it);

            if (totalNow > bestScore + 1e-12){ bestScore=totalNow; bestG=g; if (verbose) cerr<<"[*] new best "<<bestScore<<"\n"; }
            if (verbose){
                string t = (chosen.type==Move::ADD?"ADD": chosen.type==Move::REMOVE?"REM":"REV");
                cerr << "[it "<<it+1<<"] "<<t<<" "<<chosen.u<<"->"<<chosen.v
                     <<"  delta="<<chosen.delta<<"  cur="<<totalNow<<"\n";
            }
        }
        return {bestG, bestScore, it};
    }
};

//==================== 引数処理 / 初期構造ロード / MI 前処理 ====================

static ScoreType parseScore(const string& s){
    string t=s; for (auto& c:t) c=(char)tolower(c);
    if (t=="bic") return ScoreType::BIC;
    if (t=="k2")  return ScoreType::K2;
    if (t=="bdeu")return ScoreType::BDeu;
    throw runtime_error("Unknown score: "+s+" (use: bic|k2|bdeu)");
}

// 初期構造（エッジリスト）を読み込む
// ・タブ区切り（TSV）およびスペース区切りに対応
// ・コメント行 (#で始まる) は無視
static DAG loadInitEdges(int D, const string& path) {
    DAG g(D);
    if (path.empty()) return g;

    ifstream fin(path);
    if (!fin) throw runtime_error("Failed to open init edge list: " + path);

    string line;
    int line_no = 0;
    while (getline(fin, line)) {
        ++line_no;
        if (line.empty()) continue;
        if (line[0] == '#') continue; // コメント行スキップ

        // 区切り文字自動判定（タブ or スペース）
        char delim = (line.find('\t') != string::npos) ? '\t' : ' ';

        stringstream ss(line);
        string u_str, v_str;
        if (!getline(ss, u_str, delim)) continue;
        if (!getline(ss, v_str, delim)) continue;

        // 前後空白除去
        auto trim = [](string& s) {
            s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char c){return !isspace(c);} ));
            s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char c){return !isspace(c);} ).base(), s.end());
        };
        trim(u_str);
        trim(v_str);

        if (u_str.empty() || v_str.empty()) continue;

        int u = stoi(u_str);
        int v = stoi(v_str);
        if (u < 0 || u >= D || v < 0 || v >= D || u == v)
            throw runtime_error("Invalid edge at line " + to_string(line_no) + ": " + u_str + " " + v_str);

        if (!g.adj[u][v]) g.addEdge(u, v);
    }

    cerr << "[info] loaded init edges from " << path << " (" << g.edges().size() << " edges)\n";
    return g;
}

// ================================================
// 最終構造に基づく全ノード分のカウントを 1 つの TSV に出力
// - outfile: 出力パス（TSV）
// - ds     : データセット（基数・変数名に使用）
// - g      : 最終 DAG（親集合取得に使用）
// - counts : 各ノードの最終カウント（nodeCounts[v]）
// 形式：タブ区切り（TSV）
//   ・各ノードの前にメタ情報を # で出力（親の並び・基数・変数名など）
//   ・本体の行は 2 種類：
//        v <tab> j <tab> k <tab> n_ijk
//        v <tab> j <tab> * <tab> n_ij     （* は合計を表す）
//   ・n=0 の行は省略（ファイルサイズ節約）
// ================================================
void saveAllCountsTSV(const std::string& path, const Dataset& ds, const DAG& g) {
    std::ofstream fout(path);
    if (!fout) throw std::runtime_error("open failed: " + path);

    for (int v = 0; v < ds.D; ++v) {
        auto pa = g.parents(v);

        // メタ情報（親が空なら空行を出す）
        fout << "# --- node " << v << " ---\n";
        fout << "# node_name\t" << (v < (int)ds.var_names.size() ? ds.var_names[v] : ("X"+std::to_string(v))) << "\n";
        fout << "# parents_indices\t";
        for (size_t t=0; t<pa.size(); ++t) fout << (t? ",":"") << pa[t];
        fout << "\n# parents_names\t";
        for (size_t t=0; t<pa.size(); ++t) {
            int p = pa[t];
            fout << (t? ",":"") << (p < (int)ds.var_names.size() ? ds.var_names[p] : ("X"+std::to_string(p)));
        }
        fout << "\n# parents_cardinalities\t";
        for (size_t t=0; t<pa.size(); ++t) fout << (t? ",":"") << ds.r[pa[t]];
        fout << "\n# child_cardinality\t" << ds.r[v] << "\n";

        // 再度カウント（q_i は親配置数）
        Counts C = computeCountsForNode_full(v, pa, ds);
        const int r_i = C.r_i;         // 子の取りうる値の数
        int q_i = C.q_i;               // 親配置数
        if (pa.empty()) q_i = 1;       // 念のため明示

        // j は 0..q_i-1 で回すのが正しい
        for (int j = 0; j < q_i; ++j) {
            long long nij = C.n_ij[j];

            // 各 k 行（ゼロは省略
            for (int k = 0; k < r_i; ++k) {
                long long nijk = C.n_ijk[(size_t)j*r_i + k];
                if (nijk == 0) continue;  // 省略したくない場合は消す
                fout << v << "\t" << j << "\t" << k << "\t" << nijk << "\n";
            }
            // 合計行は必ず出す
            fout << v << "\t" << j << "\t*\t" << nij << "\n";
        }
    }
}

// ============ all_counts.tsv ロード ============
// all_counts.tsv（単一ファイル）から、各ノード v の Counts を復元する。
// フォーマット：データ行 "v \t j \t k \t n" （k='*' のとき n_ij）
static vector<Counts> loadAllCountsTSV(const string& path, int expected_D, vector<int>* out_child_r = nullptr) {
    ifstream fin(path);
    if (!fin) throw runtime_error("Failed to open all-counts file: " + path);

    // 一時的に「ノードごと」に (j,k,n_ijk) と (j, n_ij) を保持
    struct KEntry { int j; int k; long long n; };
    struct JEntry { int j; long long n; };
    vector<vector<KEntry>> tmpK; // tmpK[v] に (j,k,n)
    vector<vector<JEntry>> tmpJ; // tmpJ[v] に (j,n)
    tmpK.resize(expected_D);
    tmpJ.resize(expected_D);

    // サイズ推定用
    vector<int> maxJ(expected_D, -1), maxK(expected_D, -1);

    string line;
    long long data_lines = 0;
    while (getline(fin, line)) {
        if (line.empty() || line[0]=='#') continue;

        // v \t j \t k \t n
        // k は整数か '*'（合計）を取る
        string vstr, jstr, kstr, nstr;
        {
            stringstream ss(line);
            if (!getline(ss, vstr, '\t')) continue;
            if (!getline(ss, jstr, '\t')) continue;
            if (!getline(ss, kstr, '\t')) continue;
            if (!getline(ss, nstr, '\t')) continue;
        }
        // trim 簡易
        auto trim = [](string& s){
            s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char c){return !isspace(c);} ));
            s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char c){return !isspace(c);} ).base(), s.end());
        };
        trim(vstr); trim(jstr); trim(kstr); trim(nstr);

        int v = stoi(vstr);
        int j = stoi(jstr);
        if (v<0 || v>=expected_D) throw runtime_error("counts file: node index out of range: "+to_string(v));
        long long n = stoll(nstr);

        if (kstr == "*" || kstr == "'*'") {
            // n_ij
            tmpJ[v].push_back({j, n});
            maxJ[v] = max(maxJ[v], j);
        } else {
            int k = stoi(kstr);
            tmpK[v].push_back({j, k, n});
            maxJ[v] = max(maxJ[v], j);
            maxK[v] = max(maxK[v], k);
        }
        ++data_lines;
    }
    if (data_lines==0) throw runtime_error("counts file has no data lines: " + path);

    // 復元（配列確保→詰め込み）
    vector<Counts> C(expected_D);
    vector<int> child_r_detected(expected_D, 0);

    for (int v=0; v<expected_D; ++v){
        int q = (maxJ[v] >= 0 ? (maxJ[v]+1) : 1); // 観測なしなら q=1 とみなす
        int r = (maxK[v] >= 0 ? (maxK[v]+1) : 1); // 観測なしなら r=1 とみなす

        C[v].q_i = q;
        C[v].r_i = r;
        C[v].n_ij.assign(q, 0);
        C[v].n_ijk.assign((size_t)q * r, 0);

        for (auto &e : tmpJ[v]) {
            if (e.j < 0 || e.j >= q) throw runtime_error("n_ij j out of range at node "+to_string(v));
            C[v].n_ij[e.j] += e.n;
        }
        for (auto &e : tmpK[v]) {
            if (e.j < 0 || e.j >= q) throw runtime_error("n_ijk j out of range at node "+to_string(v));
            if (e.k < 0 || e.k >= r) throw runtime_error("n_ijk k out of range at node "+to_string(v));
            C[v].n_ijk[(size_t)e.j * r + e.k] += e.n;
        }
        child_r_detected[v] = r;
    }

    if (out_child_r) *out_child_r = move(child_r_detected);
    return C;
}

// ============ 新データセットに対する対数尤度 ============
// alpha_ij: 平滑化（Dirichlet 事前）パラメータ a
//   P(x_i=k | pa=j) = (n_ijk + a/r_i) / (n_ij + a)
//   a=0 なら MLE（ゼロ割・確率0→ -inf に注意）
static double computeLogLikelihoodOnDataset(const Dataset& ds_new,
                                            const DAG& g,
                                            const vector<Counts>& C,
                                            double alpha_ij,
                                            double* out_avg_per_sample = nullptr,
                                            double* out_avg_per_var = nullptr,
                                            long long* out_zero_hits = nullptr)
{
    if (g.D != ds_new.D) throw runtime_error("D mismatch: graph vs dataset");
    if ((int)C.size() != g.D) throw runtime_error("Counts vector size != D");

    const int N = ds_new.N;
    const int D = ds_new.D;

    // 親集合は g.parents(v) の順序を使用（counts の j もこの順序で作っている前提）
    double LL = 0.0;
    long long zero_hits = 0;

    // 事前のための定数（r_i は C[v].r_i を使う）
    for (int n=0; n<N; ++n){
        for (int v=0; v<D; ++v){
            const auto parents = g.parents(v);
            int j = 0;
            if (!parents.empty()){
                // 混合基数（右端の親が最下位桁）で j を計算
                int mult = 1;
                for (int t=(int)parents.size()-1; t>=0; --t){
                    int p = parents[t];
                    j += ds_new.x(n,p) * mult;
                    mult *= ds_new.r[p];
                }
            }
            int k = ds_new.x(n,v);
            const Counts& Cv = C[v];
            if (j<0 || j>=Cv.q_i) {
                // counts に存在しない親配置（例えば子や親の基数が異なる等）
                // → smoothing がゼロなら -inf、a>0 なら一様事前で評価
                if (alpha_ij<=0.0) return -INFINITY;
                double prob = (1.0 / (double)Cv.r_i); // n_ij=0 前提で (0+a/r)/(0+a)
                LL += log(prob);
                continue;
            }

            long long nij = Cv.n_ij[j];
            long long nijk = (k>=0 && k<Cv.r_i) ? Cv.n_ijk[(size_t)j*Cv.r_i + k] : 0;

            double prob;
            if (alpha_ij > 0.0){
                prob = ( (double)nijk + alpha_ij / (double)Cv.r_i ) / ( (double)nij + alpha_ij );
            } else {
                if (nij == 0 || nijk == 0) {
                    ++zero_hits;
                    return -INFINITY; // MLE で確率0に当たった
                }
                prob = (double)nijk / (double)nij;
            }
            LL += log(prob);
        }
    }

    if (out_avg_per_sample) *out_avg_per_sample = LL / (double)N;
    if (out_avg_per_var)    *out_avg_per_var    = LL / (double)(N * D);
    if (out_zero_hits)      *out_zero_hits      = zero_hits;
    return LL;
};

struct DeltaStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double stdev = std::numeric_limits<double>::quiet_NaN();
    long long finite_count = 0;  // 有限値として集計できたサンプル数
};

// 親集合に対する j_index（各サンプルの親配置インデックス）を構築
static vector<int> buildJIndexForParents(const Dataset& ds,
                                         const vector<int>& parents)
{
    vector<int> jidx(ds.N, 0);
    if (parents.empty()) return jidx;
    // 右端の親が最下位桁になる混合基数
    std::vector<int> radix;
    build_mixed_radix(parents, ds.r, radix);
    for (int n=0; n<ds.N; ++n) {
	const int j = mixed_radix_index_row(ds, n, parents, radix);
        jidx[n]=j;
    }
    return jidx;
}


// 親集合 before/after のみが異なるノード v について、
// 新データ ds_new 上の「サンプルごとの logL 差分 (after - before)」の
// 平均・標準偏差を計算して返す。
// alpha_ij>0 なら (n_ijk + a/r_i)/(n_ij + a) を用いた平滑化付き確率、0ならMLE。
// ※ MLEで0確率に当たるサンプルは統計から除外（finite のみ集計）。
DeltaStats perSampleDeltaLogLStats(const Dataset& ds_new,
                                   int v,
                                   const std::vector<int>& parents_before,
                                   const std::vector<int>& parents_after,
                                   double alpha_ij)
{
    // before/after の counts と j_index を用意
    Counts C_before = computeCountsForNode_full(v, parents_before, ds_new);
    Counts C_after  = computeCountsForNode_full(v, parents_after , ds_new);

    auto buildJ = [&](const std::vector<int>& pa){
        return buildJIndexForParents(ds_new, pa);
    };
    std::vector<int> j_before = buildJ(parents_before);
    std::vector<int> j_after  = buildJ(parents_after);

    auto prob_from_counts = [&](const Counts& C, int j, int k)->double{
        if (j<0 || j>=C.q_i || k<0 || k>=C.r_i) return 0.0;
        double nij  = (double)C.n_ij[j];
        double nijk = (double)C.n_ijk[(size_t)j*C.r_i + k];
        if (alpha_ij > 0.0) {
            return (nijk + alpha_ij/(double)C.r_i) / (nij + alpha_ij);
        } else {
            if (nij<=0.0 || nijk<=0.0) return 0.0;
            return nijk / nij;
        }
    };

    long long cnt = 0;
    long double sum = 0.0L, sum2 = 0.0L;

    for (int n=0; n<ds_new.N; ++n) {
        int kx = ds_new.x(n,v);
        double p_b = prob_from_counts(C_before, j_before[n], kx);
        double p_a = prob_from_counts(C_after , j_after [n], kx);
        double lb = (p_b>0.0 ? std::log(p_b) : -INFINITY);
        double la = (p_a>0.0 ? std::log(p_a) : -INFINITY);
        double d  = la - lb; // after - before
        if (std::isfinite(d)) {
            ++cnt;
            long double dl = (long double)d;
            sum  += dl;
            sum2 += dl*dl;
        }
    }

    DeltaStats st;
    st.finite_count = cnt;
    if (cnt > 0) {
        long double mu = sum / (long double)cnt;
        st.mean = (double)mu;
        if (cnt > 1) {
            long double var = (sum2 - (long double)cnt*mu*mu) / (long double)(cnt-1);
            st.stdev = (var>0 ? std::sqrt((double)var) : 0.0);
        } else {
            st.stdev = 0.0;
        }
    }
    return st;
}

// ============================================================
// 各エッジ (u->v) を1本ずつ削除したときのスコア変化を計算
// - データセット ds_new : 新しい評価データ
// - g_base : 学習済み構造（init_edges.tsv）
// - counts : all_counts.tsv から読み込んだ n_ijk / n_ij
// - alpha_ij : 平滑化 (MLE時=0.0)
// - ess : BDeu スコア用の等価サンプルサイズ
// - outfile : 結果TSVファイル (u, v, ΔlogL, ΔBIC, ΔK2, ΔBDeu)
//
// 計算の流れ：
//   1. 元スコアを計算（logL, BIC, K2, BDeu）
//   2. 各エッジ u->v を削除した DAG を作成
//   3. そのノード v のみ再スコア（親集合からuを除く）
//   4. 各スコア差 Δ = (after - before) を出力
// ============================================================
void computeEdgeImportanceScores(const Dataset& ds_new,
                                 const DAG& g_base,
                                 const std::vector<Counts>& counts,
                                 double alpha_ij,
                                 double ess,
                                 const std::string& outfile)
{
    std::ofstream fout(outfile);
    if (!fout) throw std::runtime_error("Failed to open edge-importance file: " + outfile);
    fout << "u\tv\tΔlogL\tΔBIC\tΔK2\tΔBDeu\n";

    const int D = ds_new.D;
    const int N = ds_new.N;
    const double logN = std::log((double)std::max(1, N));

    // ========== スコア計算関数群 ==========
    auto nodeLogLikelihood = [&](int v, const DAG& g) -> double {
        const auto pa = g.parents(v);
        Counts C = computeCountsForNode_full(v, pa, ds_new);
        double ll = 0.0;
        for (int j = 0; j < C.q_i; ++j) {
            double nij = (double)C.n_ij[j];
            if (nij <= 0) continue;
            for (int k = 0; k < C.r_i; ++k) {
                long long nijk = C.n_ijk[j*C.r_i + k];
                if (nijk == 0) continue;
                ll += nijk * (std::log((double)nijk) - std::log(nij));
            }
        }
        return ll;
    };

    auto nodeBIC = [&](int v, const DAG& g) -> double {
        const auto pa = g.parents(v);
        Counts C = computeCountsForNode_full(v, pa, ds_new);
        double ll = 0.0;
        for (int j = 0; j < C.q_i; ++j) {
            double nij = (double)C.n_ij[j];
            if (nij <= 0) continue;
            for (int k = 0; k < C.r_i; ++k) {
                long long nijk = C.n_ijk[j*C.r_i + k];
                if (nijk == 0) continue;
                ll += nijk * (std::log((double)nijk) - std::log(nij));
            }
        }
        int d = (C.r_i - 1) * C.q_i;
        double pen = 0.5 * d * logN;
        return ll - pen;
    };

    auto nodeK2 = [&](int v, const DAG& g) -> double {
        const auto pa = g.parents(v);
        Counts C = computeCountsForNode_full(v, pa, ds_new);
        double s = 0.0;
        for (int j = 0; j < C.q_i; ++j) {
            double nij = (double)C.n_ij[j];
            s += std::lgamma((double)C.r_i) - std::lgamma(nij + (double)C.r_i);
            for (int k = 0; k < C.r_i; ++k) {
                double nijk = (double)C.n_ijk[j*C.r_i + k];
                s += std::lgamma(nijk + 1.0);
            }
        }
        return s;
    };

    auto nodeBDeu = [&](int v, const DAG& g) -> double {
        const auto pa = g.parents(v);
        Counts C = computeCountsForNode_full(v, pa, ds_new);
        double s = 0.0;
        if (C.q_i == 0) return -INFINITY;
        double alpha_ij_local = ess / (double)C.q_i;
        double alpha_ijk_base = alpha_ij_local / (double)C.r_i;
        for (int j = 0; j < C.q_i; ++j) {
            double nij = (double)C.n_ij[j];
            s += std::lgamma(alpha_ij_local) - std::lgamma(nij + alpha_ij_local);
            for (int k = 0; k < C.r_i; ++k) {
                double nijk = (double)C.n_ijk[j*C.r_i + k];
                s += std::lgamma(nijk + alpha_ijk_base) - std::lgamma(alpha_ijk_base);
            }
        }
        return s;
    };

    // ========== 元スコア計算 ==========
    double baseLL = 0.0, baseBIC = 0.0, baseK2 = 0.0, baseBDeu = 0.0;
    for (int v = 0; v < D; ++v) {
        baseLL  += nodeLogLikelihood(v, g_base);
        baseBIC += nodeBIC(v, g_base);
        baseK2  += nodeK2(v, g_base);
        baseBDeu+= nodeBDeu(v, g_base);
    }

    // ========== 各エッジ削除時のスコア差 ==========
    for (int u = 0; u < D; ++u) {
        for (int v = 0; v < D; ++v) {
            if (!g_base.hasEdge(u, v)) continue;

            DAG g_mod = g_base;
            g_mod.removeEdge(u, v);

            // 対象ノード v のみ再スコア
            double ll_new   = nodeLogLikelihood(v, g_mod);
            double bic_new  = nodeBIC(v, g_mod);
            double k2_new   = nodeK2(v, g_mod);
            double bdeu_new = nodeBDeu(v, g_mod);

            // 元スコアとの差分（Δ）
            double deltaLL   = ll_new   - nodeLogLikelihood(v, g_base);
            double deltaBIC  = bic_new  - nodeBIC(v, g_base);
            double deltaK2   = k2_new   - nodeK2(v, g_base);
            double deltaBDeu = bdeu_new - nodeBDeu(v, g_base);

            auto pa_before = g_base.parents(v);
            std::vector<int> pa_after; pa_after.reserve(pa_before.size());
            for (int p : pa_before) if (p!=u) pa_after.push_back(p);
            DeltaStats st = perSampleDeltaLogLStats(ds_new, v, pa_before, pa_after, alpha_ij);

            fout << u << "\t" << v << "\t"
                 << deltaLL   << "\t"
                 << deltaBIC  << "\t"
                 << deltaK2   << "\t"
                 << deltaBDeu << "\t"
                 << std::setprecision(12) << st.mean << '\t'
                 << std::setprecision(12) << st.stdev << '\n';
        }
    }

    fout.close();
    std::cerr << "[info] wrote edge importance results to " << outfile << "\n";
}
// ================================================================
// ブートストラップ学習を B 回まわして、エッジ出現回数を TSV で保存
//
// - 入力: 元データ ds（ヘッダ/基数付き）, スコア種別, ESS(BDeu用),
//         初期構造ファイルパス（空文字なら空DAGで開始OK）,
//         既存HCのパラメタ（tabu/iters/max-parents/max-children/topK 等）,
//         reach モード, jindex キャッシュ, ブートストラップ回数/seed,
//         MI の sample/budget 制御（各バッチで再計算）
//
// - 出力: save_path に TSV で保存（u v count prob）
//         include_zero==true のとき、全 (u!=v) の 0 含む一覧を出力（大規模Dでは注意）
//
// 注意: メモリ節約のためカウントは unordered_map で疎管理
// ================================================================
void runBootstrapStructureCounts(const Dataset& ds,
                                 ScoreType score_type,
                                 double ess_for_bdeu,
                                 const std::string& init_path,
                                 int tabu_tenure,
                                 int iters,
                                 int max_parents,
                                 int max_children,
				 CandMetric cand_metric,
                                 int topK,
                                 int mi_sample,    // 0: 全行
                                 int mi_budget,    // 0: 全変数
				 double mi_threshold,
				 double chi2_p_threshold,
                                 Reachability::Mode reach_mode,
                                 int jindex_cache_cap,
                                 int B,            // ブートストラップ反復回数
                                 uint64_t seed,
                                 const std::string& save_path)
{
    if (B <= 0) throw std::runtime_error("bootstrap B must be > 0");
    std::mt19937_64 rng(seed);

    // 疎カウント：エッジ (u->v) の出現回数
    std::unordered_map<uint64_t, uint32_t> edge_counts;
    edge_counts.reserve((size_t)ds.D * 8); // 適当に初期予約（出現が疎な想定）

    for (int b = 0; b < B; ++b) {
        // 1) ブートストラップ標本の生成
        Dataset ds_b = bootstrapResampleDataset(ds, rng);

        // 2) 初期構造のロード
        DAG init = loadInitEdges(ds_b.D, init_path);

        // 3) 候補親Kの前処理（サンプル/バジェットは引数で制御）
        vector<vector<int>> topKlist;
        if (topK >0 || mi_threshold>0.0 || chi2_p_threshold<1.0){
            vector<int> rows;
            if (mi_sample > 0 && mi_sample < ds_b.N) {
                rows.resize(ds_b.N); iota(rows.begin(), rows.end(), 0);
                shuffle(rows.begin(), rows.end(), rng);
                rows.resize(mi_sample);
                sort(rows.begin(), rows.end());
            } else {
                rows.resize(ds_b.N); iota(rows.begin(), rows.end(), 0);
            }
            int budget = (mi_budget > 0 ? mi_budget : (ds_b.D - 1));
            //topKlist = MICandidates::compute(ds_b, topK, budget, rows, rng, mi_threshold);
	    topKlist = AssocCandidates::compute(
                ds, topK, budget, rows, rng, cand_metric, mi_threshold, chi2_p_threshold
            );
        }
        
        // 4) 探索器の構築 & 実行（学習本体）
        HillClimber hc(ds_b, score_type, ess_for_bdeu, init, reach_mode, jindex_cache_cap);
        hc.max_iter = iters;
        hc.tabu_tenure = tabu_tenure;
        hc.max_parents = max_parents;
        hc.max_children = max_children;
        if (topK > 0) { hc.candParents = std::move(topKlist); hc.topK = topK; }
        hc.verbose = false; // ブートストラップ時は静かに

        auto [g_learned, score, it] = hc.run(/*use_tabu=*/(tabu_tenure > 0));

        // 5) 出現エッジを蓄積
        for (auto &e : g_learned.edges()) {
            uint64_t key = edgeKeyUV((uint32_t)e.first, (uint32_t)e.second);
            auto it = edge_counts.find(key);
            if (it == edge_counts.end()) edge_counts.emplace(key, 1u);
            else ++(it->second);
        }

        // 進捗表示（必要なら）
        if ( (b % 10) == 0 ) {
            std::cerr << "[bootstrap] " << (b+1) << "/" << B
                      << " edges_seen=" << edge_counts.size() << "\n";
        }
    }

    // 6) 出力：TSV（u, v, count, prob）
    //    include_zero_edges が true なら（u!=v）の全組み合わせで 0 も出力（巨大Dでは注意）
    namespace fs = std::filesystem;
    fs::path base_path(save_path);
    std::string stem = base_path.stem().string();   // 例: "boot_edges"
    std::string ext  = base_path.extension().string(); // 例: ".tsv"
    fs::path parent  = base_path.parent_path();

    std::ostringstream seed_str;
    seed_str << std::setfill('0') << std::setw(4) << (seed % 10000);
    // ファイル名: <stem>_seed<seed><ext>
    std::ostringstream oss;
    oss << stem << "_seed" << seed_str.str() << ext;
    fs::path out_path = parent / oss.str();

    std::ofstream fout(out_path);

    if (!fout) throw std::runtime_error("Failed to open --save-bootstrap-counts: " + save_path);

    fout << "u\tv\tcount\tprob\n";

    // 出現したエッジのみ
    // キーを (u,v) に分解して出力
    for (auto &kv : edge_counts) {
        uint64_t key = kv.first;
        uint32_t c = kv.second;
        uint32_t u = (uint32_t)(key >> 32);
        uint32_t v = (uint32_t)(key & 0xffffffffULL);
        double p = (double)c / (double)B;
        fout << u << "\t" << v << "\t" << c << "\t" << std::setprecision(12) << p << "\n";
    }

    fout.close();
    std::cerr << "[info] wrote bootstrap edge counts to " << save_path
              << " (B=" << B << ", unique_edges=" << edge_counts.size() << ")\n";
}

static void print_help() {
    const char* env_lang = std::getenv("LANG");
    std::string lang = (env_lang ? std::string(env_lang) : "ja");

    if (lang.substr(0,2) == "en") {
	    std::cout <<
R"(Bayesian Network Structure Learning Tool (Hill-Climb + Tabu + Bootstrap)
---------------------------------------------------------------------
Usage:
  fast_bn --input <data.tsv> --score <bic|k2|bdeu> [options]

General options:
  -h, --help
      Show this help message and exit.
  --input <file>
      Input dataset file (CSV/TSV with header, separator auto-detected).
  --score <bic|k2|bdeu>
      Scoring function (default=bic).
  --init <file>
      Initial network structure (tab-separated u v).
  --iters <n>
      Max iterations for Hill-Climb (default=1000).
  --tabu <len>
      Tabu search tenure. 0 disables tabu (default=0).
  --max-parents <k>
      Maximum number of parents per node (default=5).
  --max-children <k>
      Maximum number of children per node (default=unlimited).
  --ess <val>
      Equivalent sample size for BDeu (default=1.0).
  --alpha <a>
      Smoothing coefficient for likelihood (default=1.0).
  --verbose
      Enable verbose output (default=on).
  --quiet
      Suppress detailed logging.

Candidate parent selection:
  --cand-metric mi|chi2
      Measure for candidate parent selection: mi=Mutual Information, chi2=Chi-square (default=mi).
  --topk <K>
      Keep top-K strongest candidates (0=no limit, default=50).
  --mi-threshold <t>
      Keep variables with MI ≥ t (in nats, default=0.0).
  --chi2-p-threshold <p>
      Keep variables with p ≤ threshold (default=0.05).
  --mi-sample <n>
      Number of samples for MI/Chi2 calculation (0=use all, default=0).
  --mi-budget <n>
      Candidate variable limit per node (0=all, default=0).
  --reach <dense|lazy>
      Reachability check mode (default=lazy).
  --jindex-cache <cap>
      Parent configuration cache size (default=32).

Bootstrap mode:
  --bootstrap <B>
      Run bootstrap sampling B times (default=0).
  --seed <s>
      Random seed (default=2025).
  --save-bootstrap-counts <path>
      Output base filename; "_seed####.tsv" will be appended automatically.
  --bootstrap-include-zero
      Include edges that never appeared (default=off).

Edge importance mode:
  --edge-importance
      Evaluate score changes by removing each edge.
  --score-dataset <file>
      New dataset for evaluation.
  --init <file>
      Trained network structure (tab-separated u v).
  --counts <file>
      Trained all_counts.tsv file.
  --save-edge-importance <path>
      Output file for edge importance (TSV).

Output & runtime:
  --verbose
      Enable detailed logging (default=on).
  --quiet
      Suppress detailed logging.
---------------------------------------------------------------------
Example:
./fast_bn --input data.tsv --score bic --iters 2000 --tabu 20 --max-parent 3
./fast_bn --input gene.tsv --score bdeu --ess 10 --iters 3000

Example: bootstrap
./fast_bn --input data.tsv --score bic \
  --bootstrap 100 \
  --save-bootstrap-counts results/boot_edges.tsv

 
Example: edge importance
./fast_bn --input test.tsv --score bic \
  --edge-importance \
  --init init_edges.tsv \
  --counts all_counts.tsv \
  --save-edge-importance edge_imp.tsv
---------------------------------------------------------------------
	)";
    }else{
        std::cout <<
R"(Bayesian Network Structure Learning Tool (Hill-Climb + Tabu + Bootstrap)
---------------------------------------------------------------------
Usage:
  fast_bn --input <data.tsv> --score <score> [options]


General options:
  -h, --help
      このヘルプを表示して終了します。
  --input <file>
      入力データファイル（CSV/TSV, ヘッダ付き, 区切り自動判定）。
  --score <bic|k2|bdeu>
      スコア関数を指定します（default=bic）。
  --init <file>
      初期構造ファイル (タブ区切り u v)。
  --iters <n>
      Hill-Climb の最大反復数（default=1000）。
  --tabu <len>
      Tabu サーチの禁制期間。0で無効（default=0）。
  --max-parents <k>
      各ノードの最大親数（default=5）。
  --max-children <k>
      各ノードの最大子数（default=無制限）。
  --ess <val>
      BDeu スコアの等価サンプルサイズ（default=1.0）。
  --alpha <a>
      スムージング係数（default=1.0）。
  --verbose
      詳細ログを出力（default=有効）。
  --quiet
      詳細ログを抑制。

Candidate parent selection:
  --cand-metric mi|chi2
      候補親の関連度指標。mi=相互情報量, chi2=カイ二乗p値（default=mi）。
  --topk <K>
      候補親を上位K個に制限。0で制限なし（default=50）。
  --mi-threshold <t>
      MI >= t (nats単位) のみ採用（default=0.0）。
  --chi2-p-threshold <p>
      p <= p_threshold のみ採用（default=0.05）。
  --mi-sample <n>
      MI/chi2 計算時のサンプル数（0=全行, default=0）。
  --mi-budget <n>
      MI/chi2 計算時の候補数上限（0=全変数, default=0）。
  --reach <dense|lazy>
      到達可能性チェックモード（default=lazy）。
  --jindex-cache <cap>
      親配置キャッシュサイズ（default=32）。

Bootstrap mode:
  --bootstrap <B>
      ブートストラップ試行回数。B>0で実行モード（default=0）。
  --seed <s>
      乱数シード値（default=2025）。
  --save-bootstrap-counts <path>
      出力ファイルベース名。自動で "_seed####.tsv" が付与されます。
  --bootstrap-include-zero
      出現しなかったエッジも出力（default=off）。

Edge importance mode:
  --edge-importance
      エッジ除去によるスコア変化を評価。
  --score-dataset <file>
      評価用データセット（新しいCSV/TSV）。
  --init <file>
      学習済み構造（タブ区切り u v）。
  --counts <file>
      学習済み all_counts.tsv。
  --save-edge-importance <path>
      エッジ重要度を出力するTSVファイル。

Output & runtime:
  --save-bootstrap-counts <path>
      ブートストラップの出力ファイル（ベース名）。
  --seed <n>
      ブートストラップの乱数シード（default=2025）。
  --bootstrap-include-zero
      出現しないエッジも出力。
  --verbose
      詳細ログ出力（default=有効）。
  --quiet
      詳細ログを抑制。

---------------------------------------------------------------------
例:
./fast_bn --input data.tsv --score bic --iters 2000 --tabu 20 --max-parent 3
./fast_bn --input gene.tsv --score bdeu --ess 10 --iters 3000

例: bootstrap
./fast_bn --input data.tsv --score bic \
  --bootstrap 100 \
  --save-bootstrap-counts results/boot_edges.tsv

 
例: エッジ貢献
./fast_bn --input test.tsv --score bic \
  --edge-importance \
  --init init_edges.tsv \
  --counts all_counts.tsv \
  --save-edge-importance edge_imp.tsv
---------------------------------------------------------------------
)";
    }
}
int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc == 1){
        cerr << "Usage:\n"
             << "  "<<argv[0]<<" --input <csv_path> --score <score(bic|k2|bdeu)>\n"
             << "    [--ess E] [--init init_edges.txt]\n"
             << "    [--tabu T] [--iters N] [--quiet]\n"
             << "    [--max-parents M] [--max-children M]\n"
             << "    [--topk K] [--mi-sample S] [--mi-budget B]\n"
             << "    [--reach dense|lazy] [--jindex-cache C]\n";
        return 1;
    }

    std::string input_path;
    std::string score_name;
    ScoreType sc = ScoreType::BIC;

    // --- help オプション ---
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help();
            return 0;
        }
    }

    // ===== 既定値（D=10k〜30k を想定した省メモリ寄り） =====
    double ess = 1.0;
    string init_path;
    int tabu_tenure = 20;       // 0でタブー無効（純HC）
    int iters = 5000;
    bool verbose = true;
    int max_parents = 3, max_children = 16;
    int K = 50;                 // 候補親K
    int mi_sample = 8000;       // MI 行サンプル数（0で全行）
    int mi_budget = 1500;       // 1ノードあたり評価する相手変数数（0で全変数）
    string reach_mode_str = "lazy"; // 省メモリ既定
    int jindex_cache_cap = 1024;  // j_index LRU 容量（ノード数）
    string save_init_path;        // インデックスで保存
    string save_init_names_path;  // 変数名で保存（任意）
    string save_counts_path;   // カウントを保存する単一ファイル
    string score_dataset_path;  // 尤度を計算する新データセット
    string counts_in_path;      // all_counts.tsv の入力
    double alpha_ij = 0.0;      // 平滑化パラメータ（既定0=MLE）
    bool edge_importance_mode=false;
    string save_edge_importance_path;
    //
    int bootstrap_B = 0;                 // >0 でブートストラップ実行
    uint64_t seed = 2025;  // 既定シード
    std::string save_bootstrap_counts;   // 出力TSV
    bool bootstrap_include_zero = false; // 0カウントも出すか
    std::string cand_metric_str = "mi";   // "mi" | "chi2"
    double mi_threshold = 0.0;           // MIの閾値：単位は ナチュラル対数（nats）
    double chi2_p_threshold = 1.0;       // 採用は p <= しきい値

    //
    //
    // 追加オプションのパース
    for (int i=1;i<argc;++i){
        string a=argv[i];
        auto need=[&](const char* f){return a==f && i+1<argc;};
        if (need("--input")) input_path = argv[++i];
	else if (need("--score")) score_name = argv[++i];
	else if (need("--ess")) ess = stod(argv[++i]);
        else if (need("--init")) init_path = argv[++i];
        else if (need("--tabu")) tabu_tenure = stoi(argv[++i]);
        else if (need("--iters")) iters = stoi(argv[++i]);
        else if (need("--max-parents")) max_parents = stoi(argv[++i]);
        else if (need("--max-children")) max_children = stoi(argv[++i]);
        else if (need("--topk")) K = stoi(argv[++i]);
        else if (need("--mi-sample")) mi_sample = stoi(argv[++i]);
        else if (need("--mi-budget")) mi_budget = stoi(argv[++i]);
	else if (need("--mi-threshold")) mi_threshold = stod(argv[++i]);
	else if (need("--cand-metric"))        cand_metric_str = argv[++i];  // "mi" or "chi2"
        else if (need("--chi2-p-threshold"))   chi2_p_threshold = stod(argv[++i]);

        else if (need("--reach")) reach_mode_str = argv[++i];
        else if (need("--jindex-cache")) jindex_cache_cap = stoi(argv[++i]);
        else if (need("--save"))        save_init_path = argv[++i];
        else if (need("--save-names"))  save_init_names_path = argv[++i];
        else if (need("--save-counts")) save_counts_path = argv[++i];
        else if (need("--score-dataset"))  score_dataset_path = argv[++i];
        else if (need("--counts"))         counts_in_path = argv[++i];
        else if (need("--alpha"))          alpha_ij = stod(argv[++i]);
        else if (need("--edge-importance")) edge_importance_mode = true;
        else if (need("--save-edge-importance")) save_edge_importance_path = argv[++i];
        else if (need("--bootstrap"))             bootstrap_B = stoi(argv[++i]);
        else if (need("--seed"))        seed = stoull(argv[++i]);
        else if (need("--save-bootstrap-counts")) save_bootstrap_counts = argv[++i];
        else if (a=="--verbose") verbose=true;   // 既定で有効（明示指定も受け付ける）
        else if (a=="--quiet") verbose=false;
        else { cerr<<"Unknown option: "<<a<<"\n"; return 1; }
    }
    if(!score_name.empty()){
    	sc=parseScore(score_name);
    }

    if (sc!=ScoreType::BDeu && ess!=1.0) cerr << "[warn] --ess is ignored for non-BDeu.\n";
    using clock = std::chrono::steady_clock;
    auto t_start = clock::now();
    
    // metric 解釈
    CandMetric cand_metric = CandMetric::MI;
    if (cand_metric_str == "chi2" || cand_metric_str == "CHI2") cand_metric = CandMetric::CHI2;
    // ===== ブートストラップ・モード（学習は行わず、B回の再学習＋集計のみ） =====
    if (bootstrap_B > 0) {
        if (save_bootstrap_counts.empty()) {
            std::cerr << "Error: --save-bootstrap-counts <path> is required with --bootstrap\n";
            return 2;
        }
        // 入力データ（CSV/TSV・ヘッダあり）
        Dataset ds = Dataset::fromCSV(input_path);

        // 既存のオプション（score, ess, init, tabu, iters, max-parents/children, topK, mi_sample, mi_budget, reach_mode, jindex_cache_cap）
        runBootstrapStructureCounts(
          ds, sc, ess, init_path,
          tabu_tenure, iters,
          max_parents, max_children,
	  cand_metric,
          K, mi_sample, mi_budget, mi_threshold, chi2_p_threshold,
          (reach_mode_str=="dense"? Reachability::DENSE : Reachability::LAZY),
          jindex_cache_cap,
          bootstrap_B, seed,
          save_bootstrap_counts
        );
        // --- 終了時間計測 ---
        auto t_end = clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;
        double sec = elapsed.count();
        std::cerr << "[info] total runtime: " << sec << " seconds\n";

        return 0;
    }else if (edge_importance_mode) {
        try {
            Dataset ds_new = Dataset::fromCSV(score_dataset_path);
            DAG g_base = loadInitEdges(ds_new.D, init_path);
            std::vector<Counts> C_loaded = loadAllCountsTSV(counts_in_path, ds_new.D);

            computeEdgeImportanceScores(ds_new, g_base, C_loaded, alpha_ij, ess, save_edge_importance_path);

            // --- 終了時間計測 ---
            auto t_end = clock::now();
            std::chrono::duration<double> elapsed = t_end - t_start;
            double sec = elapsed.count();
            std::cerr << "[info] total runtime: " << sec << " seconds\n";
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "Error (edge-importance): " << e.what() << "\n";
            return 2;
        }
    }
    // 予測モード
    if (!score_dataset_path.empty()) {
        try {
            // 1) 新データセットのロード（CSV/TSV自動判定・ヘッダあり）
            Dataset ds_new = Dataset::fromCSV(score_dataset_path);

            // 2) 構造（DAG）のロード（TSV/スペース; 既存の loadInitEdges を利用）
            if (init_path.empty())
                throw runtime_error("--init <init_edges.tsv> is required for scoring mode");
            DAG g_new = loadInitEdges(ds_new.D, init_path);

            // 3) all_counts.tsv のロード（Counts 配列を復元）
            if (counts_in_path.empty())
                throw runtime_error("--counts <all_counts.tsv> is required for scoring mode");
            vector<int> child_r_from_counts;
            vector<Counts> C_loaded = loadAllCountsTSV(counts_in_path, ds_new.D, &child_r_from_counts);

            // 子の基数の簡易チェック（新データと counts の r_i が食い違っていないか）
            for (int v=0; v<ds_new.D; ++v){
                if (child_r_from_counts[v] > 0 && child_r_from_counts[v] != ds_new.r[v]) {
                    cerr << "[warn] child cardinality differs at node " << v
                         << " (counts r_i=" << child_r_from_counts[v]
                         << ", new-data r_i=" << ds_new.r[v] << ")\n";
                }
            }

            // 4) 尤度計算
            double avgN=0.0, avgND=0.0; long long zero_hits=0;
            double LL = computeLogLikelihoodOnDataset(ds_new, g_new, C_loaded, alpha_ij, &avgN, &avgND, &zero_hits);

            cout << fixed << setprecision(6);
            cout << "# scoring_mode\n";
            cout << "log_likelihood_total\t" << LL << "\n";
            cout << "log_likelihood_per_sample\t" << avgN << "\n";
            cout << "log_likelihood_per_variable\t" << avgND << "\n";
            cout << "alpha_ij\t" << alpha_ij << "\n";
            // --- 終了時間計測 ---
            auto t_end = clock::now();
            std::chrono::duration<double> elapsed = t_end - t_start;
            double sec = elapsed.count();
            std::cerr << "[info] total runtime: " << sec << " seconds\n";
            return 0;
        } catch (const exception& e){
            cerr << "Error (scoring): " << e.what() << "\n";
            return 2;
        }
    }
    // 単純な探索モード
    try{
        Dataset ds = Dataset::fromCSV(input_path);
        DAG init = loadInitEdges(ds.D, init_path);


        // --- 候補親Kの前処理（必要なときのみ） ---
        vector<vector<int>> topKlist;
        auto t_start_klist = clock::now();
        if (K>0 || mi_threshold>0.0 || chi2_p_threshold<1.0){
            mt19937_64 rng(seed);
	    // 行サンプルを作る（mi_sample==0 なら全行）
            vector<int> rows;
            if (mi_sample>0 && mi_sample < ds.N){
                rows.resize(ds.N); iota(rows.begin(), rows.end(), 0);
                shuffle(rows.begin(), rows.end(), rng);
                rows.resize(mi_sample); sort(rows.begin(), rows.end());
            } else {
                rows.resize(ds.N); iota(rows.begin(), rows.end(), 0);
            }
            int budget = (mi_budget>0? mi_budget : (ds.D-1));
            //topKlist = MICandidates::compute(ds, K, budget, rows, rng, mi_threshold);
	    topKlist = AssocCandidates::compute(
                ds, K, budget, rows, rng, cand_metric, mi_threshold, chi2_p_threshold
            );
            if (verbose) cerr << "[mi] topK="<<K<<" sample="<<rows.size()<<" budget="<<budget<<"\n";
        }
        {
	    for(int i=0;i<topKlist.size();++i){
		std::cout << i << " ";
	    	for(int j=0;j<topKlist[i].size();++j){
		    std::cout << topKlist[i][j] << " ";
		}
	    	std::cout<<std::endl;
	    }
	}
        // --- 終了時間計測 ---
        auto t_end_klist = clock::now();
        std::chrono::duration<double> elapsed_klist = t_end_klist - t_start_klist;
        double sec_klist = elapsed_klist.count();
        std::cerr << "[info] top-K list runtime: " << sec_klist << " seconds\n";

        // 到達性モード確定
        Reachability::Mode rmode = (reach_mode_str=="dense"? Reachability::DENSE : Reachability::LAZY);

        // 探索器のセットアップ
        auto t_start_iter = clock::now();
        HillClimber hc(ds, sc, ess, init, rmode, jindex_cache_cap);
        hc.max_iter = iters;
        hc.verbose = verbose;
        hc.tabu_tenure = tabu_tenure;
        hc.max_parents = max_parents;
        hc.max_children = max_children;
        if (K>0){ hc.candParents = move(topKlist); hc.topK = K; }

        auto [g, score, it] = hc.run(/*use_tabu=*/(tabu_tenure>0));

        // --- イテレーションあたり時間計測 ---
        auto t_end_iter = clock::now();
        std::chrono::duration<double> elapsed_iter = t_end_iter - t_start_iter;
        double sec_iter = elapsed_iter.count();
        std::cerr << "[info] time / iter: " << sec_iter/it << " seconds\n";
        
	// --- トータル終了時間計測 ---
        auto t_end = clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;
        double sec = elapsed.count();
        std::cerr << "[info] total time: " << sec << " seconds\n";


        // 結果出力
        cout << fixed << setprecision(6);
        cout << "# learned_score="<<score<<"\n";
        // グラフ出力
        /*
        cout << "# edges (u v as u->v):\n";
        for (auto &e : g.edges()) {
            if (!ds.var_names.empty()) {
                cout << ds.var_names[e.first] << "\t" << ds.var_names[e.second] << "\n";
            } else {
                cout << e.first << "\t" << e.second << "\n";
            }
        }
        */
        
        // ===== 次回初期化用にファイル保存（インデックス版） =====
        if (!save_init_path.empty()) {
            ofstream fout(save_init_path);
            if (!fout) {
                cerr << "Error: failed to open --save path: " << save_init_path << "\n";
                return 2;
            }
            for (auto &e : g.edges()) {
                fout << e.first << "\t" << e.second << "\n";  // 0始まり u v 形式
            }
            cerr << "[info] wrote init edges (indices) to " << save_init_path << "\n";
        }

        // =====  参考用に変数名でも保存（再利用不可・人間可読） =====
        if (!save_init_names_path.empty()) {
            ofstream foutn(save_init_names_path);
            if (!foutn) {
                cerr << "Error: failed to open --save-names path: " << save_init_names_path << "\n";
                return 2;
            }
            for (auto &e : g.edges()) {
                if (!ds.var_names.empty()) {
                    foutn << ds.var_names[e.first] << "\t" << ds.var_names[e.second] << "\n";
                } else {
                    // ヘッダが無い場合はインデックスでフォールバック
                    foutn << e.first << "\t" << e.second << "\n";
                }
            }
            cerr << "[info] wrote init edges (names) to " << save_init_names_path << "\n";
        }
        // 全ノード分のカウントを 1 ファイル（TSV）にまとめて出力
        if (!save_counts_path.empty()) {
            try {
                //saveAllCountsTSV(save_counts_path, ds, g, hc.nodeCounts);
                saveAllCountsTSV(save_counts_path, ds, g);
                std::cerr << "[info] wrote all CPT counts (TSV) to " << save_counts_path << "\n";
            } catch (const std::exception& e) {
                std::cerr << "Error: cannot write --save-counts: " << e.what() << "\n";
                return 2;
            }
        }

    } catch (const exception& e){
        cerr << "Error: " << e.what() << "\n";
        return 2;
    }
    return 0;
}


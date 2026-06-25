#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <iostream>

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
    std::vector<int> X_flat;
    
    int N=0, D=0;
    std::vector<int> r;             // 各列の基数
    std::vector<std::string> var_names;  // 変数名（CSV/TSV のヘッダ）
    // accessor
    inline int x(int n, int d) const noexcept {
        return X_flat[(size_t)n * D + d];
    }
    ~Dataset(){
        acc_delete();
    }
    void acc_copyin(void){
        const int* ds_x_ptr = X_flat.data();
        const int* ds_r_ptr = r.data();
        #pragma acc enter data copyin(ds_x_ptr[0:N*D],ds_r_ptr[0:D])
    }
    void acc_delete(void){
        const int* ds_x_ptr = X_flat.data();
        const int* ds_r_ptr = r.data();
        #pragma acc exit data delete(ds_x_ptr[0:N*D],ds_r_ptr[0:D])
    }
    // CSV または TSV を読み込み（区切り文字を自動判定）
    static Dataset fromCSV(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) throw std::runtime_error("Failed to open file: " + path);

        std::string line;
        std::vector<std::string> headers;
        std::vector<std::vector<int>> rows;
        bool first = true;
        char delim = ','; // デフォルトはカンマ

        while (getline(fin, line)) {
            if (line.empty()) continue;

            // 区切り文字の自動判定：最初の行の '\t' の有無で決める
            if (first) {
                if (line.find('\t') != std::string::npos)
                    delim = '\t';
            }

            std::stringstream ss(line);
            std::string tmp;

            if (first) {
                // --- 1行目はヘッダ ---
                while (getline(ss, tmp, delim)) {
                    // 前後の空白除去
                    tmp.erase(tmp.begin(), std::find_if(tmp.begin(), tmp.end(), [](unsigned char c){return !isspace(c);} ));
                    tmp.erase(std::find_if(tmp.rbegin(), tmp.rend(), [](unsigned char c){return !isspace(c);} ).base(), tmp.end());
                    headers.push_back(tmp);
                }
                first = false;
                continue;
            }

            // --- データ行 ---
            std::vector<int> row;
            while (getline(ss, tmp, delim)) {
                tmp.erase(tmp.begin(), std::find_if(tmp.begin(), tmp.end(), [](unsigned char c){return !isspace(c);} ));
                tmp.erase(std::find_if(tmp.rbegin(), tmp.rend(), [](unsigned char c){return !isspace(c);} ).base(), tmp.end());
                if (tmp.empty()) throw std::runtime_error("Empty cell in data row.");
                try {
                    row.push_back(stoi(tmp));
                } catch (...) {
                    throw std::runtime_error("Non-integer cell detected: \"" + tmp + "\"");
                }
            }
            if (!row.empty())
                rows.push_back(move(row));
        }

        if (rows.empty()) throw std::runtime_error("Empty dataset (no data rows).");
        int D = (int)rows[0].size();
        for (auto& r: rows)
            if ((int)r.size()!=D)
                throw std::runtime_error("Inconsistent column count in data rows.");

        if (!headers.empty() && (int)headers.size()!=D)
            throw std::runtime_error("Header count does not match column count.");

        // 各列の基数（最大値+1）
        std::vector<int> rcard(D,0);
        for (int j=0;j<D;++j){
            int mx=0;
            for (auto& row: rows) mx=std::max(mx, row[j]);
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

        std::cerr << "[info] Loaded " << ds.N << " samples, "
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



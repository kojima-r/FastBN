#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_dummy_expr.py
==================
バルク RNA-seq 風の**ダミー発現量データ**を生成する。パイプラインの動作確認・
使用例 (example_bulk/) 用で、**正解のネットワーク構造が既知**なので学習結果の
評価にも使える。

生成モデル:
  1. 遺伝子上の疎な DAG (真の構造) をランダムに生成する
     (各遺伝子は自分より前の遺伝子から最大 --max-parents 個の親を選ぶ)
  2. 各サンプルについて、DAG の順に潜在的な log2 発現量を生成する
         x_j = base_j + Σ_i w_ij (x_i - base_i) + group_effect_j(g) + ε_j
     一部の遺伝子 (--frac-responsive) は群 (条件) ごとのシフトを持ち、
     その影響は DAG の下流にも伝播する
  3. さらに DAG に参加しない無相関ノイズ遺伝子を --n-noise 個加える
     (前処理の分散フィルタの動作確認用)
  4. 2^x を発現量の期待値、ライブラリサイズ変動を掛けて
     負の二項分布からリードカウントを生成する

出力 (--outdir 以下):
  counts.tsv        : 行=遺伝子, 列=サンプル の生カウント行列
                      (先頭列 gene_id, gene_name, gene_length)
  sample_meta.tsv   : sample_id / group / replicate / library_size
  true_edges.tsv    : 真の DAG (parent / child / weight) ※遺伝子名表記
  target_genes.txt  : 注目遺伝子リスト (ハブ遺伝子と群応答遺伝子)

使用例:
  python3 make_dummy_expr.py --outdir data --n-genes 80 --n-noise 40 \\
      --groups Control,TreatA,TreatB,Combo --n-replicates 8 --seed 12345
"""

import argparse
import os
import sys

import numpy as np


def log(*args):
    print("[dummy]", *args, file=sys.stderr)


def build_dag(n_genes, max_parents, edge_prob, rng, w_min=0.5, w_max=1.1):
    """遺伝子 0..n-1 の順序に沿った疎な DAG を作り、親リストと重みを返す。

    返り値: parents[j] = [(i, w), ...]  (i < j)
    """
    parents = [[] for _ in range(n_genes)]
    for j in range(1, n_genes):
        # 候補は自分より前の遺伝子。近い遺伝子ほど選ばれやすくして
        # モジュール的な (局所的に密な) 構造にする。
        cand = np.arange(max(0, j - 25), j)
        if cand.size == 0:
            continue
        n_max = min(max_parents, cand.size)
        # 期待親数 = edge_prob * cand.size を n_max で打ち切る
        k = min(n_max, rng.binomial(cand.size, edge_prob))
        if k == 0:
            continue
        chosen = rng.choice(cand, size=k, replace=False)
        for i in chosen:
            w = rng.uniform(w_min, w_max) * rng.choice([-1.0, 1.0])
            parents[j].append((int(i), float(w)))
    n_edges = sum(len(p) for p in parents)
    log(f"真の DAG: ノード {n_genes}, エッジ {n_edges} "
        f"(平均入次数 {n_edges / n_genes:.2f})")
    return parents


def simulate_latent(parents, base, gene_sd, group_effect, group_of_sample,
                    signal_frac, rng):
    """DAG の順に潜在 log2 発現量 x[gene, sample] を生成する。

    各遺伝子の変動 (ベースラインからのずれ) を
        gene_sd * ( sqrt(h) * 親からの寄与(標準化) + sqrt(1-h) * 独立ノイズ )
    と構成する (h = signal_frac)。こうすると
      * どの遺伝子も周辺分散が gene_sd^2 に揃う (カスケードで分散が爆発しない)
      * 親子の相関が sqrt(h) 程度に制御される (h=0.6 なら約 0.77)
    となり、依存の強さを 1 つのパラメータで調整できる。
    """
    n_genes = len(parents)
    n_samples = len(group_of_sample)
    dev = np.zeros((n_genes, n_samples))     # 標準化された変動 (sd≈1)
    h = float(np.clip(signal_frac, 0.0, 0.95))
    for j in range(n_genes):
        eps = rng.normal(0.0, 1.0, size=n_samples)
        if parents[j]:
            c = np.zeros(n_samples)
            for i, w in parents[j]:
                c = c + w * dev[i]
            sd = c.std()
            c = c / sd if sd > 1e-12 else c
            dev[j] = np.sqrt(h) * c + np.sqrt(1.0 - h) * eps
        else:
            dev[j] = eps
    x = base[:, None] + gene_sd * dev + group_effect[:, group_of_sample]
    return x


def main():
    ap = argparse.ArgumentParser(
        description="バルク RNA-seq 風ダミーデータ (真の構造つき) を生成",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--outdir", default="data", help="出力ディレクトリ")
    ap.add_argument("--out-counts", default="counts.tsv", help="カウント行列のファイル名")
    ap.add_argument("--out-meta", default="sample_meta.tsv", help="サンプル情報のファイル名")
    ap.add_argument("--out-edges", default="true_edges.tsv", help="真の DAG のファイル名")
    ap.add_argument("--out-targets", default="target_genes.txt",
                    help="注目遺伝子リストのファイル名 (--outdir の親に置く場合は"
                         "相対パスを指定)")
    ap.add_argument("--n-genes", type=int, default=80,
                    help="DAG に参加する遺伝子数")
    ap.add_argument("--n-noise", type=int, default=40,
                    help="DAG に参加しないノイズ遺伝子数 (フィルタ動作確認用)")
    ap.add_argument("--groups", default="Control,TreatA,TreatB,Combo",
                    help="群 (条件) 名のカンマ区切り")
    ap.add_argument("--n-replicates", type=int, default=8, help="各群のサンプル数")
    ap.add_argument("--max-parents", type=int, default=3, help="真の DAG の最大親数")
    ap.add_argument("--edge-prob", type=float, default=0.10,
                    help="親候補が親になる確率 (構造の密度)")
    ap.add_argument("--frac-responsive", type=float, default=0.25,
                    help="群ごとの発現シフトを持つ遺伝子の割合")
    ap.add_argument("--effect-size", type=float, default=1.2,
                    help="群効果の大きさ (log2 スケールの標準偏差)")
    ap.add_argument("--base-mean", type=float, default=5.0,
                    help="log2 発現量のベースライン平均 (遺伝子ごとの発現水準)")
    ap.add_argument("--base-sd", type=float, default=2.0,
                    help="ベースライン発現水準の遺伝子間ばらつき")
    ap.add_argument("--signal-frac", type=float, default=0.65,
                    help="遺伝子の変動のうち親で説明される割合 h "
                         "(親子の相関は概ね sqrt(h))")
    ap.add_argument("--gene-sd", type=float, default=1.2,
                    help="構造に参加する遺伝子の log2 発現量の標準偏差")
    ap.add_argument("--noise-gene-sd", type=float, default=0.5,
                    help="ノイズ遺伝子の log2 発現量の標準偏差 "
                         "(構造遺伝子より小さくして分散フィルタで落ちるようにする)")
    ap.add_argument("--weight-min", type=float, default=0.5,
                    help="真の DAG の重み |w| の下限")
    ap.add_argument("--weight-max", type=float, default=1.1,
                    help="真の DAG の重み |w| の上限")
    ap.add_argument("--dispersion", type=float, default=0.15,
                    help="負の二項分布の過分散パラメータ (var = mu + d*mu^2)")
    ap.add_argument("--lib-sd", type=float, default=0.2,
                    help="ライブラリサイズ (総リード数) の log 正規ばらつき")
    ap.add_argument("--total-reads", type=float, default=2.0e7,
                    help="1 サンプルあたりの総リード数の基準値")
    ap.add_argument("--n-targets", type=int, default=6,
                    help="注目遺伝子リストに載せる遺伝子数")
    ap.add_argument("--seed", type=int, default=12345, help="乱数シード")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    groups = [g.strip() for g in args.groups.split(",") if g.strip()]
    if not groups:
        sys.exit("[dummy] エラー: --groups が空です")
    n_g, n_n = args.n_genes, args.n_noise
    n_total = n_g + n_n
    n_samples = len(groups) * args.n_replicates
    log(f"遺伝子 {n_total} (構造 {n_g} + ノイズ {n_n}), "
        f"サンプル {n_samples} = {len(groups)} 群 x {args.n_replicates} 反復")

    # --- 真の DAG ----------------------------------------------------------
    parents = build_dag(n_g, args.max_parents, args.edge_prob, rng,
                        args.weight_min, args.weight_max)

    # --- ベースライン発現量と群効果 ----------------------------------------
    base = rng.normal(args.base_mean, args.base_sd, size=n_g)
    base = np.clip(base, 0.5, None)
    n_resp = int(round(args.frac_responsive * n_g))
    responsive = rng.choice(n_g, size=n_resp, replace=False)
    group_effect = np.zeros((n_g, len(groups)))
    for j in responsive:
        # 対照群 (先頭) は 0 とし、他群にシフトを与える
        group_effect[j, 1:] = rng.normal(0.0, args.effect_size, size=len(groups) - 1)
    log(f"群応答遺伝子: {n_resp} 件 (効果量 sd={args.effect_size})")

    group_of_sample = np.repeat(np.arange(len(groups)), args.n_replicates)

    # --- 潜在 log2 発現量 -> カウント --------------------------------------
    x = simulate_latent(parents, base, args.gene_sd, group_effect,
                        group_of_sample, args.signal_frac, rng)
    log(f"親子の相関 ≈ {np.sqrt(min(args.signal_frac, 0.95)):.2f} "
        f"(signal-frac={args.signal_frac})")
    # ノイズ遺伝子: 構造も群効果も持たず、変動も小さい (分散フィルタで落ちる)
    if n_n > 0:
        nbase = rng.normal(args.base_mean - 1.0, args.base_sd, size=n_n)
        nbase = np.clip(nbase, 0.2, None)
        xn = nbase[:, None] + rng.normal(0.0, args.noise_gene_sd,
                                         size=(n_n, n_samples))
        x = np.vstack([x, xn])

    mu_rel = np.power(2.0, np.clip(x, -5.0, 16.0))          # 相対発現量
    lib = rng.lognormal(mean=0.0, sigma=args.lib_sd, size=n_samples)
    # 各サンプル内の相対量を「総リード数 x ライブラリサイズ変動」にスケールする
    # (実際のシーケンスと同様、サンプルごとの総リード数がばらつく)
    col_sum = np.maximum(mu_rel.sum(axis=0), 1e-9)
    depth = args.total_reads * lib / col_sum
    mu = mu_rel * depth[None, :]
    r = 1.0 / max(args.dispersion, 1e-9)                     # NB の形状パラメータ
    counts = rng.negative_binomial(n=r, p=r / (r + mu)).astype(int)

    # --- 名前・ファイル出力 ------------------------------------------------
    names = [f"G{j + 1:03d}" for j in range(n_g)] + [f"NOISE{j + 1:03d}" for j in range(n_n)]
    gene_ids = [f"DUMMYG{j + 1:05d}" for j in range(n_total)]
    gene_len = rng.integers(500, 8000, size=n_total)
    sample_ids = [f"S{i + 1:02d}_{groups[group_of_sample[i]]}" for i in range(n_samples)]

    os.makedirs(args.outdir, exist_ok=True)
    counts_path = os.path.join(args.outdir, args.out_counts)
    with open(counts_path, "w", encoding="utf-8") as fp:
        fp.write("gene_id\tgene_name\tgene_length\t" + "\t".join(sample_ids) + "\n")
        for j in range(n_total):
            fp.write(f"{gene_ids[j]}\t{names[j]}\t{gene_len[j]}\t"
                     + "\t".join(str(v) for v in counts[j]) + "\n")
    log(f"出力: {counts_path} ({n_total} 遺伝子 x {n_samples} サンプル, 生カウント)")

    meta_path = os.path.join(args.outdir, args.out_meta)
    with open(meta_path, "w", encoding="utf-8") as fp:
        fp.write("sample_id\tgroup\treplicate\tlibrary_size\n")
        for i in range(n_samples):
            rep = i % args.n_replicates + 1
            fp.write(f"{sample_ids[i]}\t{groups[group_of_sample[i]]}\t{rep}\t"
                     f"{int(counts[:, i].sum())}\n")
    log(f"出力: {meta_path} (サンプル情報)")

    edges_path = os.path.join(args.outdir, args.out_edges)
    n_edges = 0
    with open(edges_path, "w", encoding="utf-8") as fp:
        fp.write("parent\tchild\tweight\n")
        for j, plist in enumerate(parents):
            for i, w in plist:
                fp.write(f"{names[i]}\t{names[j]}\t{w:.4f}\n")
                n_edges += 1
    log(f"出力: {edges_path} (真の DAG, {n_edges} エッジ)")

    # 注目遺伝子: 次数の高いハブ + 群応答遺伝子を混ぜる
    degree = np.zeros(n_g, dtype=int)
    for j, plist in enumerate(parents):
        degree[j] += len(plist)
        for i, _ in plist:
            degree[i] += 1
    hub_order = np.argsort(degree)[::-1]
    n_hub = max(1, args.n_targets // 2)
    picks = list(hub_order[:n_hub])
    resp_sorted = [j for j in hub_order if j in set(responsive.tolist())]
    for j in resp_sorted:
        if len(picks) >= args.n_targets:
            break
        if j not in picks:
            picks.append(j)
    # 注目遺伝子リストは解析ディレクトリ直下に置くことが多いため、パス区切りを
    # 含む指定 (例 ./target_genes.txt) はカレント基準、単なるファイル名は
    # --outdir 基準として扱う。
    tgt_path = (args.out_targets
                if (os.path.isabs(args.out_targets) or os.sep in args.out_targets)
                else os.path.join(args.outdir, args.out_targets))
    os.makedirs(os.path.dirname(os.path.abspath(tgt_path)) or ".", exist_ok=True)
    with open(tgt_path, "w", encoding="utf-8") as fp:
        fp.write("# ダミーデータの注目遺伝子 (ハブ + 群応答遺伝子)\n")
        fp.write("# 1 行 1 遺伝子。分散フィルタを免除して必ず解析対象に残す。\n")
        for j in picks:
            fp.write(f"{names[j]}\n")
    log(f"出力: {tgt_path} (注目遺伝子 {len(picks)} 件: "
        f"{', '.join(names[j] for j in picks)})")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute CPT from all_counts.tsv produced by --save-counts.

Input format (TSV):
  - Comment/meta lines start with '#'
  - Data lines: v \t j \t k|* \t n
      k='*' means n_ij (sum over child states)

Output:
  - Single TSV (--out) with rows: v j k prob
  - Or per-node TSV files (--out-dir), named cpt_<v>.tsv
Smoothing:
  P = (n_ijk + alpha/r_i) / (n_ij + alpha), alpha>=0
"""

import argparse
import math
import os
import sys
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute CPTs from all_counts.tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--counts", required=True, help="Input all_counts.tsv path")
    out_grp = p.add_mutually_exclusive_group(required=False)
    out_grp.add_argument("--out", help="Output single TSV (v j k prob). If omitted and --out-dir not set, prints to stdout")
    out_grp.add_argument("--out-dir", help="Output directory for per-node TSV files (cpt_<v>.tsv)")
    p.add_argument("--alpha", type=float, default=0.0, help="Dirichlet smoothing alpha (0.0 = MLE)")
    p.add_argument("--precision", type=int, default=12, help="Floating output precision")
    p.add_argument("--skip-nan", action="store_true",
                   help="Skip rows where nij=0 and alpha=0 (would yield NaN); default prints NaN")
    p.add_argument("--quiet", action="store_true", help="Suppress progress messages")
    return p.parse_args()

def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)

def main():
    args = parse_args()
    counts_path = args.counts
    alpha = float(args.alpha)
    prec = int(args.precision)
    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)

    # Sparse accumulators
    # n_ijk[(v,j,k)] = count ; n_ij[(v,j)] = count
    n_ijk = defaultdict(int)
    n_ij  = defaultdict(int)

    # Track per node maximum k and j to infer r_i, q_i
    max_k_for_v = defaultdict(int)    # r_i = max_k+1
    max_j_for_v = defaultdict(int)    # q_i = max_j+1

    total_lines = 0
    data_lines = 0

    # Stream read
    with open(counts_path, "r", encoding="utf-8") as f:
        for line in f:
            total_lines += 1
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                # tolerant: try splitting by whitespace
                parts = line.split()
                if len(parts) != 4:
                    continue
            vstr, jstr, kstr, nstr = parts
            try:
                v = int(vstr); j = int(jstr); n = int(nstr)
            except Exception:
                # non-data line
                continue

            if kstr == "*" or kstr == "'*'":
                n_ij[(v, j)] += n
                if j > max_j_for_v[v]: max_j_for_v[v] = j
                data_lines += 1
                continue

            try:
                k = int(kstr)
            except Exception:
                # malformed
                continue

            n_ijk[(v, j, k)] += n
            if j > max_j_for_v[v]: max_j_for_v[v] = j
            if k > max_k_for_v[v]: max_k_for_v[v] = k
            data_lines += 1

    if not args.quiet:
        eprint(f"[info] read {data_lines} data lines (total {total_lines}) from {counts_path}")

    # Helper to compute CPT rows for a node v
    def cpt_rows_for_node(v):
        # infer r_i, q_i from maxima (fallback to 1 if absent)
        r_i = (max_k_for_v[v] + 1) if v in max_k_for_v else 1
        q_i = (max_j_for_v[v] + 1) if v in max_j_for_v else 1

        # For each j in [0, q_i):
        for j in range(q_i):
            # nij: prefer explicit n_ij; otherwise sum over k present
            nij = n_ij.get((v, j), None)
            if nij is None:
                # compute from n_ijk sparse entries
                s = 0
                # iterate only known ks; to avoid O(r_i) scan over all k, collect from dict
                # but dict is sparse; generate keys present for this (v,j,*)
                # This is still O(#present k for (v,j,*))
                for kk in range(r_i):
                    s += n_ijk.get((v, j, kk), 0)
                nij = s

            denom = nij + alpha
            for k in range(r_i):
                nijk = n_ijk.get((v, j, k), 0)
                numer = nijk + alpha / r_i if alpha > 0.0 else nijk

                if denom == 0.0:
                    # alpha==0 and nij==0 → undefined (NaN) or skip
                    if alpha == 0.0:
                        if args.skip_nan:
                            continue
                        prob = float('nan')
                    else:
                        prob = 1.0 / r_i  # (0 + a/r) / (0 + a)
                else:
                    prob = numer / denom

                yield (v, j, k, prob)

    # Output
    fmt = f"{{:.{prec}g}}"  # compact with given precision

    if args.out_dir:
        # per-node files
        written = 0
        # nodes present in either dict:
        vs_present = set(v for (v, _) in n_ij.keys()) | set(v for (v, _, _) in n_ijk.keys())
        if not vs_present:
            # fallback to 0..max_v
            vs_present = set(list(max_j_for_v.keys()) + list(max_k_for_v.keys()))
        for v in sorted(vs_present):
            path = os.path.join(args.out_dir, f"cpt_{v}.tsv")
            with open(path, "w", encoding="utf-8") as g:
                g.write("v\tj\tk\tprob\n")
                for (vv, j, k, p) in cpt_rows_for_node(v):
                    g.write(f"{vv}\t{j}\t{k}\t{fmt.format(p)}\n")
            written += 1
        if not args.quiet:
            eprint(f"[info] wrote CPTs for {written} nodes into: {args.out_dir}")
    else:
        # single TSV to file or stdout
        out_stream = open(args.out, "w", encoding="utf-8") if args.out else sys.stdout
        close_needed = (out_stream is not sys.stdout)
        try:
            out_stream.write("v\tj\tk\tprob\n")
            # nodes present (as above)
            vs_present = set(v for (v, _) in n_ij.keys()) | set(v for (v, _, _) in n_ijk.keys())
            if not vs_present:
                vs_present = set(list(max_j_for_v.keys()) + list(max_k_for_v.keys()))
            for v in sorted(vs_present):
                for (vv, j, k, p) in cpt_rows_for_node(v):
                    out_stream.write(f"{vv}\t{j}\t{k}\t{fmt.format(p)}\n")
        finally:
            if close_needed:
                out_stream.close()
        if not args.quiet:
            target = args.out if args.out else "<stdout>"
            eprint(f"[info] wrote CPTs to {target}")

if __name__ == "__main__":
    main()


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute CPTs (with node/parent names and parent variable values)
from all_counts.tsv produced by --save-counts.

Output columns:
  node_id  node_name  parent_ids  parent_names  parent_values  j  k  prob
"""

import argparse
import os
import sys
import math
from collections import defaultdict
from itertools import product

def parse_args():
    p = argparse.ArgumentParser(description="Compute CPT with parent values from all_counts.tsv")
    p.add_argument("--counts", required=True, help="Input all_counts.tsv")
    p.add_argument("--out", help="Output single TSV file (default: stdout)")
    p.add_argument("--out-dir", help="Output directory for per-node TSVs (cpt_<v>.tsv)")
    p.add_argument("--alpha", type=float, default=0.0, help="Dirichlet smoothing (default=0.0)")
    p.add_argument("--skip-nan", action="store_true", help="Skip rows where nij=0 and alpha=0")
    p.add_argument("--precision", type=int, default=12, help="Decimal precision")
    p.add_argument("--quiet", action="store_true", help="Suppress info logs")
    return p.parse_args()

def eprint(*a, **k): print(*a, file=sys.stderr, **k)

def decode_parent_values(j, parent_card):
    """Return tuple of parent values corresponding to index j."""
    if not parent_card:
        return []
    vals = []
    for r in reversed(parent_card):
        vals.append(j % r)
        j //= r
    return list(reversed(vals))

def main():
    args = parse_args()
    alpha = args.alpha
    prec = args.precision

    n_ijk = defaultdict(int)
    n_ij = defaultdict(int)
    max_k_for_v = defaultdict(int)
    max_j_for_v = defaultdict(int)
    node_name = {}
    parent_ids = defaultdict(list)
    parent_names = defaultdict(list)
    parent_card = defaultdict(list)   # v -> [ri of each parent]

    current_node = None

    with open(args.counts, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                s = line[1:].strip()
                if s.startswith("--- node"):
                    try:
                        current_node = int(s.split()[2])
                    except Exception:
                        current_node = None
                elif s.startswith("node_name") and current_node is not None:
                    _, val = s.split("\t", 1)
                    node_name[current_node] = val.strip()
                elif s.startswith("parents_indices") and current_node is not None:
                    try:
                        _, val = s.split("\t", 1)
                        ids = [int(x.strip()) for x in val.split(",") if x.strip()]
                        parent_ids[current_node] = ids
                    except Exception:
                        parent_ids[current_node] = []
                elif s.startswith("parents_names") and current_node is not None:
                    try:
                        _, val = s.split("\t", 1)
                        names = [x.strip() for x in val.split(",") if x.strip()]
                        parent_names[current_node] = names
                    except Exception:
                        parent_names[current_node] = []
                elif s.startswith("parents_cardinalities") and current_node is not None:
                    # optional metadata if available
                    try:
                        _, val = s.split("\t", 1)
                        parent_card[current_node] = [int(x.strip()) for x in val.split(",") if x.strip()]
                    except Exception:
                        pass
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            vstr, jstr, kstr, nstr = parts
            try:
                v = int(vstr); j = int(jstr); n = int(nstr)
            except:
                continue
            if kstr == "*" or kstr == "'*'":
                n_ij[(v, j)] += n
                max_j_for_v[v] = max(max_j_for_v[v], j)
            else:
                try:
                    k = int(kstr)
                except:
                    continue
                n_ijk[(v, j, k)] += n
                max_j_for_v[v] = max(max_j_for_v[v], j)
                max_k_for_v[v] = max(max_k_for_v[v], k)

    vs = set(v for (v, _) in n_ij.keys()) | set(v for (v, _, _) in n_ijk.keys())
    fmt = f"{{:.{prec}g}}"

    def cpt_rows(v):
        r_i = (max_k_for_v[v] + 1) if v in max_k_for_v else 1
        q_i = (max_j_for_v[v] + 1) if v in max_j_for_v else 1
        pids = parent_ids.get(v, [])
        pnames = parent_names.get(v, [])
        pcards = parent_card.get(v, [])

        for j in range(q_i):
            nij = n_ij.get((v, j), sum(n_ijk.get((v, j, k), 0) for k in range(r_i)))
            denom = nij + alpha
            pvals = decode_parent_values(j, pcards) if pcards else []
            for k in range(r_i):
                nijk = n_ijk.get((v, j, k), 0)
                numer = nijk + alpha / r_i if alpha > 0 else nijk
                if denom == 0:
                    if alpha == 0:
                        if args.skip_nan:
                            continue
                        prob = float("nan")
                    else:
                        prob = 1.0 / r_i
                else:
                    prob = numer / denom
                yield (v, node_name.get(v, f"X{v}"),
                       ",".join(map(str, pids)),
                       ",".join(pnames),
                       ",".join(map(str, pvals)) if pvals else "",
                       j, k, prob)

    header = "node_id\tnode_name\tparent_ids\tparent_names\tparent_values\tj\tk\tprob\n"

    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)
        for v in sorted(vs):
            with open(os.path.join(args.out_dir, f"cpt_{v}.tsv"), "w", encoding="utf-8") as g:
                g.write(header)
                for row in cpt_rows(v):
                    g.write("\t".join(map(str, row[:-1])) + f"\t{fmt.format(row[-1])}\n")
        if not args.quiet:
            eprint(f"[info] wrote CPTs with parent values into {args.out_dir}")
    else:
        out = open(args.out, "w", encoding="utf-8") if args.out else sys.stdout
        close = (out is not sys.stdout)
        try:
            out.write(header)
            for v in sorted(vs):
                for row in cpt_rows(v):
                    out.write("\t".join(map(str, row[:-1])) + f"\t{fmt.format(row[-1])}\n")
        finally:
            if close:
                out.close()
        if not args.quiet:
            eprint(f"[info] CPTs written to {args.out or '<stdout>'}")

if __name__ == "__main__":
    main()


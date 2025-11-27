#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute CPTs (with node/parent names and parent variable values)
from all_counts.tsv produced by --save-counts.

python compute_bs_prob.py --input example/bs_bak/edges_seed*.tsv --threshold-prob 0.2 --remove-cycle --out out.tsv

"""

import argparse
import os
import sys

from collections import defaultdict

def has_cycle(edges):
    # 隣接リストを作成
    graph = defaultdict(list)
    nodes = set()
    for u, v in edges:
        graph[u].append(v)
        nodes.add(u)
        nodes.add(v)

    visited = set()
    in_stack = set()

    # DFS（Graph may be disconnected）
    def dfs(node):
        visited.add(node)
        in_stack.add(node)

        for nei in graph[node]:
            # 未訪問 → 探索
            if nei not in visited:
                if dfs(nei):
                    return True
            # 探索中のノードに再度到達 → サイクル
            elif nei in in_stack:
                return True

        in_stack.remove(node)
        return False

    # 全ノードから開始（非連結グラフ対応）
    for node in nodes:
        if node not in visited:
            if dfs(node):
                return True

    return False

def parse_args():
    p = argparse.ArgumentParser(description="Compute bootstrap probability from edge_seed****.tsv")
    p.add_argument("--input", required=True, nargs="*", type=str, help="directory .tsv")
    p.add_argument("--out", type=str, default="", help="Output single TSV file (default: stdout)")
    p.add_argument("--threshold-prob", type=float, default=0.2, help="bootstrap probability threshold (default=0.1)")
    p.add_argument("--threshold-count", type=int, default=2, help="bootstrap count threshold (default=2)")
    p.add_argument("--remove-cycle", action="store_true", help="Suppress info logs")
    p.add_argument("--sort-by-prob", action="store_true", help="Suppress info logs")
    p.add_argument("--quiet", action="store_true", help="Suppress info logs")
    return p.parse_args()

def load_graph(args):
    graph = {}
    m = len(args.input)
    for filename in args.input:
        if not args.quiet:
            print("[LOAD]",filename)
        edge_count=0
        fp = open(filename)
        head_s = next(fp)
        head = head_s.strip().split("\t")
        for line in fp:
            arr=line.strip().split("\t")
            if len(arr)>=4:
                u=int(arr[0])
                v=int(arr[1])
                cnt=int(arr[2])
                prob=float(arr[3])
                key=(u,v)
                if key not in graph:
                    graph[key]=[]
                graph[key].append((cnt, prob))
                edge_count+=1
        if not args.quiet:
            print("#edge:",edge_count)
    return graph

def summarize_graph(args, graph, m):
    threshold_prob = args.threshold_prob
    threshold_count = args.threshold_count
    new_graph={}
    for k, v in graph.items():
        v+=[(0,0)]*(m-len(v))
        all_cnt=sum([cnt for cnt, prob in v])
        all_prob=sum([prob for cnt, prob in v])/len(v)
        if all_prob >= threshold_prob and all_cnt >= threshold_count:
            new_graph[k]= (all_cnt, all_prob)
    return new_graph

def remove_cycle(args, sum_graph):
    pre_sorted_graph=[]
    for e,v in sum_graph.items():
        sort_key = v[1] # all_prob
        pre_sorted_graph.append((sort_key, e))
    sorted_graph=sorted(pre_sorted_graph, reverse=True)
    cycle_detect=False
    new_graph=[]
    for edge in sorted_graph:
        temp_g=new_graph+[edge]
        if has_cycle([e for key, e in temp_g]):
            cycle_detect=True
        else:
            new_graph.append(edge)

    if cycle_detect:
        if not args.quiet:
            print("removed cycle graph: #edge =", len(new_graph), "/ #original =", len(sorted_graph))
    return new_graph

def main():
    args = parse_args()
    graph = load_graph(args)
    m = len(args.input)
    sum_graph = summarize_graph(args, graph, m)
    if args.remove_cycle:
        g=remove_cycle(args, sum_graph)
    else:
        g=[(0,e) for e,v in sum_graph.items()]
        print("output graph: #edge =", len(sum_graph))

    out_sort_by_prob=args.sort_by_prob
    new_graph=[]
    for _,e in g:
        val=sum_graph[e]
        all_cnt = val[0]
        all_prob = val[1]
        obj =list(map(str, [e[0], e[1], all_cnt, all_prob]))

        if out_sort_by_prob:
            sort_key = all_prob
        else:
            sort_key = e
        new_graph.append((sort_key, obj))
    
    if args.out!="":
        with open(args.out,"w") as ofp:
            for _, v in sorted(new_graph, reverse=out_sort_by_prob):
                s="\t".join(v)
                ofp.write(s)
                ofp.write("\n")
    else:
        for _, v in sorted(new_graph, reverse=out_sort_by_prob):
            s="\t".join(v)
            print(s)

if __name__ == "__main__":
    main()


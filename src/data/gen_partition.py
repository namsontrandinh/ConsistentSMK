#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate random balanced partition of graph nodes from an edge-list file.

Input edge list format (text):
    u v
(one edge per line, whitespace separated; lines starting with # are ignored)

Usage:
    python gen_part.py <edge_list.txt> <p> <seed>

Output:
    <edge_list>-part.txt
    n lines (n = number of distinct nodes), each line:
        i j
    where i = node id (original), j = group id in [0, p-1]
"""

import os
import sys
import random
from typing import List, Set


def read_nodes(edge_path: str) -> List[int]:
    nodes: Set[int] = set()
    with open(edge_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                raise ValueError(f"Line {line_no}: expected two integers 'u v'")
            u = int(parts[0])
            v = int(parts[1])
            nodes.add(u)
            nodes.add(v)
    return list(nodes)


def balanced_random_partition(nodes: List[int], p: int, seed: int) -> List[tuple[int, int]]:
    if p <= 0:
        raise ValueError("p must be a positive integer")
    if not nodes:
        return []

    rng = random.Random(seed)
    rng.shuffle(nodes)

    n = len(nodes)
    base = n // p
    rem = n % p

    out: List[tuple[int, int]] = []
    idx = 0
    for g in range(p):
        size_g = base + (1 if g < rem else 0)
        for _ in range(size_g):
            out.append((nodes[idx], g))
            idx += 1
    return out


def make_out_path(edge_path: str) -> str:
    base, _ = os.path.splitext(edge_path)
    return f"{base}-part.txt"


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage: python gen_part.py <edge_list.txt> <p> <seed>")
        return 1

    edge_path = sys.argv[1]
    try:
        p = int(sys.argv[2])
        seed = int(sys.argv[3])
    except ValueError:
        print("Error: <p> and <seed> must be integers.")
        return 1

    try:
        nodes = read_nodes(edge_path)
        pairs = balanced_random_partition(nodes, p, seed)
    except Exception as e:
        print(f"Error: {e}")
        return 1

    out_path = make_out_path(edge_path)
    with open(out_path, "w", encoding="utf-8") as fo:
        for u, g in pairs:
            fo.write(f"{u} {g}\n")

    print(f"Input : {edge_path}")
    print(f"|V|   : {len(nodes)}")
    print(f"p     : {p}")
    print(f"seed  : {seed}")
    print(f"Output: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

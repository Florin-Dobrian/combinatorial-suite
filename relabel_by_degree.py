#!/usr/bin/env python3
"""Relabel graph vertices by degree.
Usage: python3 relabel_by_degree.py input.txt output.txt [--asc|--desc]

  --asc   Low degree first (satellites before hubs) [default]
  --desc  High degree first (hubs before satellites)
"""
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 relabel_by_degree.py input.txt output.txt [--asc|--desc]")
        sys.exit(1)

    descending = False
    for arg in sys.argv[3:]:
        if arg == "--desc":
            descending = True
        elif arg == "--asc":
            descending = False

    with open(sys.argv[1]) as f:
        n, m = map(int, f.readline().split())
        edges = []
        deg = [0] * n
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                u, v = int(parts[0]), int(parts[1])
                edges.append((u, v))
                deg[u] += 1
                deg[v] += 1

    # Sort vertices by degree
    if descending:
        order = sorted(range(n), key=lambda v: (-deg[v], v))
    else:
        order = sorted(range(n), key=lambda v: (deg[v], v))

    new_id = [0] * n
    for i, old in enumerate(order):
        new_id[old] = i

    with open(sys.argv[2], 'w') as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{new_id[u]} {new_id[v]}\n")

    mode = "descending (hubs first)" if descending else "ascending (satellites first)"
    print(f"Relabeled {n} vertices, {len(edges)} edges ({mode})")
    print(f"Degree range: {min(deg)} - {max(deg)}")

if __name__ == "__main__":
    main()

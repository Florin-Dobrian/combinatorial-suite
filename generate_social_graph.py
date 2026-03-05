#!/usr/bin/env python3
"""
Generate social-network-style graphs: a k-clique hub with l satellite vertices.

Structure:
  - k hub vertices forming a complete graph (k*(k-1)/2 edges)
  - l satellite vertices, each connected to all k hub vertices (l*k edges)
  - No edges among satellites
  - Total edges: k*(k-1)/2 + l*k
  - Maximum matching: k (each hub matched to one satellite)
  - Matched vertices: 2k out of k+l total

Vertex labeling:
  --hub-first: hubs are 0..k-1, satellites are k..k+l-1  (default)
  --hub-last:  satellites are 0..l-1, hubs are l..l+k-1

Output: V E header, then u v edges (general format for GO/MV).

Usage:
  python generate_social_graph.py -k 100 -l 10000 -o social_100_10000.txt
  python generate_social_graph.py -k 100 -l 10000 --hub-last -o social_100_10000_hl.txt
"""

import argparse
import sys


def generate_social_graph(k, l, hub_first=True):
    """Generate edges for hub-clique + satellite graph.
    
    Returns (V, edges) where edges is list of (u, v) with u < v.
    """
    V = k + l
    edges = []

    if hub_first:
        # Hubs: 0..k-1, Satellites: k..k+l-1
        hub_start = 0
        sat_start = k
    else:
        # Satellites: 0..l-1, Hubs: l..l+k-1
        sat_start = 0
        hub_start = l

    # Hub clique edges
    for i in range(k):
        for j in range(i + 1, k):
            u, v = hub_start + i, hub_start + j
            edges.append((min(u, v), max(u, v)))

    # Satellite-to-hub edges
    for s in range(l):
        for h in range(k):
            u, v = sat_start + s, hub_start + h
            edges.append((min(u, v), max(u, v)))

    edges.sort()
    return V, edges


def main():
    parser = argparse.ArgumentParser(
        description="Generate social-network-style hub-clique + satellite graphs"
    )
    parser.add_argument("-k", type=int, required=True,
                        help="Number of hub vertices (clique size)")
    parser.add_argument("-l", type=int, required=True,
                        help="Number of satellite vertices")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output filename")
    parser.add_argument("--hub-last", action="store_true", default=False,
                        help="Place hub vertices after satellites (default: hub first)")
    args = parser.parse_args()

    k, l = args.k, args.l
    if l <= k:
        print(f"Warning: l ({l}) <= k ({k}). Matching will be limited by satellites.", 
              file=sys.stderr)

    hub_first = not args.hub_last
    V, edges = generate_social_graph(k, l, hub_first=hub_first)
    E = len(edges)

    clique_edges = k * (k - 1) // 2
    satellite_edges = l * k
    match_pct = 100.0 * (2 * k) / V if V > 0 else 0
    avg_deg = 2.0 * E / V if V > 0 else 0

    with open(args.output, "w") as f:
        f.write(f"{V} {E}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")

    label = "hub-first" if hub_first else "hub-last"
    print(f"Generated: k={k} hubs, l={l} satellites, V={V}, E={E} "
          f"(clique={clique_edges}, sat={satellite_edges})")
    print(f"  Avg degree: {avg_deg:.1f}, Max matching: {k}, "
          f"Matched vertices: {2*k}/{V} ({match_pct:.2f}%)")
    print(f"  Labeling: {label}")
    print(f"  Written to: {args.output}")


if __name__ == "__main__":
    main()

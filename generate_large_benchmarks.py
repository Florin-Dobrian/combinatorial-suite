#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# ///
"""
Large-scale sparse graph generator for the combinatorial-suite benchmarks.

Generates both general and bipartite graphs at the sizes and sparsities
defined in large_scale_benchmarks.md. All graphs use the same file
formats as the existing data/ directory:

  General:   "<V> <E>\n" then one "u v\n" per edge, 0-indexed.
  Bipartite: "<L> <R> <E>\n" then one "l r\n" per edge, 0-indexed.

Usage:
    python generate_large_benchmarks.py                  # generate all
    python generate_large_benchmarks.py --sizes 100k 1m  # specific sizes
    python generate_large_benchmarks.py --type general    # general only
    python generate_large_benchmarks.py --list            # show plan, don't generate

Output directory: data/large-benchmarks/  (created if needed)
"""

import argparse
import os
import random
import sys
import time


# ── Benchmark specification ──────────────────────────────────────────

BENCHMARKS = [
    # (label,       V,         c)    where E = c * V
    ("10k_2",        10_000,   2),
    ("10k_3",        10_000,   3),
    ("10k_5",        10_000,   5),
    ("10k_10",       10_000,  10),
    ("10k_15",       10_000,  15),
    ("50k_2",        50_000,   2),
    ("50k_3",        50_000,   3),
    ("50k_5",        50_000,   5),
    ("50k_10",       50_000,  10),
    ("50k_15",       50_000,  15),
    ("100k_2",      100_000,   2),
    ("100k_3",      100_000,   3),
    ("100k_5",      100_000,   5),
    ("100k_10",     100_000,  10),
    ("100k_15",     100_000,  15),
    ("500k_2",      500_000,   2),
    ("500k_3",      500_000,   3),
    ("500k_5",      500_000,   5),
    ("500k_10",     500_000,  10),
    ("500k_15",     500_000,  15),
    ("1m_2",      1_000_000,   2),
    ("1m_3",      1_000_000,   3),
    ("1m_5",      1_000_000,   5),
    ("1m_10",     1_000_000,  10),
    ("1m_15",     1_000_000,  15),
    ("5m_2",      5_000_000,   2),
    ("5m_3",      5_000_000,   3),
    ("5m_5",      5_000_000,   5),
    ("5m_10",     5_000_000,  10),
    ("5m_15",     5_000_000,  15),
    ("10m_2",    10_000_000,   2),
    ("10m_3",    10_000_000,   3),
    ("10m_5",    10_000_000,   5),
    ("10m_10",   10_000_000,  10),
    ("10m_15",   10_000_000,  15),
]

SEED_BASE = 42


# ── Graph generators ─────────────────────────────────────────────────

def generate_general_sparse(V, E, seed):
    """
    Generate a simple undirected graph with V vertices and exactly E edges.
    Uses rejection sampling on random vertex pairs. For sparse graphs
    (E << V^2) collisions are rare, so this runs in O(E) expected time.
    Guarantees: no self-loops, no multi-edges.
    """
    rng = random.Random(seed)
    edges = set()
    while len(edges) < E:
        u = rng.randrange(V)
        v = rng.randrange(V)
        if u == v:
            continue
        edge = (min(u, v), max(u, v))
        edges.add(edge)
    return sorted(edges)


def generate_bipartite_sparse(L, R, E, seed):
    """
    Generate a simple bipartite graph with L left vertices, R right
    vertices, and exactly E edges. Left vertices are 0..L-1, right
    vertices are 0..R-1 (separate index spaces, as in the repo format).
    """
    rng = random.Random(seed)
    max_edges = L * R
    if E > max_edges:
        raise ValueError(f"E={E} exceeds maximum L*R={max_edges}")
    edges = set()
    while len(edges) < E:
        u = rng.randrange(L)
        v = rng.randrange(R)
        edges.add((u, v))
    return sorted(edges)


def generate_general_bipartite_sparse(V, E, edge_seed, shuffle_seed):
    """
    Generate a bipartite graph stored in general format with shuffled
    vertex labels. Produces the SAME edges as generate_bipartite_sparse
    (using the same edge_seed), then randomly permutes all vertex IDs
    (using shuffle_seed) so the bipartite structure is hidden.

    Left vertices: 0..V/2-1, Right vertices: V/2..V-1 (before shuffle).
    Output: sorted (u, v) edges with u < v, for V E header format.
    """
    L = V // 2
    R = V - L

    # Generate edges with same seed as bipartite format
    # bipartite format: left 0..L-1, right 0..R-1 (separate spaces)
    # here: left 0..L-1, right L..L+R-1 (single space)
    rng = random.Random(edge_seed)
    edges_bip = set()
    while len(edges_bip) < E:
        u = rng.randrange(L)
        v = rng.randrange(R)
        edges_bip.add((u, v))

    # Map to single vertex space: left u stays u, right v becomes L + v
    edges_single = [(u, L + v) for u, v in edges_bip]

    # Shuffle vertex labels
    perm = list(range(V))
    shuffle_rng = random.Random(shuffle_seed)
    shuffle_rng.shuffle(perm)

    shuffled = set()
    for u, v in edges_single:
        a, b = perm[u], perm[v]
        if a > b:
            a, b = b, a
        shuffled.add((a, b))

    return sorted(shuffled)


def generate_stacked_bipartite_sparse(V, E, edge_seed):
    """
    Generate a bipartite graph stored in general format with stacked
    vertex labels. Left vertices: 0..L-1, right vertices: L..L+R-1.
    Same edges as generate_bipartite_sparse (same edge_seed), just
    with right vertices offset by L. No shuffling.

    This preserves the vertex ordering so greedy-md produces the same
    matching quality as the bipartite format, enabling fair comparison
    between bipartite-aware (HK) and general (MV) algorithms.
    """
    L = V // 2
    R = V - L

    rng = random.Random(edge_seed)
    edges_bip = set()
    while len(edges_bip) < E:
        u = rng.randrange(L)
        v = rng.randrange(R)
        edges_bip.add((u, v))

    # Map to single vertex space: left u stays u, right v becomes L + v
    edges = sorted((u, L + v) for u, v in edges_bip)
    return edges


# ── File I/O ─────────────────────────────────────────────────────────

def save_general(filepath, V, edges):
    """Save in repo general format: 'V E' header, then 'u v' per edge."""
    with open(filepath, 'w') as f:
        f.write(f"{V} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")


def save_bipartite(filepath, L, R, edges):
    """Save in repo bipartite format: 'L R E' header, then 'l r' per edge."""
    with open(filepath, 'w') as f:
        f.write(f"{L} {R} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")


# ── Main logic ───────────────────────────────────────────────────────

def build_plan(sizes=None, graph_type=None):
    """Build list of (label, V, E, type, seed, filename) to generate."""
    plan = []
    for label, V, c in BENCHMARKS:
        # Filter by requested sizes
        if sizes is not None:
            size_tag = label.split("_")[0]  # "100k", "500k", "1m", etc.
            if size_tag not in sizes:
                continue

        E = c * V

        if graph_type is None or graph_type == "general":
            seed = SEED_BASE + V + c
            fname = f"general_sparse_{label}.txt"
            plan.append((label, V, E, "general", seed, fname))

        if graph_type is None or graph_type == "bipartite":
            seed = SEED_BASE + V + c + 1
            fname = f"bipartite_sparse_{label}.txt"
            plan.append((label, V, E, "bipartite", seed, fname))

        if graph_type is None or graph_type == "general_bipartite":
            edge_seed = SEED_BASE + V + c + 1  # same as bipartite
            shuffle_seed = SEED_BASE + V + c + 2  # separate seed for permutation
            fname = f"general_bipartite_sparse_{label}.txt"
            plan.append((label, V, E, "general_bipartite", (edge_seed, shuffle_seed), fname))

        if graph_type is None or graph_type == "stacked_bipartite":
            edge_seed = SEED_BASE + V + c + 1  # same as bipartite
            fname = f"stacked_bipartite_sparse_{label}.txt"
            plan.append((label, V, E, "stacked_bipartite", edge_seed, fname))

    return plan


def human_size(n):
    if n >= 1_000_000:
        return f"{n/1_000_000:.0f}M"
    elif n >= 1_000:
        return f"{n/1_000:.0f}K"
    return str(n)


def estimate_disk(E):
    """Rough estimate: ~16 bytes per edge as text."""
    return E * 16


def main():
    parser = argparse.ArgumentParser(
        description="Generate large-scale sparse benchmark graphs."
    )
    parser.add_argument(
        "--sizes", nargs="+",
        help="Generate only these sizes (e.g., 100k 1m 5m). Default: all."
    )
    parser.add_argument(
        "--type", choices=["general", "bipartite", "general_bipartite", "stacked_bipartite"],
        help="Generate only this type. Default: all four."
    )
    parser.add_argument(
        "--outdir", default="data/large-benchmarks",
        help="Output directory (default: data/large-benchmarks)."
    )
    parser.add_argument(
        "--list", action="store_true",
        help="Show generation plan without creating files."
    )
    args = parser.parse_args()

    plan = build_plan(sizes=args.sizes, graph_type=args.type)

    if not plan:
        print("No graphs match the requested filters.")
        sys.exit(1)

    # Show plan
    print(f"{'File':<40s} {'V':>10s} {'E':>10s} {'Type':<10s} {'~Disk':>10s}")
    print("-" * 85)
    total_disk = 0
    for label, V, E, gtype, seed, fname in plan:
        disk = estimate_disk(E)
        total_disk += disk
        print(f"{fname:<40s} {human_size(V):>10s} {human_size(E):>10s} "
              f"{gtype:<10s} {disk / 1e6:>8.0f} MB")
    print("-" * 85)
    print(f"{'Total':<40s} {'':>10s} {'':>10s} {'':>10s} {total_disk / 1e6:>8.0f} MB")
    print()

    if args.list:
        return

    # Generate
    os.makedirs(args.outdir, exist_ok=True)

    for i, (label, V, E, gtype, seed, fname) in enumerate(plan):
        filepath = os.path.join(args.outdir, fname)
        print(f"[{i+1}/{len(plan)}] Generating {fname} "
              f"(V={human_size(V)}, E={human_size(E)}, {gtype}) ...",
              end="", flush=True)

        t0 = time.time()

        if gtype == "general":
            edges = generate_general_sparse(V, E, seed)
            save_general(filepath, V, edges)
        elif gtype == "bipartite":
            L = V // 2
            R = V - L
            edges = generate_bipartite_sparse(L, R, E, seed)
            save_bipartite(filepath, L, R, edges)
        elif gtype == "general_bipartite":
            edge_seed, shuffle_seed = seed
            edges = generate_general_bipartite_sparse(V, E, edge_seed, shuffle_seed)
            save_general(filepath, V, edges)
        elif gtype == "stacked_bipartite":
            edges = generate_stacked_bipartite_sparse(V, E, seed)
            save_general(filepath, V, edges)

        t1 = time.time()
        size_mb = os.path.getsize(filepath) / 1e6
        print(f" done in {t1-t0:.1f}s ({size_mb:.0f} MB)")

    print(f"\nAll files written to {args.outdir}/")


if __name__ == "__main__":
    main()

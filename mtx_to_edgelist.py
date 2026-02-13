#!/usr/bin/env python3
"""Convert Matrix Market (.mtx) sparse matrix to edge-list format.

Input:  .mtx file (symmetric, from SuiteSparse Matrix Collection)
Output: our format — first line "V E", then E lines of "u v" (0-indexed, u < v)

Usage:
    python3 mtx_to_edgelist.py input.mtx output.txt
    python3 mtx_to_edgelist.py input.mtx              # writes to stdout
    python3 mtx_to_edgelist.py *.mtx --outdir graphs/  # batch mode

Self-loops are dropped (matching is on simple graphs).
Duplicate edges are dropped.
Vertices are renumbered to 0..V-1 if the matrix has 1-indexed rows/cols.
"""

import sys
import os
import argparse


def convert(mtx_path, out_file=None):
    """Convert a single .mtx file. Returns (V, E) on success."""
    edges = set()
    n = None
    header_read = False

    with open(mtx_path) as f:
        for line in f:
            line = line.strip()
            # Skip comments
            if line.startswith('%'):
                continue
            parts = line.split()
            # First non-comment line is the header: rows cols nnz
            if not header_read:
                n = int(parts[0])
                # parts[1] = cols (should equal rows for symmetric)
                # parts[2] = number of stored entries
                header_read = True
                continue
            # Data line: row col [value]
            i, j = int(parts[0]) - 1, int(parts[1]) - 1  # 0-indexed
            if i == j:
                continue  # skip self-loops
            u, v = min(i, j), max(i, j)
            edges.add((u, v))

    edges = sorted(edges)
    E = len(edges)

    # Write output
    close_after = False
    if out_file is None:
        f = sys.stdout
    elif isinstance(out_file, str):
        f = open(out_file, 'w')
        close_after = True
    else:
        f = out_file

    f.write(f"{n} {E}\n")
    for u, v in edges:
        f.write(f"{u} {v}\n")

    if close_after:
        f.close()

    name = os.path.basename(mtx_path)
    avg_deg = 2 * E / n if n > 0 else 0
    print(f"{name}: {n:,} vertices, {E:,} edges, "
          f"avg degree {avg_deg:.1f}", file=sys.stderr)
    return n, E


def main():
    parser = argparse.ArgumentParser(
        description="Convert Matrix Market (.mtx) to edge-list format")
    parser.add_argument('files', nargs='+', help='.mtx files to convert')
    parser.add_argument('--outdir', '-d', default=None,
                        help='Output directory (batch mode)')
    args = parser.parse_args()

    # Single file, second arg is output path
    if len(args.files) == 2 and not args.outdir:
        inp, out = args.files
        if not out.endswith('.mtx'):
            convert(inp, out)
            return

    # Batch mode or stdout
    for mtx_path in args.files:
        if not mtx_path.endswith('.mtx'):
            print(f"Skipping {mtx_path} (not .mtx)", file=sys.stderr)
            continue
        if args.outdir:
            os.makedirs(args.outdir, exist_ok=True)
            base = os.path.splitext(os.path.basename(mtx_path))[0]
            out_path = os.path.join(args.outdir, f"{base}.txt")
            convert(mtx_path, out_path)
        elif len(args.files) == 1:
            convert(mtx_path)  # stdout
        else:
            base = os.path.splitext(os.path.basename(mtx_path))[0]
            out_path = f"{base}.txt"
            convert(mtx_path, out_path)


if __name__ == '__main__':
    main()

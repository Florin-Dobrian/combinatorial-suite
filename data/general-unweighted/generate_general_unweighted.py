#!/usr/bin/env python3
"""
Data generator for general (non-bipartite) unweighted graphs.

This script generates random general graphs with controlled parameters
for benchmarking combinatorial algorithms.
"""

import random
import argparse
from pathlib import Path


def generate_random_graph(vertex_count, edge_probability, seed=None):
    """
    Generate a random general (undirected) graph.
    
    Args:
        vertex_count: Number of vertices
        edge_probability: Probability of edge between any two vertices (0-1)
        seed: Random seed for reproducibility
    
    Returns:
        List of (vertex1, vertex2) tuples
    """
    if seed is not None:
        random.seed(seed)
    
    edges = []
    for i in range(vertex_count):
        for j in range(i + 1, vertex_count):
            if random.random() < edge_probability:
                edges.append((i, j))
    
    return edges


def save_graph(filename, vertex_count, edges):
    """
    Save graph to file in standard format.
    
    Format:
        <vertex_count> <edge_count>
        <vertex1> <vertex2>
        ...
    """
    with open(filename, 'w') as f:
        f.write(f"{vertex_count} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")
    
    print(f"Generated {filename}:")
    print(f"  {vertex_count} vertices, {len(edges)} edges")
    max_edges = vertex_count * (vertex_count - 1) // 2
    density = len(edges) / max_edges if max_edges > 0 else 0
    print(f"  Edge density: {density:.2%}")


def main():
    parser = argparse.ArgumentParser(description='Generate general unweighted graph test data')
    parser.add_argument('--vertices', type=int, required=True, help='Number of vertices')
    parser.add_argument('--output', type=str, required=True, help='Output filename')
    parser.add_argument('--probability', type=float, required=True, help='Edge probability (0-1)')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    if not 0 <= args.probability <= 1:
        parser.error("Probability must be between 0 and 1")
    
    if args.vertices < 1:
        parser.error("Number of vertices must be at least 1")
    
    edges = generate_random_graph(args.vertices, args.probability, args.seed)
    save_graph(args.output, args.vertices, edges)


if __name__ == '__main__':
    main()

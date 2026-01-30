#!/usr/bin/env python3
"""
Data generator for bipartite matching test cases.

This script generates random bipartite graphs with controlled parameters
for benchmarking combinatorial algorithms.
"""

import random
import argparse
from pathlib import Path


def generate_random_graph(left_count, right_count, edge_probability, seed=None):
    """
    Generate a random bipartite graph.
    
    Args:
        left_count: Number of nodes in left partition
        right_count: Number of nodes in right partition
        edge_probability: Probability of edge between any left-right pair (0-1)
        seed: Random seed for reproducibility
    
    Returns:
        List of (left_node, right_node) tuples
    """
    if seed is not None:
        random.seed(seed)
    
    edges = []
    for i in range(left_count):
        for j in range(right_count):
            if random.random() < edge_probability:
                edges.append((i, j))
    
    return edges


def generate_regular_graph(left_count, right_count, degree):
    """
    Generate a random regular bipartite graph where each left node has exactly
    'degree' edges.
    
    Args:
        left_count: Number of nodes in left partition
        right_count: Number of nodes in right partition
        degree: Number of edges per left node
    
    Returns:
        List of (left_node, right_node) tuples
    """
    edges = []
    for i in range(left_count):
        # Sample 'degree' random right nodes
        right_nodes = random.sample(range(right_count), min(degree, right_count))
        for j in right_nodes:
            edges.append((i, j))
    
    return edges


def save_graph(filename, left_count, right_count, edges):
    """
    Save graph to file in standard format.
    
    Format:
        <left_count> <right_count> <edge_count>
        <left_node> <right_node>
        ...
    """
    with open(filename, 'w') as f:
        f.write(f"{left_count} {right_count} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")
    
    print(f"Generated {filename}:")
    print(f"  {left_count} left nodes, {right_count} right nodes, {len(edges)} edges")
    density = len(edges) / (left_count * right_count) if left_count * right_count > 0 else 0
    print(f"  Edge density: {density:.2%}")


def main():
    parser = argparse.ArgumentParser(description='Generate bipartite graph test data')
    parser.add_argument('--left', type=int, required=True, help='Number of left nodes')
    parser.add_argument('--right', type=int, required=True, help='Number of right nodes')
    parser.add_argument('--output', type=str, required=True, help='Output filename')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--probability', type=float, help='Edge probability (0-1)')
    group.add_argument('--degree', type=int, help='Edges per left node (regular graph)')
    
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')
    
    args = parser.parse_args()
    
    if args.probability is not None:
        if not 0 <= args.probability <= 1:
            parser.error("Probability must be between 0 and 1")
        edges = generate_random_graph(args.left, args.right, args.probability, args.seed)
    else:
        if args.degree < 0:
            parser.error("Degree must be non-negative")
        edges = generate_regular_graph(args.left, args.right, args.degree)
    
    save_graph(args.output, args.left, args.right, edges)


if __name__ == '__main__':
    main()

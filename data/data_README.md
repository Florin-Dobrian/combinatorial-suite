# Test Data and Benchmarks

This directory contains test cases and benchmark datasets for evaluating combinatorial algorithms.

## Overview

The data is organized by graph type and problem category:

- **Bipartite Unweighted Graphs** - For maximum cardinality matching in bipartite graphs (Hopcroft-Karp, etc.)
- **General Unweighted Graphs** - For maximum cardinality matching in general graphs (Edmonds' Blossom, Micali-Vazirani, etc.)

Future additions may include:
- Bipartite Weighted Graphs (Hungarian algorithm, minimum cost matching, etc.)
- General Weighted Graphs
- Directed Graphs
- Specialized structures (Trees, DAGs, etc.)

## Directory Structure

```
data/
├── data_README.md                       # This file
│
├── bipartite-unweighted/                # Bipartite unweighted graph data
│   ├── bipartite_unweighted_README.md   # Bipartite-specific documentation
│   ├── small/                           # Small test cases for correctness
│   ├── medium/                          # Medium benchmarks (100-1000 nodes)
│   ├── large/                           # Large benchmarks (10000+ nodes)
│   └── generate_bipartite_unweighted.py # Generator for bipartite graphs
│
└── general-unweighted/                  # General (non-bipartite) unweighted
    ├── general_unweighted_README.md     # General graph documentation
    ├── small/                           # Small test cases for correctness
    ├── medium/                          # Medium benchmarks (20-100 vertices)
    ├── large/                           # Large benchmarks (1000+ vertices)
    └── generate_general_unweighted.py   # Generator for general graphs
```

## Graph Type Comparison

### Bipartite Graphs

**Structure:**
- Vertices divided into two disjoint sets (left and right)
- Edges only connect vertices from different sets
- Contains only even-length cycles

**File Format:**
```
<left_count> <right_count> <edge_count>
<left_node> <right_node>
...
```

**Algorithms:**
- Hopcroft-Karp (O(E√V))
- Ford-Fulkerson on bipartite graphs

**Example:** Job assignment (workers to jobs), course scheduling

### General Graphs

**Structure:**
- No natural partition
- Can contain odd-length cycles (triangles, pentagons, etc.)
- More complex structure

**File Format:**
```
<vertex_count> <edge_count>
<vertex1> <vertex2>
...
```

**Algorithms:**
- Edmonds' Blossom (O(V⁴) simple, O(V²E) optimized)
- Micali-Vazirani (O(E√V))

**Example:** Social networks, molecular structures, general matching problems

## Quick Start

### For Bipartite Graphs
See [bipartite_unweighted_README.md](bipartite-unweighted/bipartite_unweighted_README.md) for:
- Detailed file format
- Hand-crafted test cases with visualizations
- Generation commands
- Usage examples in Python, C++, and Rust

### For General Graphs
See [general_unweighted_README.md](general-unweighted/general_unweighted_README.md) for:
- Detailed file format
- Classic graphs (Petersen, triangle, pentagon)
- Generation commands
- Usage examples in Python, C++, and Rust

## Summary Statistics

### Bipartite Unweighted Graphs
| Category | Files | Size Range | Purpose |
|----------|-------|------------|---------|
| Small | 4 | 3-5 nodes | Correctness verification |
| Medium | 4 | 100-1000 nodes | Performance testing |
| Large | 2 | 10000 nodes | Stress testing |

### General Unweighted Graphs
| Category | Files | Size Range | Purpose |
|----------|-------|------------|---------|
| Small | 4 | 3-10 vertices | Correctness, blossom testing |
| Medium | 3 | 20-100 vertices | Performance testing |
| Large | 2 | 1000-5000 vertices | Stress testing |

## Why Different Formats?

**Bipartite format** requires partition information (`<left_count> <right_count>`) because:
- The two-partition structure is fundamental to bipartite graphs
- Algorithms need to know which vertices belong to which partition
- Makes the graph type explicit

**General format** is simpler (`<vertex_count>`) because:
- No natural partition exists
- Vertices are all equivalent in structure
- Algorithms discover structure through search

This separation:
- Makes each format semantically correct for its graph type
- Prevents accidentally using wrong data with wrong algorithm
- Keeps files clean and purpose-specific

## Contributing Data

When adding new datasets:

1. Choose the appropriate graph type directory
2. Follow the format specification for that type
3. Use descriptive filenames with appropriate prefix
4. Place in correct size category (small/medium/large)
5. Update the relevant README with dataset details
6. Document expected results where applicable

## Data Sources

- **Small test cases**: Hand-crafted for algorithm verification and educational purposes
- **Medium/Large benchmarks**: Randomly generated with controlled parameters
- **Classic graphs**: Well-known graphs from graph theory literature (Petersen, etc.)

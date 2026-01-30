# General Unweighted Graph Data

This directory contains test cases and benchmark datasets for general (non-bipartite) unweighted graphs.

## Overview

General unweighted graphs can contain odd-length cycles, making them suitable for testing algorithms like:
- **Edmonds' Blossom Algorithm** - Maximum cardinality matching
- **Micali-Vazirani Algorithm** - Faster maximum cardinality matching  
- Other general graph algorithms

## Directory Structure

```
general-unweighted/
├── general_unweighted_README.md         # This file
├── small/                               # Small test cases for correctness
├── medium/                              # Medium benchmarks (20-100 vertices)
├── large/                               # Large benchmarks (1000+ vertices)
└── generate_general_unweighted.py       # Generator for general unweighted
                                         # graphs
```

## File Format

All general unweighted graph data files use a simple text format:

```
<vertex_count> <edge_count>
<vertex1> <vertex2>
<vertex1> <vertex2>
...
```

### Example: `general_unweighted_triangle.txt`
```
3 3
0 1
1 2
2 0
```

This represents:
- 3 vertices (0, 1, 2)
- 3 edges forming a triangle (odd cycle)
- Edges: (0-1), (1-2), (2-0)

### Vertex Naming Convention

- Vertices are represented as integers starting from 0
- Vertex range: 0 to (vertex_count - 1)
- In implementations, these are typically prefixed as "V0", "V1", ...

## General Unweighted Dataset Categories

### Small Test Cases (`general-unweighted/small/`)

Hand-crafted graphs for correctness testing:

#### `general_unweighted_triangle.txt` - Simple odd cycle (3 vertices)
```
Graph (Triangle - simplest odd cycle):

  0 ─── 1
   \   /
    \ /
     2

Maximum cardinality matching: 1 edge
Possible matchings: (0-1), (1-2), or (2-0)
```

#### `general_unweighted_pentagon.txt` - Pentagon (5-cycle)
```
Graph (Pentagon - 5-cycle):

    0
   / \
  4   1
  |   |
  3 - 2

Maximum cardinality matching: 2 edges
Example: (0-1), (3-4)
```

#### `general_unweighted_petersen.txt` - Petersen graph (classic example)
```
The Petersen graph is a well-known graph in graph theory:
- 10 vertices
- 15 edges
- 3-regular (each vertex has degree 3)
- Contains many odd cycles
- Maximum matching: 5 edges (perfect matching)

ASCII representation (outer pentagon + inner pentagram):

      0
     /|\
    / | \
   /  |  \
  4   |   1
  |\  |  /|
  | \ | / |
  |  \|/  |
  9   5   6
   \  |  /
    \ | /
     \|/
      7
      |
      3---2---8

Outer pentagon: 0-1-2-3-4-0
Inner pentagram (star): 5-7-9-6-8-5
Spokes connect outer to inner:
  0-5, 1-6, 2-7, 3-8, 4-9
```

#### `general_unweighted_blossom_example.txt` - Graph that creates a blossom
```
Graph designed to trigger blossom contraction:

  0 ─── 1 ─── 2
         \   /
          \ /
           3 ─── 4

Contains a triangle (1-2-3) that becomes a blossom during matching.
```

### Medium Benchmarks (`general-unweighted/medium/`)

Generated graphs for performance testing:

- **`general_unweighted_random_20.txt`** - 20 vertices, random edges
- **`general_unweighted_random_50.txt`** - 50 vertices, random edges  
- **`general_unweighted_dense_100.txt`** - 100 vertices, high edge density

### Large Benchmarks (`general-unweighted/large/`)

Large-scale datasets for stress testing:

- **`general_unweighted_random_1000.txt`** - 1000 vertices, moderate density
- **`general_unweighted_sparse_5000.txt`** - 5000 vertices, sparse
- **`general_unweighted_sparse_10000.txt`** - 10000 vertices, very sparse

## Generation Commands

All generated general unweighted datasets (excluding hand-crafted small test cases) were created using `generate_general_unweighted.py`:

```bash
# Medium benchmarks
python3 generate_general_unweighted.py --vertices 20 --probability 0.3 --seed 42 --output general_unweighted_random_20.txt

python3 generate_general_unweighted.py --vertices 50 --probability 0.2 --seed 42 --output general_unweighted_random_50.txt

python3 generate_general_unweighted.py --vertices 100 --probability 0.3 --seed 42 --output general_unweighted_dense_100.txt

# Large benchmarks
python3 generate_general_unweighted.py --vertices 1000 --probability 0.01 --seed 42 --output general_unweighted_random_1000.txt

python3 generate_general_unweighted.py --vertices 5000 --probability 0.001 --seed 42 --output general_unweighted_sparse_5000.txt

python3 generate_general_unweighted.py --vertices 10000 --probability 0.0005 --seed 42 --output general_unweighted_sparse_10000.txt
```

Or with `uv`:
```bash
uv run generate_general_unweighted.py --vertices 20 --probability 0.3 --seed 42 --output general_unweighted_random_20.txt
# ... (same commands, just replace python3 with "uv run")
```

**Notes:**
- `--seed 42` ensures reproducible random graphs
- `--probability` is the probability of an edge existing between any two vertices
- For general graphs, edge probability creates undirected edges

## Using Test Data

### Python Example
```python
def load_graph(filename):
    with open(filename, 'r') as f:
        vertex_count, edge_count = map(int, f.readline().split())
        
        vertices = [f"V{i}" for i in range(vertex_count)]
        edges = []
        
        for _ in range(edge_count):
            u, v = map(int, f.readline().split())
            edges.append((f"V{u}", f"V{v}"))
    
    return vertices, edges

# Usage
vertices, edges = load_graph("general_unweighted_triangle.txt")
```

## Expected Matching Sizes

For reference, here are the expected maximum matching sizes for the test cases:

### Small Test Cases
| File | Vertices | Edges | Max Matching | Notes |
|------|----------|-------|--------------|-------|
| general_unweighted_triangle.txt | 3 | 3 | 1 | Odd cycle |
| general_unweighted_pentagon.txt | 5 | 5 | 2 | 5-cycle |
| general_unweighted_petersen.txt | 10 | 15 | 5 | Perfect matching |
| general_unweighted_blossom_example.txt | 5 | 5 | 2 | Contains blossom |

### Generated Benchmarks
| File | Vertices | Edges | Edge Density |
|------|----------|-------|--------------|
| general_unweighted_random_20.txt | 20 | ~67 | ~35% |
| general_unweighted_random_50.txt | 50 | ~232 | ~19% |
| general_unweighted_dense_100.txt | 100 | ~1477 | ~30% |
| general_unweighted_random_1000.txt | 1000 | ~4985 | ~1% |
| general_unweighted_sparse_5000.txt | 5000 | ~12158 | ~0.1% |
| general_unweighted_sparse_10000.txt | 10000 | ~24907 | ~0.05% |

## Key Differences from Bipartite Graphs

**General graphs can have:**
- Odd-length cycles (triangles, pentagons, etc.)
- No natural partition into two sets
- More complex matching algorithms needed

**Bipartite graphs:**
- Only even-length cycles
- Clear left/right partition
- Simpler matching algorithms (Hopcroft-Karp)

## Contributing Data

When adding new test cases:

1. Follow the file format specification
2. Document the expected matching size
3. Name files descriptively with `general_unweighted_` prefix
4. Place in appropriate size category (small/medium/large)
5. Update this README with new entries

# Test Data and Benchmarks

This directory contains test cases and benchmark datasets for evaluating combinatorial algorithms.

## Overview

The data is organized by graph type and problem category. Currently available:
- **Bipartite Unweighted Graphs** - For maximum cardinality matching algorithms (Hopcroft-Karp, Ford-Fulkerson on bipartite graphs, etc.)

Future additions may include:
- Bipartite Weighted Graphs (for Hungarian algorithm, minimum cost matching, etc.)
- General (Non-bipartite) Graphs
- Directed Graphs
- Specialized structures (Trees, DAGs, etc.)

## Directory Structure

```
bipartite-unweighted/
├── bipartite_unweighted_README.md       # This file
├── small/                               # Small test cases for correctness
├── medium/                              # Medium benchmarks (100-1000 nodes)
├── large/                               # Large benchmarks (10000+ nodes)
└── generate_bipartite_unweighted.py     # Generator for bipartite unweighted
                                         # graphs
```

## Bipartite Unweighted Graphs

### File Naming Convention

Files follow the pattern: `bipartite_unweighted_<descriptor>.txt`

Examples:
- `bipartite_unweighted_simple.txt` - Hand-crafted simple case
- `bipartite_unweighted_random_100.txt` - Random graph, 100 nodes
- `bipartite_unweighted_dense_10000.txt` - Dense graph, 10000 nodes

## File Format (Bipartite Unweighted)

All bipartite unweighted graph data files use a simple text format:

```
<left_node_count> <right_node_count> <edge_count>
<left_node_1> <right_node_1>
<left_node_2> <right_node_2>
...
```

### Example: `bipartite_simple.txt`
```
3 3 4
0 0
0 1
1 1
2 2
```

This represents:
- 3 nodes on the left (0, 1, 2)
- 3 nodes on the right (0, 1, 2)
- 4 edges: (0→0), (0→1), (1→1), (2→2)

### Node Naming Convention

- Nodes are represented as integers starting from 0
- Left partition: 0 to (left_count - 1)
- Right partition: 0 to (right_count - 1)
- In implementations, these are typically prefixed as "L0", "L1", ... and "R0", "R1", ...

## Bipartite Unweighted Dataset Categories

### Small Test Cases (`bipartite-unweighted/small/`)

Hand-crafted graphs for correctness testing:

#### `bipartite_unweighted_simple.txt` - Basic matching with 3 nodes per side
```
Input Graph:                Maximum Matching (size 3):
Left    Right              Left    Right
 0 ───── 0                  0 ═════ 0
 │                          │
 └────── 1                  └────── 1
                            
 1 ───── 1                  1 ═════ 1
                            
 2 ───── 2                  2 ═════ 2

Legend: ═══ matched edge, ─── unmatched edge
```

#### `bipartite_unweighted_perfect.txt` - Perfect matching exists (all nodes matched)
```
Input Graph:                Maximum Matching (size 4):
Left    Right              Left    Right
 0 ───── 0                  0 ═════ 0
                            
 1 ───── 1                  1 ═════ 1
                            
 2 ───── 2                  2 ═════ 2
                            
 3 ───── 3                  3 ═════ 3

All nodes are matched (perfect matching).
```

#### `bipartite_unweighted_bottleneck.txt` - More left nodes than right (bottleneck)
```
Input Graph:                Maximum Matching (size 3):
Left    Right              Left    Right
 0 ───── 0                  0 ─────  0
 │                          │
 └────── 1                  └──────  1
                            
 1 ───── 1                  1 ═════  1
 │                          │
 └────── 2                  └────── 2
                            
 2 ───── 2                  2 ═════  2
                            
 3 ───── 0                  3 ═════  0
 │                          │
 └────── 1                  └────── 1
                            
 4 ───── 2                  4 ─────  2

Only 3 left nodes can be matched (limited by 3 right nodes).
Possible matching: (1,1), (2,2), (3,0) shown above.
Other valid maximum matchings exist.
```

#### `bipartite_unweighted_no_match.txt` - Cases where some nodes cannot be matched
```
Input Graph:                Maximum Matching (size 2):
Left    Right              Left    Right
 0 ───── 0                  0 ═════ 0
                            
 1 ───── 1                  1 ═════ 1
                            
 2       2                  2       2
         
Node 2 (left) has no edges, so it cannot be matched.
Node 2 (right) has no edges, so it cannot be matched.
Maximum matching size: 2 (out of 3 possible).
```

### Medium Benchmarks (`bipartite-unweighted/medium/`)

Generated graphs for performance testing:

- **`bipartite_unweighted_random_100.txt`** - 100 nodes per side, 5 edges per left node (500 total edges)
- **`bipartite_unweighted_dense_100.txt`** - 100 nodes per side, high edge density (~50%)
- **`bipartite_unweighted_random_500.txt`** - 500 nodes per side, ~2% edge density
- **`bipartite_unweighted_sparse_1000.txt`** - 1000 nodes per side, low edge density (~1%)

### Large Benchmarks (`bipartite-unweighted/large/`)

Large-scale datasets for stress testing:

- **`bipartite_unweighted_dense_10000.txt`** - 10,000 nodes per side, ~100k edges
- **`bipartite_unweighted_sparse_10000_degree5.txt`** - 10,000 nodes per side, 5 edges per left node (50k edges)

## Generation Commands

All generated bipartite unweighted datasets (excluding hand-crafted small test cases) were created using `generate_bipartite_unweighted.py`:

```bash
# Medium benchmarks
python3 generate_bipartite_unweighted.py --left 100 --right 100 --degree 5 --output bipartite_unweighted_random_100.txt

python3 generate_bipartite_unweighted.py --left 100 --right 100 --probability 0.5 --seed 42 --output bipartite_unweighted_dense_100.txt

python3 generate_bipartite_unweighted.py --left 500 --right 500 --probability 0.02 --seed 42 --output bipartite_unweighted_random_500.txt

python3 generate_bipartite_unweighted.py --left 1000 --right 1000 --probability 0.01 --seed 42 --output bipartite_unweighted_sparse_1000.txt

# Large benchmarks
python3 generate_bipartite_unweighted.py --left 10000 --right 10000 --probability 0.001 --seed 42 --output bipartite_unweighted_dense_10000.txt

python3 generate_bipartite_unweighted.py --left 10000 --right 10000 --degree 5 --output bipartite_unweighted_sparse_10000_degree5.txt
```

Or with `uv`:
```bash
uv run generate_bipartite_unweighted.py --left 100 --right 100 --degree 5 --output bipartite_unweighted_random_100.txt
# ... (same commands, just replace python3 with "uv run")
```

**Notes:**
- `--seed 42` ensures reproducible random graphs
- `--probability` creates random edges with given probability
- `--degree` creates regular graphs with fixed edges per left node

## Using Test Data

### Python Example
```python
def load_graph(filename):
    with open(filename, 'r') as f:
        left_count, right_count, edge_count = map(int, f.readline().split())
        
        left = [f"L{i}" for i in range(left_count)]
        right = [f"R{i}" for i in range(right_count)]
        edges = []
        
        for _ in range(edge_count):
            u, v = map(int, f.readline().split())
            edges.append((f"L{u}", f"R{v}"))
    
    return left, right, edges

# Usage
left, right, edges = load_graph("data/small/bipartite_simple.txt")
```

### C++ Example
```cpp
#include <fstream>
#include <sstream>

void load_graph(const string& filename, 
                vector<string>& left, 
                vector<string>& right,
                vector<pair<string, string>>& edges) {
    ifstream file(filename);
    int left_count, right_count, edge_count;
    file >> left_count >> right_count >> edge_count;
    
    for (int i = 0; i < left_count; i++) {
        left.push_back("L" + to_string(i));
    }
    for (int i = 0; i < right_count; i++) {
        right.push_back("R" + to_string(i));
    }
    
    for (int i = 0; i < edge_count; i++) {
        int u, v;
        file >> u >> v;
        edges.push_back({"L" + to_string(u), "R" + to_string(v)});
    }
}
```

### Rust Example
```rust
use std::fs::File;
use std::io::{BufRead, BufReader};

fn load_graph(filename: &str) -> (Vec<String>, Vec<String>, Vec<(String, String)>) {
    let file = File::open(filename).unwrap();
    let mut lines = BufReader::new(file).lines();
    
    let first_line = lines.next().unwrap().unwrap();
    let parts: Vec<usize> = first_line.split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    let (left_count, right_count, edge_count) = (parts[0], parts[1], parts[2]);
    
    let left: Vec<String> = (0..left_count).map(|i| format!("L{}", i)).collect();
    let right: Vec<String> = (0..right_count).map(|i| format!("R{}", i)).collect();
    
    let mut edges = Vec::new();
    for line in lines.take(edge_count) {
        let line = line.unwrap();
        let nums: Vec<usize> = line.split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        edges.push((format!("L{}", nums[0]), format!("R{}", nums[1])));
    }
    
    (left, right, edges)
}
```

## Generating Custom Datasets

You can generate custom test data using the provided scripts (coming soon) or manually create files following the format above.

### Random Graph Generation

For random bipartite graphs with controlled density:

```python
import random

def generate_random_graph(left_count, right_count, edge_probability, filename):
    edges = []
    for i in range(left_count):
        for j in range(right_count):
            if random.random() < edge_probability:
                edges.append((i, j))
    
    with open(filename, 'w') as f:
        f.write(f"{left_count} {right_count} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")

# Generate a sparse graph: 1000 nodes per side, 1% edge probability
generate_random_graph(1000, 1000, 0.01, "sparse_1000.txt")

# Generate a dense graph: 100 nodes per side, 50% edge probability
generate_random_graph(100, 100, 0.5, "dense_100.txt")
```

## Expected Matching Sizes

For reference, here are the expected maximum matching sizes for the test cases:

### Small Test Cases
| File | Left Nodes | Right Nodes | Edges | Max Matching |
|------|-----------|-------------|-------|--------------|
| bipartite_unweighted_simple.txt | 3 | 3 | 4 | 3 |
| bipartite_unweighted_perfect.txt | 4 | 4 | 4 | 4 |
| bipartite_unweighted_bottleneck.txt | 5 | 3 | 8 | 3 |
| bipartite_unweighted_no_match.txt | 3 | 3 | 2 | 2 |

### Generated Benchmarks
| File | Left Nodes | Right Nodes | Edges | Edge Density |
|------|-----------|-------------|-------|--------------|
| bipartite_unweighted_random_100.txt | 100 | 100 | 500 | 5.00% |
| bipartite_unweighted_dense_100.txt | 100 | 100 | ~4990 | ~49.90% |
| bipartite_unweighted_random_500.txt | 500 | 500 | ~4936 | ~1.97% |
| bipartite_unweighted_sparse_1000.txt | 1000 | 1000 | ~9962 | ~1.00% |
| bipartite_unweighted_dense_10000.txt | 10000 | 10000 | ~99887 | ~0.10% |
| bipartite_unweighted_sparse_10000_degree5.txt | 10000 | 10000 | 50000 | 0.05% |

## Contributing Data

When adding new test cases:

1. Follow the file format specification
2. Document the expected matching size
3. Name files descriptively
4. Place in appropriate size category (small/medium/large)
5. Update this README with new entries

## Data Sources

- Small test cases: Hand-crafted for algorithm verification
- Medium/Large benchmarks: Randomly generated with controlled parameters
- Real-world datasets: (To be added - sources will be cited here)

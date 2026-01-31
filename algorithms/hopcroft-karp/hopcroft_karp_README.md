# Hopcroft-Karp Algorithm

## Overview

The Hopcroft-Karp algorithm finds the maximum cardinality matching in a bipartite graph. It was discovered by John Hopcroft and Richard Karp in 1973 and represents a significant improvement over earlier augmenting path algorithms.

## Implementation Features

**✅ Integer Vertices Only**
- All vertices are integers (0, 1, 2, ...), not strings
- Better performance, better cache locality

**✅ Deterministic Behavior**
- C++: Uses `vector`, NOT `unordered_map`/`unordered_set`
- Rust: Uses `Vec`, NOT `HashMap`/`HashSet`  
- Python: Uses `list`, NOT `set()` for adjacency
- Sorted adjacency lists guarantee same output every time

**✅ Comprehensive Validation**
- Checks all edges exist in original graph
- Verifies no vertex appears twice
- Reports detailed errors if validation fails

## Input File Format

```
<left_count> <right_count> <edge_count>
<left_node> <right_node>
...
```

**Example:**
```
3 3 5
0 0
0 1
1 1
1 2
2 2
```

All vertices are integers starting from 0.

## Building and Running

### Python
```bash
python3 hopcroft_karp.py <filename>
```

Or with `uv` (faster, modern Python package manager):
```bash
uv run hopcroft_karp.py <filename>
```

### C++
```bash
g++ -O3 -std=c++17 hopcroft_karp.cpp -o hopcroft_karp
./hopcroft_karp <filename>
```

### Rust
```bash
rustc -O hopcroft_karp.rs -o hopcroft_karp
./hopcroft_karp <filename>
```

## Example Output

```
Hopcroft-Karp Algorithm - C++ Implementation
=============================================

Loading graph from: test.txt
File: test.txt
Graph: 1000 left nodes, 1000 right nodes, 10000 edges

=== Validation Report ===
Matching size (claimed): 995
Number of edges in matching: 995
Left nodes matched: 995
Right nodes matched: 995
VALIDATION PASSED: Matching is valid
=========================

Matching size: 995
Execution time: 7 ms
```

## Complexity

- **Time**: O(E√V)
- **Space**: O(V + E)

## Performance

Test: 1000 × 1000 nodes, 10,000 edges

| Language | Time | Speed |
|----------|------|-------|
| Python | ~10-15 ms | 1.0x |
| C++ | ~5-8 ms | ~2x |
| Rust | ~5-8 ms | ~2x |

## See Also

- Edmonds' Blossom for general (non-bipartite) graphs
- Hungarian algorithm for weighted matching

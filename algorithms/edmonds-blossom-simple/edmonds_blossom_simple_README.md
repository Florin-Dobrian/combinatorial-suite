# Edmonds' Blossom Algorithm (Simple)

## Overview

Edmonds' Blossom algorithm finds the maximum cardinality matching in general (non-bipartite) graphs by detecting and contracting odd cycles ("blossoms").

Discovered by Jack Edmonds in 1965, this was a breakthrough in combinatorial optimization.

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

## What's a Blossom?

Unlike bipartite graphs (only even cycles), general graphs have **odd cycles** that require special handling:

```
Triangle (odd cycle):
  0 --- 1
   \   /
    \ /
     2

The algorithm "contracts" this triangle into a single vertex,
finds paths in the contracted graph, then expands back.
```

## Input File Format

```
<vertex_count> <edge_count>
<vertex1> <vertex2>
...
```

**Example:**
```
5 7
0 1
1 2
2 0
2 3
3 4
```

All vertices are integers starting from 0.

## Building and Running

### Python
```bash
# Standard Python
python3 edmonds_blossom_simple.py <filename>

# Or with uv (faster, modern Python)
uv run edmonds_blossom_simple.py <filename>
```

### C++
```bash
g++ -O3 -std=c++17 edmonds_blossom_simple.cpp -o edmonds_blossom_simple_cpp
./edmonds_blossom_simple_cpp <filename>
```

### Rust
```bash
rustc -O edmonds_blossom_simple.rs -o edmonds_blossom_simple_rust
./edmonds_blossom_simple_rust <filename>
```

## Example Output

### Python
```
$ python3 edmonds_blossom_simple.py <filename>

Edmonds' Blossom Algorithm (Simple) - Python Implementation
============================================================

Loading graph from: <filename>
File: <filename>
Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size (claimed): 4962
Number of edges in matching: 4962
Number of unique vertices: 9924
VALIDATION PASSED: Matching is valid
=========================

Matching size: 4962
Execution time: 12647.77 ms
```

### C++
```
$ ./edmonds_blossom_simple_cpp <filename>

Edmonds' Blossom Algorithm (Simple) - C++ Implementation
============================================================

Loading graph from: <filename>
File: <filename>
Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size (claimed): 4962
Number of edges in matching: 4962
Number of unique vertices: 9924
VALIDATION PASSED: Matching is valid
=========================

Matching size: 4962
Execution time: 384 ms
```

### Rust
```
$ ./edmonds_blossom_simple_rust <filename>

Edmonds' Blossom Algorithm (Simple) - Rust Implementation
==============================================================

Loading graph from: <filename>
File: <filename>
Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size (claimed): 4962
Number of edges in matching: 4962
Number of unique vertices: 9924
VALIDATION PASSED: Matching is valid
=========================

Matching size: 4962
Execution time: 367 ms
```

## Complexity

### Simple Version
- **Time**: O(V⁴) 
- **Space**: O(V + E)

### Optimized Version
- **Time**: O(V²E) - see `edmonds_blossom_optimized`
- **Space**: O(V + E)

## Performance

Test: 10,000 vertices, 24,907 edges

| Language | Time | Speed |
|----------|------|-------|
| Python | ~12.6 sec | 1.0x |
| C++ | ~384 ms | ~33x |
| Rust | ~367 ms | ~34x |

**Note**: For large graphs (>1000 vertices), use the optimized version (O(V²E)) or Gabow/Micali-Vazirani (O(√VE)) for significantly better performance.

## When to Use

**Use Edmonds' Blossom when:**
- Graph has odd cycles (non-bipartite)
- Need maximum matching in general graphs
- Graph structure is arbitrary

**Use Hopcroft-Karp instead when:**
- Graph is bipartite (only even cycles)
- Much faster: O(√VE) vs O(V⁴)

## Comparison

| Algorithm | Complexity | Graph Type |
|-----------|-----------|------------|
| Hopcroft-Karp | O(√VE) | Bipartite only |
| Edmonds Simple | O(V⁴) | General graphs |
| Edmonds Optimized | O(V²E) | General graphs |

## See Also

- `edmonds_blossom_optimized` for O(V²E) version (faster for large graphs)
- `hopcroft_karp` for bipartite graphs (much faster)

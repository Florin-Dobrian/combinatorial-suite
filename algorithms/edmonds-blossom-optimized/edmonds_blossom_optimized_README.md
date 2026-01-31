# Edmonds' Blossom Algorithm (Optimized)

## Overview

This is the **O(V²E) optimized implementation** of Edmonds' Blossom algorithm for maximum cardinality matching in general graphs.

Significantly faster than the simple O(V⁴) version for large graphs.

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

## Complexity Improvement

| Version | Complexity | Graph Size |
|---------|-----------|------------|
| Simple | O(V⁴) | Good for <100 vertices |
| Optimized | O(V²E) | Good for all sizes |

**Practical speedup**: 5-10× faster on graphs with 1000+ vertices

## Key Optimizations

1. **Label-based search** - O(1) vertex classification (outer/inner/unlabeled)
2. **Efficient blossom contraction** - Updates only blossom vertices, not entire graph
3. **Path compression** - Fast base vertex lookups
4. **Smart queue management** - Avoids redundant work

## Input File Format

```
<vertex_count> <edge_count>
<vertex1> <vertex2>
...
```

Same as simple version.

## Building and Running

### Python
```bash
python3 edmonds_blossom_optimized.py <filename>
```

Or with `uv` (faster, modern Python package manager):
```bash
uv run edmonds_blossom_optimized.py <filename>
```

### C++
```bash
g++ -O3 -std=c++17 edmonds_blossom_optimized.cpp -o edmonds_blossom_optimized
./edmonds_blossom_optimized <filename>
```

### Rust
```bash
rustc -O edmonds_blossom_optimized.rs -o edmonds_blossom_optimized
./edmonds_blossom_optimized <filename>
```

## Example Output

```
Edmonds' Blossom Algorithm (Optimized) - C++ Implementation
============================================================

Loading graph from: test.txt
File: test.txt
Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size (claimed): 4962
Number of edges in matching: 4962
Number of unique vertices: 9924
VALIDATION PASSED: Matching is valid
=========================

Matching size: 4962
Execution time: 380 ms
```

## Performance Comparison

Test: 10,000 vertices, ~25,000 edges

### Simple vs Optimized (C++)

| Version | Time | Speedup |
|---------|------|---------|
| Simple | ~2-3 sec | 1.0x |
| Optimized | ~380-400 ms | **5-7× faster** |

### Language Comparison (Optimized)

| Language | Time | Speed |
|----------|------|-------|
| Python | ~5-10 sec | 1.0x |
| C++ | ~380-400 ms | ~15-20x |
| Rust | ~380-400 ms | ~15-20x |

## When to Use Optimized vs Simple

**Use Optimized:**
- Graphs with >1000 vertices
- Production use
- Performance matters
- Benchmarking

**Use Simple:**
- Learning the algorithm
- Small graphs (<100 vertices)
- Code readability priority
- Educational purposes

## Algorithm Details

The optimized version uses:
- **Labels**: Track vertex type (0=unlabeled, 1=outer, 2=inner) in O(1)
- **Base tracking**: Path compression for fast blossom base lookups
- **Incremental updates**: Only update affected vertices during blossom contraction

This reduces per-path complexity from O(V³) to O(VE), giving overall O(V²E).

## See Also

- `edmonds_blossom_simple` for O(V⁴) version (easier to understand)
- `hopcroft_karp` for bipartite graphs (even faster)

## References

- Edmonds, J. (1965). "Paths, trees, and flowers". *Canadian Journal of Mathematics*.
- Galil, Z. (1986). "Efficient algorithms for finding maximum matching in graphs". *ACM Computing Surveys*.

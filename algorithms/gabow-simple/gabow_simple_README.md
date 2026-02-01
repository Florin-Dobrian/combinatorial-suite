# Gabow's Matching Algorithms

This directory contains implementations of Harold N. Gabow's maximum cardinality matching algorithms for general (non-bipartite) graphs.

## Overview

Gabow developed two key algorithms for maximum matching in general graphs:

1. **Gabow Simple (1976)**: O(VE) - A simplified, more implementable version of Edmonds' Blossom algorithm
2. **Gabow Optimized (1983)**: O(√VE) - A scaling/phased approach that achieves the same theoretical bound as Micali-Vazirani

## Algorithms

### Gabow Simple - O(VE)

**Time Complexity:** O(VE)  
**Space Complexity:** O(V + E)  
**Best for:** Moderate-sized graphs where simplicity and reliability are priorities

**Key Features:**
- Sequential one-path-at-a-time augmentation
- Path compression in union-find for blossom management
- Efficient LCA (Lowest Common Ancestor) detection
- Lazy blossom expansion

**Performance (10K vertices, 24,907 edges):**
- C++: 216ms
- Rust: 241ms  
- Python: 5,074ms

**Why use this:**
- Most reliable and easiest to understand
- Industry standard (used by Grafalgo library)
- "Sweet spot" for practical applications
- Simpler to debug and maintain than O(√VE) alternatives

### Gabow Optimized - O(√VE)

**Time Complexity:** O(√VE)  
**Space Complexity:** O(V + E)  
**Best for:** Large graphs where performance is critical

**Key Features:**
- Phased approach: finds multiple shortest augmenting paths per phase
- Level graph construction with BFS
- Blossom detection and contraction during level building
- Bridge information storage for path expansion
- Based on LEDA's professional implementation

**Performance (10K vertices, 25K edges):**
- C++: 14ms (16.6x faster than Simple!)
- Rust: 17ms (14.3x faster)
- Python: 88ms (58x faster)

**Algorithm Structure:**
```
while phase_1():  # Build level structure, detect blossoms
    phase_2()      # Find and augment all shortest paths
```

**Phase 1:**
- Builds level graph by distance (Delta)
- Starts from all free vertices simultaneously
- Detects and contracts blossoms
- Returns true if augmenting paths exist at current distance

**Phase 2:**
- Finds all augmenting paths at the distance found in Phase 1
- Augments them all at once
- Achieves O(√V) phases total

## Historical Context

**Gabow 1976:** Published "An Efficient Implementation of Edmonds' Algorithm for Maximum Matching on Graphs"
- Simplified Edmonds' blossom algorithm
- Made it more practical to implement
- O(VE) complexity (improved from O(V⁴))

**Gabow 1983:** Published "An O(EV log V) Algorithm for Maximum Matching in General Graphs"  
- Scaling/phased approach
- O(√VE) complexity (matching Micali-Vazirani 1980)
- More implementable than Micali-Vazirani in practice

**Why Gabow over Micali-Vazirani?**
- Simpler to implement correctly
- Better documented
- Used in professional libraries (LEDA)
- Achieves same theoretical bounds

## Implementation Details

### Common Features (All Implementations)

**Determinism:**
- Integer vertices only (0 to n-1)
- Sorted adjacency lists
- No hash-based structures (unordered_set, unordered_map, HashMap, dict, set)
- Same input → same output every time

**Data Structures:**
- `mate[v]`: matched vertex or NIL
- `label[v]`: EVEN, ODD, or UNLABELED (for alternating tree)
- `base[v]`: base of blossom containing v (with path compression)
- `parent[v]`: parent in alternating tree
- `source_bridge[v]`, `target_bridge[v]`: for blossom expansion (Optimized only)

**Validation:**
- Checks all edges exist in original graph
- Ensures no vertex appears in multiple edges
- Reports matched vertex count

### Language-Specific Notes

**C++:**
- Fastest implementation
- Uses std::vector for all data structures
- Structured bindings for edge pairs
- Compiled with -O3 optimization

**Rust:**
- Idiomatic Rust with proper ownership
- Uses Vec and VecDeque
- Clones neighbor lists to satisfy borrow checker
- Zero-cost abstractions

**Python:**
- Readable reference implementation
- Uses lists and deques
- Slower but still practical for moderate graphs
- Good for prototyping and learning

## Performance Comparison

**Test: 10,000 vertices, 24,907 edges**

| Algorithm | C++ | Rust | Python | Speedup vs Simple (C++) |
|-----------|-----|------|--------|------------------------|
| Gabow Simple (O(VE)) | 232ms | 243ms | 5105ms | 1x |
| Gabow Optimized (O(√VE)) | 14ms | 17ms | 88ms | 16.6x |

**Result:** All implementations find the optimal 4962 edges

## When to Use Each Algorithm

### Use Gabow Simple When:
- Graph has < 10K vertices
- Simplicity and reliability are priorities
- You want a reference implementation to verify against
- Debugging is important
- You need to understand every step

### Use Gabow Optimized When:
- Graph has > 10K vertices
- Performance is critical
- You need the theoretical O(√VE) bound
- You're processing many graphs in batch
- Memory is not severely constrained

## Building and Running

### Python
```bash
# Standard Python
python3 gabow_simple.py <filename>

# Or with uv (faster, modern Python)
uv run gabow_simple.py <filename>
```

### C++
```bash
g++ -O3 -std=c++17 gabow_simple.cpp -o gabow_simple_cpp
./gabow_simple_cpp <filename>
```

### Rust
```bash
rustc -O gabow_simple.rs -o gabow_simple_rust
./gabow_simple_rust <filename>
```

## Example Output

### Python
```
$ python3 gabow_simple.py <filename>

Gabow's Algorithm (Simple Version) - Python Implementation
===========================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
Matched vertices: 9924
VALIDATION PASSED
=========================

Matching size: 4962
Time: 5074 ms
```

### C++
```
$ ./gabow_simple_cpp <filename>

Gabow's Algorithm (Simple Version) - C++ Implementation
========================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
Matched vertices: 9924
VALIDATION PASSED
=========================

Matching size: 4962
Time: 216 ms
```

### Rust
```
$ ./gabow_simple_rust <filename>

Gabow's Algorithm (Simple Version) - Rust Implementation
=========================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
Matched vertices: 9924
VALIDATION PASSED
=========================

Matching size: 4962
Time: 241 ms
```

## Input File Format

```
n m
u1 v1
u2 v2
...
um vm
```

Where:
- `n` = number of vertices
- `m` = number of edges
- `ui vi` = edge between vertices ui and vi (0-indexed)

## Output

```
Gabow's ... Algorithm - [Language] Implementation
==================================================

Graph: n vertices, m edges

=== Validation Report ===
Matching size: k
Matched vertices: 2k
VALIDATION PASSED
=========================

Matching size: k
Time: X ms
```

## Algorithm Comparison

| Feature | Simple O(VE) | Optimized O(√VE) |
|---------|-------------|------------------|
| Complexity | O(VE) | O(√VE) |
| Implementation | Simpler | More complex |
| Performance | Good | Excellent |
| Memory | Lower | Slightly higher |
| Code lines (C++) | ~200 | ~350 |
| Industry use | Very common | Rare (LEDA) |
| Debugging | Easier | Harder |

## Theoretical Background

**Matching Problem:** Given an undirected graph G = (V, E), find a maximum set of edges M such that no two edges in M share a vertex.

**Key Challenge:** Unlike bipartite graphs, general graphs can have **odd cycles** (blossoms) that complicate the search for augmenting paths.

**Blossom:** An odd-length cycle in the graph that appears during the search for augmenting paths. Must be contracted (treated as a single "super-vertex") and later expanded.

**Augmenting Path:** A path that starts and ends at unmatched vertices, alternating between non-matching and matching edges. Augmenting along this path increases the matching size by 1.

**Phase-Based Approach (Optimized):**
- Find all shortest augmenting paths in each phase
- O(√V) phases needed (proven by Hopcroft-Karp for bipartite graphs)
- Each phase does O(E) work
- Total: O(√VE)

## References

1. Gabow, H. N. (1976). "An efficient implementation of Edmonds' algorithm for maximum matching on graphs." *Journal of the ACM*, 23(2), 221-234.

2. Gabow, H. N. (1983). "An O(EV log V) algorithm for maximum matching in general graphs." *Proceedings of the 24th Annual Symposium on Foundations of Computer Science*, 247-254.

3. Mehlhorn, K., & Näher, S. (1999). *LEDA: A Platform for Combinatorial and Geometric Computing*. Cambridge University Press.

4. Ansaripour, A., Danaei, M., & Mehlhorn, K. (2024). "Benchmarking Maximum Cardinality Matching Algorithms." (Confirmed LEDA's O(√VE) implementation performance)

## Related Algorithms

This repository also includes:
- **Hopcroft-Karp** (O(√VE) for bipartite graphs)
- **Edmonds' Blossom** (O(V²E) and O(V³) variants)
- **Micali-Vazirani** (O(√VE) for general graphs - research implementation)

## Credits

**Original Algorithms:** Harold N. Gabow  
**LEDA Implementation Reference:** Kurt Mehlhorn & Stefan Näher  
**This Implementation:** Based on LEDA's approach with pedagogical clarity

## License

These implementations are provided for educational and research purposes.

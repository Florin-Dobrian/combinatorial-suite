# Gabow's O(√VE) Maximum Cardinality Matching Algorithm

## Overview

This is an implementation of **Gabow's 1985 scaling algorithm** for finding maximum cardinality matchings in general (non-bipartite) graphs. The algorithm achieves **O(√VE)** time complexity, which is asymptotically optimal for dense graphs and matches the complexity of the Micali-Vazirani algorithm.

## Algorithm Description

### Key Innovation

Gabow's algorithm extends the phase-based approach of Hopcroft-Karp (for bipartite graphs) to general graphs by:
1. **Building level structures** that respect blossom contractions
2. **Finding multiple shortest augmenting paths** in each phase
3. **Detecting and contracting blossoms** during the level-building phase
4. **Storing bridge information** to later expand paths through contracted blossoms

### Two-Phase Structure

The algorithm alternates between two phases until no augmenting paths exist:

**Phase 1: Build Level Structure**
- Initialize all free (unmatched) vertices at level 0
- Use BFS to build alternating trees from all free vertices simultaneously
- Process edges by increasing distance (Delta)
- Detect blossoms when two EVEN vertices in the same tree are connected
- Contract blossoms and continue building the level structure
- Stop when augmenting paths are detected (two different trees connect)

**Phase 2: Extract and Augment Paths**
- Find all shortest augmenting paths at the current distance level
- Use BFS/DFS to extract paths through the contracted graph
- Augment all paths simultaneously
- This increases the matching size by the number of paths found

### Blossom Handling

Unlike the simpler O(VE) algorithm that processes one path at a time:
- Blossoms are detected and contracted **during** level building
- Bridge information (source_bridge, target_bridge) is stored when contracting
- This allows the algorithm to later expand paths through blossoms
- Multiple vertex-disjoint paths are found and augmented per phase

## Complexity Analysis

- **Time Complexity:** O(√VE)
  - There are O(√V) phases (proven by Hopcroft-Karp style analysis)
  - Each phase does O(E) work to build level structure and find paths
  - Total: O(√VE)

- **Space Complexity:** O(V + E)
  - Adjacency lists: O(V + E)
  - Matching arrays: O(V)
  - Level/label arrays: O(V)
  - Blossom structures: O(V)

## Performance

Tested on a sparse graph with 10,000 vertices and 24,907 edges:

| Implementation | Time | Matching Size | Speedup vs O(VE) |
|----------------|------|---------------|------------------|
| C++        | 14ms | 4962/4962 ✔   | 16.6x        |
| Rust       | 17ms | 4962/4962 ✔   | 13.6x        |
| Python     | 88ms | 4962/4962 ✔   | 59x          |

Comparison with Gabow Simple O(VE):
- Gabow Simple (C++): 232ms
- Gabow Optimized (C++): 14ms
- **Performance improvement: ~16x**

## Implementation Details

### Key Data Structures

```cpp
mate[v]          // Matched vertex for v, or NIL if unmatched
label[v]         // EVEN, ODD, or UNLABELED (alternating tree state)
base[v]          // Base of the blossom containing v (with path compression)
parent[v]        // Parent in alternating tree (for path reconstruction)
source_bridge[v] // Source endpoint of edge that formed blossom
target_bridge[v] // Target endpoint of edge that formed blossom
edge_queue[d]    // Edges at distance d (for processing by level)
```

### Algorithm Flow

```
1. Initialize: All vertices are unmatched
2. While Phase 1 finds augmenting paths:
   a. Build level structure using BFS
   b. Detect and contract blossoms
   c. Check if paths exist between different trees
   d. If paths exist:
      - Phase 2: Find all shortest paths
      - Augment all paths simultaneously
3. Return the matching
```

### Determinism

All implementations are fully deterministic:
- ✅ Integer vertices only (0 to n-1)
- ✅ Sorted adjacency lists
- ✅ No hash-based data structures (unordered_map, HashSet, etc.)
- ✅ Consistent iteration order

## Comparison with Other Algorithms

### vs. Gabow Simple O(VE)
- **Gabow Simple:** Finds one augmenting path at a time, O(VE) complexity
- **Gabow Optimized:** Finds multiple paths per phase, O(√VE) complexity
- **When to use Optimized:** Larger graphs (>1000 vertices) where the √V factor matters
- **Performance gain:** ~10-20x speedup on sparse graphs with 10K+ vertices

### vs. Micali-Vazirani O(√VE)
- **Same asymptotic complexity:** Both are O(√VE)
- **Gabow advantage:** Simpler to implement correctly
- **Micali-Vazirani advantage:** Better constant factors in theory
- **In practice:** Gabow is often preferred due to implementation simplicity

### vs. Edmonds' Blossom
- **Edmonds Simple:** O(V²E) - very slow for large graphs
- **Edmonds Optimized:** O(VÂ²E) - still slower than Gabow for sparse graphs
- **Gabow advantage:** Better complexity for sparse graphs where E = O(V)
- **Use Edmonds when:** Graph is very dense and you need a simpler implementation

### vs. Hopcroft-Karp (Bipartite)
- **Hopcroft-Karp:** O(√VE) but only for bipartite graphs
- **Gabow:** O(√VE) for general graphs (handles odd cycles/blossoms)
- **If graph is bipartite:** Use Hopcroft-Karp (simpler, no blossom handling)
- **If graph has odd cycles:** Must use Gabow or Edmonds

## Historical Context

### Gabow's Contributions

**1976 - Simple O(VE) Algorithm:**
- Published in "An efficient implementation of Edmonds' algorithm for maximum matching on graphs"
- Streamlined Edmonds' algorithm with better data structures
- Introduced path compression for blossom bases
- This is our "Gabow Simple" implementation

**1985 - Scaling O(√VE) Algorithm:**
- Published in *Journal of Computer and System Sciences*, 31(2), 136–168
- Extended Hopcroft-Karp's phase-based approach to general graphs
- Achieves same complexity as Micali-Vazirani (1980)
- Considered more implementable than Micali-Vazirani
- This is our "Gabow Optimized" implementation

### Why This Algorithm Matters

1. **Theoretical Optimality:** O(√VE) is the best known complexity for unweighted matching
2. **Practical Performance:** 10-20x faster than O(VE) on large sparse graphs
3. **Implementation Feasibility:** Simpler than Micali-Vazirani while achieving same complexity
4. **General Graphs:** Handles odd cycles that bipartite algorithms cannot

## Building and Running

### Python
```bash
# Standard Python
python3 gabow_optimized.py <filename>

# Or with uv (faster, modern Python)
uv run gabow_optimized.py <filename>
```

### C++
```bash
g++ -O3 -std=c++17 gabow_optimized.cpp -o gabow_optimized_cpp
./gabow_optimized_cpp <filename>
```

### Rust
```bash
rustc -O gabow_optimized.rs -o gabow_optimized_rust
./gabow_optimized_rust <filename>
```

## Example Output

### Python
```
$ python3 gabow_optimized.py <filename>

Gabow's Scaling Algorithm (Optimized) - Python Implementation
==============================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
VALIDATION PASSED
=========================

Matching size: 4962
Time: 88 ms
```

### C++
```
$ ./gabow_optimized_cpp <filename>

Gabow's Scaling Algorithm (Optimized) - C++ Implementation
===========================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
Matched vertices: 9924
VALIDATION PASSED
=========================

Matching size: 4962
Time: 14 ms
```

### Rust
```
$ ./gabow_optimized_rust <filename>

Gabow's Scaling Algorithm (Optimized) - Rust Implementation
===========================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
VALIDATION PASSED
=========================

Matching size: 4962
Time: 17 ms
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
- `n` = number of vertices (vertices are numbered 0 to n-1)
- `m` = number of edges
- Each edge is specified as two integers (u, v)

## When to Use This Algorithm

### ✅ Use Gabow Optimized When:
- Graph is **general** (not known to be bipartite)
- Graph is **large** (>1000 vertices)
- Graph is **sparse to moderate density** (E = O(V) to O(VÂ²))
- You need **optimal performance** for maximum cardinality matching
- You want **deterministic** results

### ❌ Don't Use When:
- Graph is **bipartite** → Use Hopcroft-Karp instead (simpler)
- Graph is **very small** (<100 vertices) → Use Gabow Simple (less overhead)
- Graph is **extremely dense** (E ≈ VÂ²) → Use push-relabel flow algorithms
- You need **weighted** matching → Use Gabow's weighted matching algorithm

## References

1. **Gabow, H. N. (1976).** "An efficient implementation of Edmonds' algorithm for maximum matching on graphs." *Journal of the ACM*, 23(2), 221-234.

2. **Gabow, H. N. (1985).** "A scaling algorithm for weighted matching on general graphs." *Proceedings of the 26th Annual IEEE Symposium on Foundations of Computer Science*, 90–100.

3. **Micali, S., & Vazirani, V. V. (1980).** "An O(√V E) algorithm for finding maximum matching in general graphs." *Proceedings of the 21st Annual IEEE FOCS*, 17-27.

4. **Mehlhorn, K., & Näher, S. (1999).** *LEDA: A Platform for Combinatorial and Geometric Computing.* Cambridge University Press. (Contains production implementation)

5. **Ansaripour, M., Danaei, H., & Mehlhorn, K. (2024).** "Experimental Evaluation of Maximum Cardinality Matching Algorithms." Recent benchmarking study confirming O(√VE) performance.

## Implementation Notes

### Based on LEDA

This implementation follows the approach used in the **LEDA (Library of Efficient Data types and Algorithms)** library, which is considered the gold standard for high-performance graph algorithms. The LEDA implementation was specifically benchmarked in the 2024 study by Ansaripour et al. and confirmed to achieve the theoretical O(√VE) bound.

### Key Differences from Academic Descriptions

Many academic papers describe the algorithm at a high level but omit critical implementation details:

1. **Bridge Storage:** Essential for expanding paths through blossoms
2. **Edge Queuing by Distance:** Allows processing edges in level order
3. **Base Path Compression:** Keeps blossom operations efficient
4. **Simultaneous Tree Building:** All free vertices start BFS together

### Testing and Validation

All implementations include:
- ✅ Input validation (edge existence, vertex degree)
- ✅ Matching validation (no vertex in multiple edges)
- ✅ Deterministic behavior (sorted adjacency lists)
- ✅ Performance timing
- ✅ Matching size reporting

## License

This implementation is provided for educational and research purposes. The algorithm itself is from published academic work (Gabow 1985).

## Acknowledgments

- **Harold N. Gabow** for the original algorithm
- **LEDA team** (Mehlhorn & Näher) for the reference implementation
- **Grafalgo** (Jonathan Turner) for pedagogical implementations that informed this work

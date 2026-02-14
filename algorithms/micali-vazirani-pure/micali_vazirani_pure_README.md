# Micali-Vazirani Algorithm (Pure)

## Overview

This is a **faithful implementation** of the original Micali-Vazirani algorithm for maximum cardinality matching in general graphs with O(√VE) time complexity.

Unlike the hybrid Micali-Vazirani in this suite (which combines MV's level building with Gabow's path finding), this implementation uses the **original MV machinery**: Double Depth-First Search (DDFS), tenacity-based bridge classification, petal contraction, and the full min-level/max-level vertex structure described in the original 1980 paper and Vazirani's 1994 correctness proof.

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

## Why a "Pure" Version?

The hybrid Micali-Vazirani (`micali-vazirani/`) simplifies the algorithm by replacing the original DDFS path-finding with a BFS approach borrowed from Gabow. This makes the code shorter (~240 lines in C++) and easier to follow, but it is not the algorithm Micali and Vazirani actually described.

The pure version (~620 lines in C++) implements the full original algorithm:

| Feature | Hybrid MV | Pure MV |
|---------|-----------|---------|
| Level building (MIN phase) | MV-style BFS | MV-style BFS |
| Path finding (MAX phase) | Gabow-style BFS | Original DDFS |
| Bridge classification | Not used | Regular + hanging bridges |
| Tenacity | Not used | Full tenacity tracking |
| Petal contraction | Standard Edmonds | Original MV petal contraction |
| Code size (C++) | ~240 lines | ~620 lines |
| Performance (10K sparse) | 7 ms | 3 ms |

The pure version is both more faithful to the original paper and faster in practice, because the DDFS mechanism avoids redundant work that the simpler BFS approach performs.

## Algorithm Approach

The algorithm alternates between two phases until no augmenting paths remain:

**Phase 1 (MIN):**
- Builds a level structure from all free vertices simultaneously using BFS
- Assigns each vertex a min-level and max-level (even and odd levels)
- Detects bridges — edges connecting vertices at complementary levels
- Classifies bridges by tenacity (sum of the endpoints' levels)
- Processes bridges in order of increasing tenacity

**Phase 2 (MAX) — Double Depth-First Search:**
- For each unprocessed bridge, launches a DDFS
- The DDFS sends two searches (green and red) simultaneously toward free vertices
- Three outcomes are possible:
  - **PATH**: Both searches reach distinct free vertices → augmenting path found
  - **PETAL**: Both searches meet at the same blossom base → petal (blossom) detected, contract and continue
  - **EMPTY**: Search is blocked → no augmenting path through this bridge

**Key Data Structures:**
- `Node`: Tracks min/max/even/odd levels, predecessor lists, bud (blossom base), DDFS state
- `bridges[t]`: Bridges bucketed by tenacity `t`
- `hanging_bridges[v]`: Bridges deferred because vertex `v` was not yet at the right level

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
python3 micali_vazirani_pure.py <filename>

# Or with uv (faster, modern Python)
uv run micali_vazirani_pure.py <filename>
```

### C++
```bash
g++ -O3 -std=c++17 micali_vazirani_pure.cpp -o micali_vazirani_pure_cpp
./micali_vazirani_pure_cpp <filename>
```

### Rust
```bash
rustc -O micali_vazirani_pure.rs -o micali_vazirani_pure_rust
./micali_vazirani_pure_rust <filename>
```

## Example Output

### Python
```
$ python3 micali_vazirani_pure.py <filename>

Micali-Vazirani Pure Algorithm - Python Implementation
=======================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
Matched vertices: 9924
VALIDATION PASSED
=========================

Matching size: 4962
Time: 56 ms
```

### C++
```
$ ./micali_vazirani_pure_cpp <filename>

Micali-Vazirani Pure Algorithm - C++ Implementation
====================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
Matched vertices: 9924
VALIDATION PASSED
=========================

Matching size: 4962
Time: 3 ms
```

### Rust
```
$ ./micali_vazirani_pure_rust <filename>

Micali-Vazirani Pure Algorithm - Rust Implementation
=====================================================

Graph: 10000 vertices, 24907 edges

=== Validation Report ===
Matching size: 4962
Matched vertices: 9924
VALIDATION PASSED
=========================

Matching size: 4962
Time: 3 ms
```

## Complexity

- **Time**: O(√VE)
  - O(√V) phases
  - O(E) per phase (level building + DDFS)
- **Space**: O(V + E)

## Performance

**Test: 10,000 vertices, 24,907 edges (sparse general graph)**

### Correctness: 100% Optimal (4962/4962 edges)

| Implementation | Time | Matching Size | Status |
|----------------|------|---------------|--------|
| Python     | 56ms | 4962/4962 | ✅ 100% |
| C++        | 3ms  | 4962/4962 | ✅ 100% |
| Rust       | 3ms  | 4962/4962 | ✅ 100% |

### Comparison with Other O(√VE) Algorithms

| Algorithm | Time (C++) | Correctness |
|-----------|-----------|-------------|
| Gabow Optimized | 7ms | 4962/4962 ✔ |
| Micali-Vazirani (Hybrid) | 7ms | 4962/4962 ✔ |
| Micali-Vazirani (Pure) | 3ms | 4962/4962 ✔ |

The pure MV is the fastest algorithm in the suite — roughly 2× faster than the hybrid MV and Gabow Optimized on this benchmark.

### Comparison with Edmonds' Blossom

| Algorithm | Time (C++) | Matching | Speedup |
|-----------|-----------|----------|---------|
| Edmonds Simple | ~427 ms | 4962/4962 | 1.0x |
| Edmonds Optimized | ~413 ms | 4962/4962 | 1.0× |
| Micali-Vazirani (Pure) | 3 ms | 4962/4962 | ~140× |

## When to Use

**Use Micali-Vazirani Pure when:**
- Need maximum matching in general (non-bipartite) graphs
- Performance is critical — this is the fastest implementation in the suite
- Want a faithful implementation of the original MV algorithm
- Studying the DDFS / tenacity / bridge machinery from the original paper

**Use the Hybrid MV instead when:**
- Code simplicity is more important than absolute speed
- Want a shorter, more readable implementation (~240 vs ~620 lines in C++)
- The simpler BFS-based path finding is easier to modify for your use case

**Use Hopcroft-Karp instead when:**
- Graph is bipartite (no odd cycles)
- Simpler implementation preferred
- Same O(√VE) complexity but simpler for bipartite case

## Implementation Notes

This implementation is ported from the production-quality Jorants MV-Matching-V2 codebase. Key differences from the hybrid version:

1. **Double DFS (DDFS)**: Two simultaneous depth-first searches (green and red) trace paths from a bridge toward free vertices, detecting augmenting paths and petals
2. **Tenacity**: Bridges are classified by the sum of their endpoints' levels, ensuring they are processed in the correct order
3. **Regular vs. Hanging Bridges**: Bridges discovered during MIN are processed immediately; bridges deferred because a vertex had not yet received its max-level are stored as hanging bridges and processed later
4. **Petal Contraction**: When DDFS discovers a petal (odd cycle reachable from a single free vertex), the blossom is contracted using the original MV mechanism rather than standard Edmonds contraction

## References

- Micali, S., & Vazirani, V. V. (1980). "An O(√V E) algorithm for finding maximum matching in general graphs". *FOCS*.
- Vazirani, V. V. (1994). "A theory of alternating paths and blossoms for proving correctness of the O(√V E) general graph maximum matching algorithm". *Combinatorica*.
- Peterson, P. A., & Loui, M. C. (1988). "The general maximum matching algorithm of Micali and Vazirani". *Algorithmica*.

## See Also

- `micali-vazirani` for the hybrid O(√VE) version (simpler, shorter code)
- `gabow-optimized` for Gabow's O(√VE) approach (different algorithm, same complexity)
- `edmonds-blossom-optimized` for O(VÂ²E) version (simpler, still good for moderate graphs)
- `hopcroft-karp` for bipartite graphs (simpler when applicable)

# Edmonds' Blossom Algorithm

## Overview

Edmonds' Blossom algorithm finds the maximum cardinality matching in general (non-bipartite) graphs. It was discovered by Jack Edmonds in 1965 and represents a major breakthrough in combinatorial optimization by efficiently handling odd-length cycles (blossoms).

## What is Maximum Cardinality Matching in General Graphs?

A **general graph** can contain any structure, including odd-length cycles. A **matching** is a set of edges with no common vertices. A **maximum cardinality matching** is the largest possible matching.

Unlike bipartite graphs (which only have even cycles), general graphs require special handling of odd cycles called "blossoms."

### Example
```
Graph with triangle:
  A --- B
   \   /
    \ /
     C --- D --- E

Maximum cardinality matching: 2 edges
Possible matchings: (B-C, D-E) or (A-B, D-E) or (A-C, D-E)
```

## Algorithm Description

### Key Idea
The algorithm repeatedly finds **augmenting paths** while detecting and contracting **blossoms** (odd cycles) into super-vertices.

### Simple Version (O(V⁴))

This implementation finds one augmenting path per iteration:

1. **Search for augmenting path** from an unmatched vertex
   - Grow alternating tree (alternate matched/unmatched edges)
   - Detect blossoms when odd cycles are found
   - Contract blossoms into super-vertices
   
2. **When blossom found:**
   - Shrink the odd cycle to a single vertex
   - Continue search in the contracted graph
   
3. **When augmenting path found:**
   - Augment the matching along the path
   - Matching grows by 1
   
4. **Repeat** until no augmenting paths exist

### Why Blossoms Matter

In bipartite graphs, all cycles have even length, so simple augmenting path algorithms work. In general graphs, odd cycles create complications:

```
Triangle (odd cycle):
  0 --- 1
   \   /
    \ /
     2

If we match 0-1, then vertex 2 cannot improve the matching.
But if we "shrink" the triangle into one node, we can find better matchings.
```

The blossom contraction allows the algorithm to handle these odd cycles correctly.

## Complexity

### Simple Version
- **Time**: O(V⁴) where V = number of vertices
  - O(V) augmenting paths needed
  - O(V³) per path search (with naive blossom handling)
- **Space**: O(V + E)

### Optimized Version (Not Yet Implemented)
- **Time**: O(V²E) 
  - O(V) augmenting paths needed
  - O(VE) per path search (with efficient data structures)
- **Space**: O(V + E)

## When to Use

**Good for**:
- Finding maximum cardinality matchings in general graphs
- Graphs with odd-length cycles
- When you need the largest matching (not weighted)
- Molecular structure matching, social networks, general optimization

**Not ideal for**:
- Bipartite graphs (use Hopcroft-Karp instead - it's faster)
- Finding maximum *weight* matching (use different algorithm)
- When perfect matching detection is sufficient (use other techniques)

## Real-World Applications

1. **Molecular Chemistry**: Matching atoms in molecular structures
2. **Social Networks**: Pairing people with complementary interests
3. **Resource Allocation**: General assignment problems without bipartite structure
4. **Graph Theory Research**: Classic problem with theoretical importance
5. **Kidney Exchange**: Matching donors to recipients in complex networks

## Implementation Details

### Data Structures
- **Adjacency list**: For the general graph
- **Matching map**: Track paired vertices (mate)
- **Parent pointers**: For alternating tree structure
- **Base map**: Track which blossom each vertex belongs to

### Augmenting Path
An augmenting path alternates between:
- Unmatched edges (not in current matching)
- Matched edges (in current matching)

It starts at an unmatched vertex and ends at another unmatched vertex.

### Blossom
A blossom is an odd-length cycle in the alternating tree. When detected:
1. Contract all vertices in the cycle into a single "super-vertex"
2. Continue the search in the contracted graph
3. When path is found, expand blossoms to find actual path

## Building and Running

### Python
```bash
cd python/
python3 edmonds_blossom_simple.py
```

Or with `uv`:
```bash
uv run edmonds_blossom_simple.py
```

No compilation needed - Python is interpreted. No external dependencies required.

### C++
```bash
cd cpp/
g++ -O3 -std=c++17 edmonds_blossom_simple.cpp -o edmonds_blossom_simple_cpp
./edmonds_blossom_simple_cpp
```

Flags:
- `-O3`: Maximum optimization
- `-std=c++17`: C++17 standard
- `-o edmonds_blossom_simple_cpp`: Output executable name

On macOS, `g++` is an alias to `clang++` (Apple's compiler). Both work identically.

### Rust
```bash
cd rust/
rustc -O edmonds_blossom_simple.rs -o edmonds_blossom_simple_rust
./edmonds_blossom_simple_rust
```

Flags:
- `-O`: Enable optimizations (equivalent to `-O3` in C++)

### Running with Test Data Files

All implementations support loading graphs from data files:

**Python:**
```bash
cd python/
python3 edmonds_blossom_simple.py ../../../data/general-unweighted/small/general_unweighted_triangle.txt
```

Or with `uv`:
```bash
uv run edmonds_blossom_simple.py ../../../data/general-unweighted/small/general_unweighted_triangle.txt
```

**C++:**
```bash
cd cpp/
g++ -O3 -std=c++17 edmonds_blossom_simple.cpp -o edmonds_blossom_simple_cpp
./edmonds_blossom_simple_cpp ../../../data/general-unweighted/small/general_unweighted_triangle.txt
```

**Rust:**
```bash
cd rust/
rustc -O edmonds_blossom_simple.rs -o edmonds_blossom_simple_rust
./edmonds_blossom_simple_rust ../../../data/general-unweighted/small/general_unweighted_triangle.txt
```

**Try different test cases:**
```bash
# Small test cases for correctness
python3 edmonds_blossom_simple.py ../../../data/general-unweighted/small/general_unweighted_pentagon.txt
python3 edmonds_blossom_simple.py ../../../data/general-unweighted/small/general_unweighted_petersen.txt

# Medium benchmarks
./edmonds_blossom_simple_cpp ../../../data/general-unweighted/medium/general_unweighted_random_20.txt
./edmonds_blossom_simple_cpp ../../../data/general-unweighted/medium/general_unweighted_dense_100.txt

# Large benchmarks
./edmonds_blossom_simple_rust ../../../data/general-unweighted/large/general_unweighted_random_1000.txt
```

See the [general unweighted data README](../../data/general-unweighted/general_unweighted_README.md) for details on available datasets and their properties.

## Comparison with Other Algorithms

| Algorithm | Time Complexity | Graph Type | Use Case |
|-----------|----------------|------------|----------|
| Ford-Fulkerson | O(VE²) | Any | General max flow |
| Hopcroft-Karp | O(E√V) | Bipartite only | Bipartite matching (faster) |
| Edmonds Blossom (simple) | O(V⁴) | General | General matching |
| Edmonds Blossom (optimized) | O(V²E) | General | General matching (faster) |
| Micali-Vazirani | O(E√V) | General | General matching (fastest) |

## Key Differences from Hopcroft-Karp

**Hopcroft-Karp (bipartite only):**
- Finds MULTIPLE augmenting paths per iteration
- No odd cycles exist (bipartite property)
- Fewer iterations: O(√V)
- Simpler algorithm

**Edmonds' Blossom (general graphs):**
- Finds ONE augmenting path per iteration
- Must handle odd cycles (blossoms)
- More iterations: O(V)
- More complex but works on any graph

## Advanced Topics

### Blossom Detection
When two vertices in the alternating tree at even distance from the root are connected by an edge, an odd cycle (blossom) is formed. The base of the blossom is their lowest common ancestor.

### Blossom Contraction
All vertices in the blossom are merged into their base vertex. The search continues with this contracted graph.

### Theoretical Bounds
- **Perfect matching**: Exists if and only if |matching| = V/2
- **Tutte's theorem**: Characterizes when perfect matching exists
- **Gallai-Edmonds decomposition**: Structure theorem for general matching

## References

- Edmonds, J. (1965). "Paths, trees, and flowers". *Canadian Journal of Mathematics*, 17, 449-467.
- Cormen, T. H., et al. (2009). *Introduction to Algorithms* (3rd ed.). MIT Press.
- Kleinberg, J., & Tardos, É. (2005). *Algorithm Design*. Addison-Wesley.
- Lovász, L., & Plummer, M. D. (1986). *Matching Theory*. North-Holland.

## Implementation Notes

This directory contains the simple O(V⁴) version:
- **Python** (`simple/python/`): Clean, readable, good for understanding
- **C++** (`simple/cpp/`): High performance implementation
- **Rust** (`simple/rust/`): Memory-safe with good performance

Future implementations:
- **Optimized** (O(V²E)): Efficient data structures, better performance
- **Micali-Vazirani** (O(E√V)): Fastest known algorithm for general matching

# Hopcroft-Karp Algorithm

## Overview

The Hopcroft-Karp algorithm finds the maximum cardinality matching in a bipartite graph. It was discovered by John Hopcroft and Richard Karp in 1973 and represents a significant improvement over earlier augmenting path algorithms.

## What is a Bipartite Matching?

A **bipartite graph** consists of two disjoint sets of vertices (left and right), where edges only connect vertices from different sets. A **matching** is a set of edges with no common vertices. A **maximum matching** is the largest possible matching.

### Example
```
Left nodes: {A, B, C}     Right nodes: {1, 2, 3}

Edges:
  A --- 1
  A --- 2
  B --- 2
  C --- 3

Maximum matching: {(A,1), (B,2), (C,3)} - size 3
```

## Algorithm Description

### Key Idea
The algorithm repeatedly finds the **shortest augmenting paths** and augments the matching along multiple disjoint paths simultaneously in each phase.

### How It Works

1. **Initialization**: Start with an empty matching
2. **Repeat until no augmenting paths exist**:
   - **BFS Phase**: Find all shortest augmenting paths from unmatched left nodes
   - **DFS Phase**: Greedily find maximal set of vertex-disjoint augmenting paths and augment the matching along all of them

### Why It's Fast

Traditional augmenting path algorithms find one path at a time in O(E) per path, requiring O(V) iterations for O(VE) total complexity.

Hopcroft-Karp finds **all shortest paths** in one BFS, then augments along **multiple disjoint paths** simultaneously. The key insight: path lengths strictly increase, and there are at most O(√V) distinct path lengths.

## Complexity

- **Time**: O(E√V) where:
  - E = number of edges
  - V = number of vertices
- **Space**: O(V + E)

## When to Use

**Good for**:
- Finding maximum matchings in bipartite graphs
- Assignment problems (jobs to workers, tasks to machines)
- When you need the largest possible matching (not weighted)
- Graphs with up to millions of vertices and edges

**Not ideal for**:
- Finding maximum *weight* matching (use Hungarian algorithm instead)
- Non-bipartite graphs (use Edmonds' blossom algorithm)
- When you only need to know if a perfect matching exists (use Hall's theorem)

## Real-World Applications

1. **Job Assignment**: Match workers to jobs based on qualifications
2. **Dating Apps**: Match users with compatible partners
3. **Course Scheduling**: Assign students to courses with capacity constraints
4. **Organ Donation**: Match donors to recipients
5. **Network Routing**: Assign flows in bipartite networks
6. **Resource Allocation**: Distribute resources to requesters

## Implementation Details

### Data Structures
- **Adjacency list**: For the bipartite graph (left → right edges)
- **Matching maps**: Track paired vertices (pair_left, pair_right)
- **Distance map**: Store BFS distances for finding shortest paths

### Augmenting Path
An augmenting path alternates between:
- Unmatched edges (not in current matching)
- Matched edges (in current matching)

It starts at an unmatched left node and ends at an unmatched right node.

### Example Trace

```
Initial: empty matching
Graph: A-1, A-2, B-2, B-3, C-3

Phase 1:
  BFS: Find shortest paths from {A, B, C} → length 1 paths exist
  DFS: Augment A-1, B-2, C-3
  Matching: {(A,1), (B,2), (C,3)}

No more augmenting paths → DONE
Maximum matching size: 3
```

## Comparison with Other Algorithms

| Algorithm | Time Complexity | Use Case |
|-----------|----------------|----------|
| Ford-Fulkerson | O(VE²) | General max flow |
| Hopcroft-Karp | O(E√V) | **Bipartite matching** |
| Hungarian | O(V³) | **Weighted** bipartite matching |
| Edmonds Blossom | O(V²E) | General (non-bipartite) matching |

## Advanced Topics

### Optimizations
1. **Dense graphs**: Use matrix representation instead of adjacency lists
2. **Parallel processing**: Independent augmenting paths can be processed in parallel
3. **Integer node IDs**: Use integers instead of strings for better cache locality
4. **Early termination**: Stop when matching size reaches min(|left|, |right|)

### Variants
- **Minimum vertex cover**: Related via König's theorem
- **Maximum independent set**: In bipartite graphs
- **Edge coloring**: In bipartite graphs

### Theoretical Bounds
- **Perfect matching**: Exists if and only if |matching| = min(|left|, |right|)
- **Hall's theorem**: Characterizes when perfect matching exists
- **König's theorem**: max matching = min vertex cover (in bipartite graphs)

## References

- Hopcroft, J. E., & Karp, R. M. (1973). "An n^(5/2) algorithm for maximum matchings in bipartite graphs". *SIAM Journal on Computing*, 2(4), 225-231.
- Cormen, T. H., et al. (2009). *Introduction to Algorithms* (3rd ed.). MIT Press. Section 26.3.
- Kleinberg, J., & Tardos, É. (2005). *Algorithm Design*. Addison-Wesley. Section 7.5.

## Implementation Notes

This directory contains three implementations:
- **Python** (`python/`): Clean, readable, good for prototyping
- **C++** (`cpp/`): High performance, ~2× faster than Python
- **Rust** (`rust/`): Memory-safe with C++-like performance

See the top-level README for benchmarking results and usage examples.

## Building and Running

### Python
```bash
cd python/
python3 hopcroft_karp.py
```

Or with `uv`:
```bash
uv run hopcroft_karp.py
```

No compilation needed - Python is interpreted. No external dependencies required.

### C++
```bash
cd cpp/
g++ -O3 -std=c++17 hopcroft_karp.cpp -o hopcroft_karp_cpp
./hopcroft_karp_cpp
```

Flags:
- `-O3`: Maximum optimization
- `-std=c++17`: C++17 standard
- `-o hopcroft_karp_cpp`: Output executable name

On macOS, `g++` is an alias to `clang++` (Apple's compiler). Both work identically.

### Rust
```bash
cd rust/
rustc -O hopcroft_karp.rs -o hopcroft_karp_rust
./hopcroft_karp_rust
```

Flags:
- `-O`: Enable optimizations (equivalent to `-O3` in C++)

If you don't have Rust installed:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

### Running All Benchmarks
From the `benchmarks/` directory:
```bash
cd ../benchmarks/
chmod +x benchmark.sh
./benchmark.sh
```

This will compile and run all three implementations, showing performance comparisons.

## Usage Examples

### Python
```python
from hopcroft_karp import HopcroftKarp

left = ['A', 'B', 'C']
right = ['1', '2', '3']
edges = [('A', '1'), ('B', '2'), ('C', '3')]

hk = HopcroftKarp(left, right, edges)
matching = hk.maximum_matching()
print(f"Matching: {matching}")
```

### C++
```cpp
#include "hopcroft_karp.cpp"

std::vector<std::string> left = {"A", "B", "C"};
std::vector<std::string> right = {"1", "2", "3"};
std::vector<std::pair<std::string, std::string>> edges = {
    {"A", "1"}, {"B", "2"}, {"C", "3"}
};

HopcroftKarp hk(left, right, edges);
auto matching = hk.maximum_matching();
```

### Rust
```rust
let left: Vec<String> = vec!["A", "B", "C"]
    .into_iter().map(|s| s.to_string()).collect();
let right: Vec<String> = vec!["1", "2", "3"]
    .into_iter().map(|s| s.to_string()).collect();
let edges: Vec<(String, String)> = vec![
    ("A", "1"), ("B", "2"), ("C", "3")
].into_iter().map(|(a,b)| (a.to_string(), b.to_string())).collect();

let mut hk = HopcroftKarp::new(&left, &right, &edges);
let matching = hk.maximum_matching();
```

## Benchmark Results

Test configuration: 1000 left nodes, 1000 right nodes, 10,000 edges

| Language | Execution Time | Relative Speed |
|----------|---------------|----------------|
| Python   | 7.32 ms       | 1.0x (baseline)|
| C++      | 4.00 ms       | 1.83x faster   |
| Rust     | ~3-5 ms*      | ~1.5-2x faster*|

*Rust benchmark to be run locally (requires Rust toolchain)

### Performance Analysis

**Python**
- **Pros**: Easy to read, quick to prototype, excellent for learning
- **Cons**: Slower execution due to interpreter overhead
- **Use case**: Prototyping, educational purposes, small to medium graphs

**C++**
- **Pros**: Fast execution, mature optimizing compilers, good STL support
- **Cons**: Manual memory management (though STL helps), more verbose
- **Use case**: Production systems, large-scale graph processing

**Rust**
- **Pros**: Memory safety without garbage collection, zero-cost abstractions, modern tooling
- **Cons**: Steeper learning curve, longer compile times
- **Use case**: Systems where both performance and safety are critical

### Key Takeaways

1. **For prototyping**: Start with Python - it's the fastest to write and debug
2. **For performance**: C++ provides ~2x speedup over Python with mature optimization
3. **For safety + performance**: Rust offers comparable speed to C++ with memory safety guarantees
4. **Real-world usage**: The performance difference matters most for graphs with >10,000 nodes

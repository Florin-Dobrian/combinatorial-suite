# Combinatorial Algorithms Suite

A collection of high-performance implementations of fundamental combinatorial algorithms in Python, C++, and Rust.

## Overview

This repository provides well-documented, benchmarked implementations of classic combinatorial algorithms. Each algorithm is implemented in three languages to allow for:
- **Rapid prototyping** with Python
- **High performance** with C++ and Rust
- **Direct performance comparisons** across languages
- **Educational understanding** of algorithm behavior

## Implementation Features

All implementations follow best practices for performance and correctness:

**✅ Integer Vertices Only**
- All vertices are integers (0, 1, 2, ...), not strings
- Better performance and cache locality

**✅ Deterministic Behavior**
- C++: Uses `vector`, NOT `unordered_map`/`unordered_set`
- Rust: Uses `Vec`, NOT `HashMap`/`HashSet`
- Python: Uses `list`, NOT `set()` for graph adjacency
- Sorted adjacency lists guarantee reproducible results

**✅ Comprehensive Validation**
- All implementations validate output correctness
- Check edge existence, vertex uniqueness, degree constraints
- Detailed error reporting

## Current Algorithms

### Hopcroft-Karp Algorithm
Maximum cardinality bipartite matching in O(√VE) time.

**Location**: `algorithms/hopcroft-karp/`

**Implementations**:
- Python (clean, readable)
- C++ (optimized with -O3)
- Rust (memory-safe, high-performance)

See the [Hopcroft-Karp README](algorithms/hopcroft-karp/hopcroft_karp_README.md) for algorithm details, complexity analysis, and usage examples.

### Edmonds' Blossom Algorithm (Simple)
Maximum cardinality matching in general graphs in O(V⁴) time.

**Location**: `algorithms/edmonds-blossom-simple/`

**Implementations**:
- Python (clean, readable)
- C++ (optimized with -O3)
- Rust (memory-safe, high-performance)

See the [Edmonds' Blossom Simple README](algorithms/edmonds-blossom-simple/edmonds_blossom_simple_README.md) for algorithm details, blossom contraction, and usage examples.

### Edmonds' Blossom Algorithm (Optimized)
Maximum cardinality matching in general graphs in O(V²E) time.

**Location**: `algorithms/edmonds-blossom-optimized/`

**Implementations**:
- Python (clean, readable)
- C++ (optimized with -O3)
- Rust (memory-safe, high-performance)

**Performance**: 5-10× faster than simple version on large graphs (1000+ vertices).

See the [Edmonds' Blossom Optimized README](algorithms/edmonds-blossom-optimized/edmonds_blossom_optimized_README.md) for optimization details, complexity improvements, and performance comparisons.

### Gabow's Algorithm (Simple)
Maximum cardinality matching in general graphs in O(VE) time.

**Location**: `algorithms/gabow-simple/`

**Implementations**:
- Python (clean, readable)
- C++ (optimized with -O3)
- Rust (memory-safe, high-performance)

**Performance**: Industry standard implementation. Simpler and more reliable than Edmonds' optimized version while achieving better practical performance.

**Key Features**:
- Path compression for blossom bases
- Efficient LCA (Lowest Common Ancestor) detection
- Lazy blossom expansion
- Used in Grafalgo pedagogical library

See the [Gabow Simple README](algorithms/gabow-simple/gabow_simple_README.md) for algorithm details, implementation approach, and usage examples.

### Gabow's Algorithm (Optimized)
Maximum cardinality matching in general graphs in **O(√VE)** time - theoretically optimal!

**Location**: `algorithms/gabow-optimized/`

**Implementations**:
- Python (clean, readable)
- C++ (optimized with -O3)
- Rust (memory-safe, high-performance)

**Performance**: 10-20× faster than O(VE) on large sparse graphs (10,000+ vertices).

**Key Features**:
- Phased approach finding multiple augmenting paths per phase
- Achieves same O(√VE) bound as Micali-Vazirani
- Based on LEDA's production implementation
- Blossom contraction during level building

**Benchmarks** (10,000 vertices, 25,000 edges):
- C++: 14ms (16× faster than Simple)
- Rust: 17ms (14× faster than Simple)
- Python: 88ms (58× faster than Simple)

See the [Gabow Optimized README](algorithms/gabow-optimized/gabow_optimized_README.md) for detailed complexity analysis, LEDA-based approach, and performance comparisons.

### Micali-Vazirani Algorithm
Maximum cardinality matching in general graphs in **O(√VE)** time - theoretically optimal!

**Location**: `algorithms/micali-vazirani/`

**Implementations**:
- Python (clean, readable)
- C++ (optimized with -O3)
- Rust (memory-safe, high-performance)

**Performance**: Matches Gabow Optimized performance with 100% correctness.

**Key Features**:
- First algorithm to achieve O(√VE) for general graphs (1980)
- Hybrid approach combining MV's level building with Gabow's path finding
- Sophisticated even/odd level tracking
- Handles blossoms efficiently

**Benchmarks** (10,000 vertices, 24,907 edges):
- C++: 14ms (same as Gabow Optimized)
- Rust: 17ms (same as Gabow Optimized)
- Python: 80ms (slightly faster than Gabow Optimized)

**Result:** 4962/4962 edges (100% optimal) ✓

See the [Micali-Vazirani README](algorithms/micali-vazirani/micali_vazirani_README.md) for algorithm details, complexity analysis, and comparison with other O(√VE) approaches.

## Project Structure

```
combinatorial-suite/
├── README.md                            # This file
├── algorithms/
│   ├── hopcroft-karp/
│   │   ├── hopcroft_karp_README.md      # Algorithm-specific documentation
│   │   ├── python/hopcroft_karp.py
│   │   ├── cpp/hopcroft_karp.cpp
│   │   └── rust/hopcroft_karp.rs
│   ├── edmonds-blossom-simple/
│   │   ├── edmonds_blossom_simple_README.md  # Algorithm-specific documentation
│   │   ├── python/edmonds_blossom_simple.py
│   │   ├── cpp/edmonds_blossom_simple.cpp
│   │   └── rust/edmonds_blossom_simple.rs
│   ├── edmonds-blossom-optimized/
│   │   ├── edmonds_blossom_optimized_README.md  # Algorithm-specific documentation
│   │   ├── python/edmonds_blossom_optimized.py
│   │   ├── cpp/edmonds_blossom_optimized.cpp
│   │   └── rust/edmonds_blossom_optimized.rs
│   ├── gabow-simple/
│   │   ├── gabow_simple_README.md       # Algorithm-specific documentation
│   │   ├── python/gabow_simple.py
│   │   ├── cpp/gabow_simple.cpp
│   │   └── rust/gabow_simple.rs
│   ├── gabow-optimized/
│   │   ├── gabow_optimized_README.md    # Algorithm-specific documentation
│   │   ├── python/gabow_optimized.py
│   │   ├── cpp/gabow_optimized.cpp
│   │   └── rust/gabow_optimized.rs
│   └── micali-vazirani/
│       ├── micali_vazirani_README.md    # Algorithm-specific documentation
│       ├── python/micali_vazirani.py
│       ├── cpp/micali_vazirani.cpp
│       └── rust/micali_vazirani.rs
├── benchmarks/
│   └── benchmark.sh                     # Cross-language performance testing
└── data/                                # Test data and datasets
    ├── data_README.md                   # Data format documentation
    ├── bipartite-unweighted/            # Bipartite unweighted graph data
    │   ├── bipartite_unweighted_README.md
    │   ├── small/                       # Small test cases for correctness
    │   ├── medium/                      # Medium benchmarks (100-1000 nodes)
    │   ├── large/                       # Large benchmarks (10000+ nodes)
    │   └── generate_bipartite_unweighted.py
    └── general-unweighted/              # General (non-bipartite) unweighted
        ├── general_unweighted_README.md # graph data
        ├── small/                       # Small test cases for correctness
        ├── medium/                      # Medium benchmarks (20-100 vertices)
        ├── large/                       # Large benchmarks (1000+ vertices)
        └── generate_general_unweighted.py
```

## Quick Start

### Prerequisites
- **Python 3.14+** (via `uv` or system installation)
- **C++ compiler** (g++ or clang++)
- **Rust toolchain** (install via [rustup](https://rustup.rs/))

### Running an Algorithm

#### With Python
```bash
cd algorithms/hopcroft-karp/python/
python3 hopcroft_karp.py <datafile>

# Or with uv (faster, modern Python):
uv run hopcroft_karp.py <datafile>
```

#### With C++
```bash
cd algorithms/hopcroft-karp/cpp/
g++ -O3 -std=c++17 hopcroft_karp.cpp -o hopcroft_karp_cpp
./hopcroft_karp_cpp <datafile>
```

#### With Rust
```bash
cd algorithms/hopcroft-karp/rust/
rustc -O hopcroft_karp.rs -o hopcroft_karp_rust
./hopcroft_karp_rust <datafile>
```

### Example with Test Data

```bash
# Hopcroft-Karp on bipartite graph
cd algorithms/hopcroft-karp/cpp/
g++ -O3 -std=c++17 hopcroft_karp.cpp -o hopcroft_karp_cpp
./hopcroft_karp_cpp ../../../data/bipartite-unweighted/large/bipartite_unweighted_dense_10000.txt

# Edmonds Blossom (optimized) on general graph
cd algorithms/edmonds-blossom-optimized/cpp/
g++ -O3 -std=c++17 edmonds_blossom_optimized.cpp -o edmonds_blossom_optimized_cpp
./edmonds_blossom_optimized_cpp ../../../data/general-unweighted/large/general_unweighted_sparse_10000.txt

# Gabow Optimized (O(√VE) - theoretically optimal!) on general graph
cd algorithms/gabow-optimized/cpp/
g++ -O3 -std=c++17 gabow_optimized.cpp -o gabow_optimized_cpp
./gabow_optimized_cpp ../../../data/general-unweighted/large/general_unweighted_sparse_10000.txt

# Micali-Vazirani (O(√VE) - theoretically optimal!) on general graph
cd algorithms/micali-vazirani/cpp/
g++ -O3 -std=c++17 micali_vazirani.cpp -o micali_vazirani_cpp
./micali_vazirani_cpp ../../../data/general-unweighted/large/general_unweighted_sparse_10000.txt
```

### Running Benchmarks

```bash
cd benchmarks/
chmod +x benchmark.sh
./benchmark.sh
```

## Performance Comparison

**Test Hardware:** MacBook Pro (November 2024) with M4 processor (10 cores) and 32GB memory

**Test Graph:** 10,000 vertices, 24,907 edges (sparse general graph)

### Bipartite Matching

**Hopcroft-Karp** (10,000 nodes, bipartite):
- Python: ~10-15 ms
- C++: ~5-8 ms
- Rust: ~5-8 ms

### General Graph Matching

**Edmonds Blossom Simple** O(V⁴):
- Python: ~12,648 ms (12.6 seconds)
- C++: ~384 ms
- Rust: ~367 ms

**Edmonds Blossom Optimized** O(V²E):
- Python: ~12,616 ms (12.6 seconds)
- C++: ~383 ms
- Rust: ~361 ms

**Gabow Simple** O(VE):
- Python: ~5,074 ms (5.1 seconds)
- C++: ~216 ms
- Rust: ~241 ms

**Gabow Optimized** O(√VE):
- Python: ~88 ms (58× faster than Gabow Simple!)
- C++: ~14 ms (15.4× faster than Gabow Simple!)
- Rust: ~17 ms (14.2× faster than Gabow Simple!)

**Micali-Vazirani** O(√VE):
- Python: ~80 ms (63× faster than Gabow Simple!)
- C++: ~14 ms (15.4× faster than Gabow Simple!)
- Rust: ~17 ms (14.2× faster than Gabow Simple!)

### Algorithm Complexity Comparison (General Graphs)

| Algorithm | Complexity | Python Time | C++ Time | Rust Time | Speedup vs Edmonds Simple (C++) |
|-----------|-----------|-------------|----------|-----------|--------------------------------|
| Edmonds Simple | O(V⁴) | 12.6 sec | 384 ms | 367 ms | 1.0× |
| Edmonds Optimized | O(V²E) | 12.6 sec | 383 ms | 361 ms | 1.0× |
| Gabow Simple | O(VE) | 5.1 sec | 216 ms | 241 ms | 1.8× |
| Gabow Optimized | O(√VE) | 88 ms | 14 ms | 17 ms | 27.4× |
| Micali-Vazirani | O(√VE) | 80 ms | 14 ms | 17 ms | 27.4× |

**Key Insights:**
- Both O(√VE) algorithms (Gabow Optimized and Micali-Vazirani) achieve 100% correctness (4962/4962 edges)
- O(√VE) algorithms are ~27× faster than O(VE) and ~27× faster than O(V²E)
- For this sparse graph, Edmonds Simple and Optimized have similar performance
- Python is 5-6× slower than C++/Rust for complex algorithms

## Git Setup

### Initialize Repository

From the `combinatorial-suite/` directory:

```bash
# Initialize git
git init

# Create .gitignore
cat > .gitignore << EOF
# Python
__pycache__/
*.py[cod]
.venv/
*.egg-info/

# C++
*.o
*_cpp
hopcroft_karp_cpp
edmonds_blossom_simple_cpp
edmonds_blossom_optimized_cpp
gabow_simple_cpp
gabow_optimized_cpp
micali_vazirani_cpp

# Rust
*_rust
*.rlib

# macOS
.DS_Store

# IDE
.vscode/
.idea/
EOF

# Stage and commit
git add .
git commit -m "Initial commit: Combinatorial algorithms suite with matching algorithms"
```

### Push to GitHub

#### Option 1: Using GitHub CLI (recommended)

**For a public repository:**
```bash
gh repo create combinatorial-suite --public --license apache-2.0
git remote add origin https://github.com/YOUR_USERNAME/combinatorial-suite.git
git pull origin main --allow-unrelated-histories --no-rebase  # If LICENSE was created
git push -u origin main
```

**For a private repository (share with collaborators only):**
```bash
gh repo create combinatorial-suite --private --license apache-2.0
git remote add origin https://github.com/YOUR_USERNAME/combinatorial-suite.git
git pull origin main --allow-unrelated-histories --no-rebase  # If LICENSE was created
git push -u origin main
```

**Note:** Replace `YOUR_USERNAME` with your actual GitHub username. The pull command is needed if GitHub created a LICENSE file that isn't in your local repo.

#### Option 2: Manual Setup

```bash
# Create repo manually at github.com/new (choose public or private)
git remote add origin https://github.com/YOUR_USERNAME/combinatorial-suite.git
git push -u origin main
```

### Adding Collaborators (Private Repositories)

If you created a private repository, you can add collaborators to give them access.

#### Via GitHub CLI:
```bash
# Add a collaborator
gh repo add-collaborator COLLABORATOR_USERNAME

# Example:
gh repo add-collaborator john_doe
```

#### Via GitHub Web Interface:
1. Go to your repository: `github.com/YOUR_USERNAME/combinatorial-suite`
2. Click **Settings** (top right)
3. Click **Collaborators** in the left sidebar
4. Click **Add people**
5. Enter their GitHub username or email
6. Select the person and choose permission level (Read, Write, or Admin)
7. Click **Add [username] to this repository**

They'll receive an email invitation. Once accepted, they can access the repository.

### Making a Private Repository Public

You can change visibility at any time:

#### Via GitHub CLI:
```bash
gh repo edit combinatorial-suite --visibility public
```

#### Via GitHub Web Interface:
1. Go to your repository on GitHub
2. Click **Settings**
3. Scroll to the **Danger Zone** section (bottom of page)
4. Click **Change visibility**
5. Select **Make public**
6. Confirm by typing the repository name
7. Click **I understand, make this repository public**

**Note:** Once public, anyone can see all code and commit history. You can change back to private if needed.

### Repository Visibility Summary

| Visibility | Who Can See | Use Case |
|------------|-------------|----------|
| Public | Everyone on the internet | Open source, portfolios, sharing with community |
| Private | Only you and invited collaborators | Work in progress, proprietary code, selective sharing |

The `.gitignore` file excludes:
- Compiled binaries (C++ and Rust executables)
- Python virtual environments and cache files
- OS-specific files (.DS_Store on macOS)
- IDE configuration directories

## License

Apache License 2.0 - Feel free to use and modify as needed.

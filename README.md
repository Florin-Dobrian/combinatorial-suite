# Combinatorial Algorithms Suite

A collection of high-performance implementations of fundamental combinatorial algorithms in Python, C++, and Rust.

## Overview

This repository provides well-documented, benchmarked implementations of classic combinatorial algorithms. Each algorithm is implemented in three languages to allow for:
- **Rapid prototyping** with Python
- **High performance** with C++ and Rust
- **Direct performance comparisons** across languages
- **Educational understanding** of algorithm behavior

## Current Algorithms

### Hopcroft-Karp Algorithm
Maximum cardinality bipartite matching in O(E√V) time.

**Location**: `algorithms/hopcroft-karp/`

**Implementations**:
- Python (clean, readable)
- C++ (optimized with -O3)
- Rust (memory-safe, high-performance)

See the [Hopcroft-Karp README](algorithms/hopcroft-karp/README.md) for algorithm details, complexity analysis, and usage examples.

## Project Structure

```
combinatorial-suite/
├── README.md                    # This file
├── algorithms/
│   └── hopcroft-karp/
│       ├── README.md            # Algorithm-specific documentation
│       ├── python/hopcroft_karp.py
│       ├── cpp/hopcroft_karp.cpp
│       └── rust/hopcroft_karp.rs
├── benchmarks/
│   └── benchmark.sh             # Cross-language performance testing
└── data/                        # Test data and datasets
```

## Quick Start

### Prerequisites
- **Python 3.14+** (via `uv` or system installation)
- **C++ compiler** (g++ or clang++)
- **Rust toolchain** (install via [rustup](https://rustup.rs/))

### Running an Algorithm

```bash
# Python
cd algorithms/hopcroft-karp/python/
python3 hopcroft_karp.py

# C++
cd algorithms/hopcroft-karp/cpp/
g++ -O3 -std=c++17 hopcroft_karp.cpp -o hopcroft_karp_cpp
./hopcroft_karp_cpp

# Rust
cd algorithms/hopcroft-karp/rust/
rustc -O hopcroft_karp.rs -o hopcroft_karp_rust
./hopcroft_karp_rust
```

### Running Benchmarks

```bash
cd benchmarks/
chmod +x benchmark.sh
./benchmark.sh
```

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
hopcroft_karp_cpp

# Rust
hopcroft_karp_rust
*.rlib

# macOS
.DS_Store

# IDE
.vscode/
.idea/
EOF

# Stage and commit
git add .
git commit -m "Initial commit: Hopcroft-Karp algorithm in Python, C++, and Rust"
```

### Push to GitHub

**Option 1: Using GitHub CLI (recommended if installed)**
```bash
gh repo create combinatorial-suite --public --license apache-2.0 --source=. --push
```

**Option 2: Manual setup**
```bash
# Create repo manually at github.com/new
git remote add origin git@github.com:YOUR_USERNAME/combinatorial-suite.git
git branch -M main
git push -u origin main
```

The `.gitignore` file excludes:
- Compiled binaries (C++ and Rust executables)
- Python virtual environments and cache files
- OS-specific files (.DS_Store on macOS)
- IDE configuration directories

## License

Apache License 2.0 - Feel free to use and modify as needed.

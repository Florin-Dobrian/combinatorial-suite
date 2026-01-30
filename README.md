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

See the [Hopcroft-Karp README](algorithms/hopcroft-karp/hopcroft_karp_README.md) for algorithm details, complexity analysis, and usage examples.

### Edmonds' Blossom Algorithm (Simple)
Maximum cardinality matching in general graphs in O(V⁴) time.

**Location**: `algorithms/edmonds-blossom/simple/`

**Implementations**:
- Python (clean, readable)
- C++ (optimized with -O3)
- Rust (memory-safe, high-performance)

See the [Edmonds' Blossom Simple README](algorithms/edmonds-blossom/simple/edmonds_blossom_simple_README.md) for algorithm details, blossom contraction, and usage examples.

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
│   └── edmonds-blossom/
│       └── simple/
│           ├── edmonds_blossom_simple_README.md  # Algorithm-specific documentation
│           ├── python/edmonds_blossom_simple.py
│           ├── cpp/edmonds_blossom_simple.cpp
│           └── rust/edmonds_blossom_simple.rs
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
git commit -m "Initial commit: Combinatorial algorithms suite with Hopcroft-Karp (Python, C++, Rust)"
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
| **Public** | Everyone on the internet | Open source, portfolios, sharing with community |
| **Private** | Only you and invited collaborators | Work in progress, proprietary code, selective sharing |

The `.gitignore` file excludes:
- Compiled binaries (C++ and Rust executables)
- Python virtual environments and cache files
- OS-specific files (.DS_Store on macOS)
- IDE configuration directories

## License

Apache License 2.0 - Feel free to use and modify as needed.

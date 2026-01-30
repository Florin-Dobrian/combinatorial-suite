#!/bin/bash

# Benchmark script for Hopcroft-Karp implementations
# Run from: combinatorial-suite/benchmarks/

echo "========================================="
echo "Hopcroft-Karp Algorithm Benchmark"
echo "Comparing Python, C++, and Rust"
echo "========================================="
echo ""

# Navigate to the algorithm directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ALGO_DIR="$SCRIPT_DIR/../algorithms/hopcroft-karp"

# Compile C++
echo "Compiling C++..."
cd "$ALGO_DIR/cpp"
g++ -O3 -std=c++17 hopcroft_karp.cpp -o hopcroft_karp_cpp
if [ $? -ne 0 ]; then
    echo "C++ compilation failed!"
    exit 1
fi
echo "C++ compiled successfully!"
echo ""

# Compile Rust
echo "Compiling Rust..."
cd "$ALGO_DIR/rust"
rustc -O hopcroft_karp.rs -o hopcroft_karp_rust
if [ $? -ne 0 ]; then
    echo "Rust compilation failed!"
    echo "Make sure Rust is installed: https://rustup.rs/"
    exit 1
fi
echo "Rust compiled successfully!"
echo ""

# Run Python
echo "========================================="
echo "Running Python implementation..."
echo "========================================="
cd "$ALGO_DIR/python"
python3 hopcroft_karp.py
echo ""

# Run C++
echo "========================================="
echo "Running C++ implementation..."
echo "========================================="
cd "$ALGO_DIR/cpp"
./hopcroft_karp_cpp
echo ""

# Run Rust
echo "========================================="
echo "Running Rust implementation..."
echo "========================================="
cd "$ALGO_DIR/rust"
./hopcroft_karp_rust
echo ""

echo "========================================="
echo "Benchmark Complete!"
echo "========================================="

# Return to original directory
cd "$SCRIPT_DIR"

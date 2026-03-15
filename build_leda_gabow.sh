#!/bin/bash
# Build LEDA Gabow MCM test harness with greedy-md support
#
# Usage:
#   bash build_leda_gabow.sh /path/to/LEDA-7
#
# Builds in current directory. Produces: ./test_leda_gabow
#
# Then run:
#   ./test_leda_gabow <graph_file> [--greedy|--greedy-md]

if [ $# -lt 1 ]; then
    echo "Usage: bash build_leda_gabow.sh /path/to/LEDA-7"
    exit 1
fi

LEDA_ROOT="$1"
INCL="$LEDA_ROOT/incl"
SRC="$LEDA_ROOT/src"

if [ ! -d "$INCL/LEDA" ]; then
    echo "ERROR: $INCL/LEDA not found."
    exit 1
fi

CXX="${CXX:-g++}"
CXXFLAGS="-O3 -std=c++17 -I $INCL"

# ---- Step 1: Build LEDA library (skip if exists) ----
if [ -f leda_obj/libleda.a ]; then
    echo "=== LEDA library exists, skipping build ==="
else
    echo "=== Building LEDA static library ==="
    mkdir -p leda_obj
    ok=0; fail=0
    for f in $(find "$SRC" -name "*.cpp"); do
        base=$(basename "$f" .cpp)
        if $CXX $CXXFLAGS -c "$f" -o "leda_obj/${base}.o" 2>/dev/null; then
            ok=$((ok+1))
        else
            fail=$((fail+1))
        fi
    done
    echo "Compiled: $ok OK, $fail failed"
    ar rcs leda_obj/libleda.a leda_obj/*.o
    echo "Library: leda_obj/libleda.a"
fi

# ---- Step 2: Patch LEDA headers ----
echo ""
echo "=== Patching LEDA headers ==="

# Patch 1: mc_matching.h circular dependency
MC_H="$INCL/LEDA/graph/mc_matching.h"
if grep -q "^inline list<edge> MAX_CARD_MATCHING_GABOW" "$MC_H"; then
    cp "$MC_H" "${MC_H}.bak"
    python3 -c "
text = open('$MC_H').read()
old = '''inline list<edge> MAX_CARD_MATCHING_GABOW(const graph& G, node_array<int>& OSC)
{ G_card_matching CM(G);
  int num_phases = 0;
  return CM.solve(OSC,num_phases);
}'''
new = '''list<edge> MAX_CARD_MATCHING_GABOW(const graph& G, node_array<int>& OSC);'''
text = text.replace(old, new)
open('$MC_H', 'w').write(text)
"
    echo "Patched mc_matching.h"
else
    echo "mc_matching.h already patched"
fi

# Patch 2: mc_matching_gabow.h — make init() public
GAB_H="$INCL/LEDA/graph/mc_matching_gabow.h"
if grep -q "int init();" "$GAB_H"; then
    # Check if init is still in private section (before public:)
    if python3 -c "
t = open('$GAB_H').read()
ip = t.find('int init();')
pp = t.find('public:', ip) if ip >= 0 else -1
# If init appears before public:, it's private
import sys
sys.exit(0 if (ip >= 0 and pp >= 0 and ip < pp) else 1)
" 2>/dev/null; then
        cp "$GAB_H" "${GAB_H}.bak"
        python3 -c "
text = open('$GAB_H').read()
# Move init declarations from private to public
old = '''  int init();
  int init(const list<edge>& M);

public: 

  G_card_matching(const graph& G);'''
new = '''public: 

  G_card_matching(const graph& G);

  int init();
  int init(const list<edge>& M);'''
text = text.replace(old, new)
open('$GAB_H', 'w').write(text)
"
        echo "Patched mc_matching_gabow.h (init now public)"
    else
        echo "mc_matching_gabow.h init already public"
    fi
else
    echo "mc_matching_gabow.h: init not found (unexpected)"
fi

# ---- Step 3: Create stubs ----
echo ""
echo "=== Creating stubs ==="
cat > leda_stubs.cpp << 'STUBS'
#include <LEDA/numbers/integer.h>
#include <LEDA/numbers/rational.h>
#include <LEDA/numbers/real.h>
#include <LEDA/numbers/vector.h>
#include <LEDA/core/string.h>

int leda_set_fpu_defaults() { return 0; }

namespace leda {
    int integer::cmp(const integer&, const integer&) { return 0; }
    int rational::cmp(const rational&, const rational&) { return 0; }
    int compare(const real&, const real&) { return 0; }
    vector::vector(int d) {}
    bool is_file(string) { return false; }
    bool delete_file(string) { return false; }
    string tmp_file_name() { return string(""); }
}
STUBS
echo "Created leda_stubs.cpp"

# ---- Step 4: Create test harness ----
echo ""
echo "=== Creating test harness ==="
cat > test_leda_gabow.cpp << 'HARNESS'
#define __exportC
#define __exportF
#include <LEDA/graph/graph.h>
#include <LEDA/graph/mc_matching_gabow.h>
#include <cstdio>
#include <chrono>
#include <climits>
#include <vector>
#include <string>
#include <algorithm>

using namespace leda;

int main(int argc, char* argv[]) {
    if (argc < 2) { printf("Usage: %s <filename> [--greedy|--greedy-md]\n", argv[0]); return 1; }

    int greedy_mode = 0;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedy_mode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedy_mode = 2;
    }

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    int n, m;
    if (fscanf(f, "%d %d", &n, &m) != 2) { fprintf(stderr, "Bad header\n"); fclose(f); return 1; }

    graph G;
    std::vector<node> nd(n);
    for (int i = 0; i < n; i++) nd[i] = G.new_node();
    for (int i = 0; i < m; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        if (u >= 0 && u < n && v >= 0 && v < n && u != v)
            G.new_edge(nd[u], nd[v]);
    }
    fclose(f);

    printf("Gabow MCM (LEDA) - C++ Implementation\n");
    printf("======================================\n\n");
    printf("Graph: %d vertices, %d edges\n", n, G.number_of_edges());

    auto t0 = std::chrono::high_resolution_clock::now();
    int greedy_count = 0;

    G_card_matching solver(G);

    if (greedy_mode > 0) {
        node_array<node> mate_ext(G, nil);

        if (greedy_mode == 1) {
            edge e;
            forall_edges(e, G) {
                node u = G.source(e), v = G.target(e);
                if (mate_ext[u] == nil && mate_ext[v] == nil) {
                    mate_ext[u] = v; mate_ext[v] = u;
                    greedy_count++;
                }
            }
        } else {
            node_array<int> deg(G, 0);
            edge e;
            forall_edges(e, G) { deg[G.source(e)]++; deg[G.target(e)]++; }
            std::vector<std::pair<int,node>> order;
            node v;
            forall_nodes(v, G) order.push_back({deg[v], v});
            std::sort(order.begin(), order.end());
            for (auto& [d, u] : order) {
                if (mate_ext[u] != nil) continue;
                int best_deg = INT_MAX;
                node best = nil;
                forall_inout_edges(e, u) {
                    node w = G.opposite(u, e);
                    if (mate_ext[w] == nil && deg[w] < best_deg) {
                        best = w; best_deg = deg[w];
                    }
                }
                if (best != nil) {
                    mate_ext[u] = best; mate_ext[best] = u;
                    greedy_count++;
                }
            }
        }

        /* Pass initial matching to LEDA */
        list<edge> M0;
        edge e;
        forall_edges(e, G) {
            node u = G.source(e), v = G.target(e);
            if (mate_ext[u] == v) {
                M0.append(e);
                mate_ext[u] = nil; mate_ext[v] = nil;
            }
        }
        solver.init(M0);
    }

    node_array<int> OSC(G);
    int num_iters = 0;
    list<edge> M = solver.solve(OSC, num_iters, 1.0, false);
    auto t1 = std::chrono::high_resolution_clock::now();

    printf("Phases: %d\n", num_iters);

    node_array<int> vdeg(G, 0);
    edge e;
    bool ok = true;
    forall(e, M) { vdeg[G.source(e)]++; vdeg[G.target(e)]++; }
    node v;
    forall_nodes(v, G) if (vdeg[v] > 1) { ok = false; break; }

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %d\n", M.size());
    printf("%s\n", ok ? "VALIDATION PASSED" : "VALIDATION FAILED");
    printf("=========================\n\n");
    printf("Matching size: %d\n", M.size());
    if (greedy_mode > 0) {
        printf("Greedy init size: %d\n", greedy_count);
        if (M.size() > 0) printf("Greedy/Final: %.2f%%\n", 100.0 * greedy_count / M.size());
    }
    printf("Time: %ld ms\n",
        (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
    return 0;
}
HARNESS
echo "Created test_leda_gabow.cpp"

# ---- Step 5: Compile ----
echo ""
echo "=== Compiling test_leda_gabow ==="
if $CXX $CXXFLAGS \
    test_leda_gabow.cpp leda_stubs.cpp \
    -L leda_obj -lleda -lpthread \
    -o test_leda_gabow; then
    echo ""
    echo "=== Done ==="
    echo "Binary: ./test_leda_gabow"
    echo "Usage:  ./test_leda_gabow <graph_file> [--greedy|--greedy-md]"
else
    echo "=== FAILED ==="
    exit 1
fi

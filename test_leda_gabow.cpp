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

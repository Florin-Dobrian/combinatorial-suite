/*
 * Gabow's Scaling Algorithm (Optimized) - O(E√V) Maximum Matching
 *
 * Phase 1: Build level structure by distance, detect blossoms.
 * Phase 2: Find and augment all shortest augmenting paths.
 * Repeat until no augmenting paths exist.
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <chrono>

static const int NIL = -1;
static const int UNLABELED = 0;
static const int EVEN = 1;
static const int ODD = 2;

struct GabowOptimized {
    int n;
    std::vector<std::vector<int>> graph;

    std::vector<int> mate;
    std::vector<int> label;          /* EVEN=1, ODD=2, UNLABELED=0 */
    std::vector<int> base;
    std::vector<int> parent;
    std::vector<int> source_bridge;
    std::vector<int> target_bridge;
    std::vector<std::vector<std::pair<int,int>>> edge_queue;

    int Delta;

    GabowOptimized(int n, const std::vector<std::pair<int,int>>& edges) : n(n) {
        graph.resize(n);
        mate.assign(n, NIL);
        label.assign(n, UNLABELED);
        base.resize(n);
        parent.assign(n, NIL);
        source_bridge.assign(n, NIL);
        target_bridge.assign(n, NIL);
        edge_queue.resize(n + 1);

        for (auto& e : edges) {
            int u = e.first, v = e.second;
            if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
                graph[u].push_back(v);
                graph[v].push_back(u);
            }
        }
        for (int i = 0; i < n; i++) std::sort(graph[i].begin(), graph[i].end());
    }

    int find_base(int v) {
        if (base[v] != v) base[v] = find_base(base[v]);
        return base[v];
    }

    int find_lca(int u, int v) {
        std::vector<bool> marked(n, false);
        int x = find_base(u);
        while (x != NIL) {
            marked[x] = true;
            if (mate[x] == NIL || parent[mate[x]] == NIL) break;
            x = find_base(parent[mate[x]]);
        }
        int y = find_base(v);
        while (y != NIL) {
            if (marked[y]) return y;
            if (mate[y] == NIL || parent[mate[y]] == NIL) break;
            y = find_base(parent[mate[y]]);
        }
        return NIL;
    }

    void shrink_path(int lca, int x, int y) {
        int v = find_base(x);
        while (v != lca) {
            base[v] = lca;
            int mv = mate[v];
            if (mv == NIL) break;
            base[mv] = lca;
            source_bridge[mv] = x;
            target_bridge[mv] = y;
            if (parent[mv] == NIL) break;
            v = find_base(parent[mv]);
        }
    }

    void scan_edge(int u, int v) {
        if (Delta < (int)edge_queue.size())
            edge_queue[Delta].push_back({u, v});
    }

    bool phase_1() {
        Delta = 0;
        bool found = false;

        for (auto& q : edge_queue) q.clear();
        for (int i = 0; i < n; i++) {
            base[i] = i;
            label[i] = (mate[i] == NIL) ? EVEN : UNLABELED;
            parent[i] = NIL;
            source_bridge[i] = NIL;
            target_bridge[i] = NIL;
        }

        for (int v = 0; v < n; v++) {
            if (mate[v] == NIL) {
                for (int u : graph[v]) scan_edge(v, u);
            }
        }

        while (Delta <= n && !found) {
            while (!edge_queue[Delta].empty()) {
                auto [x, y] = edge_queue[Delta].back();
                edge_queue[Delta].pop_back();

                int bx = find_base(x), by = find_base(y);
                if (label[bx] != EVEN) { std::swap(x, y); std::swap(bx, by); }
                if (bx == by || label[bx] != EVEN || y == mate[x] || label[by] == ODD) continue;

                if (label[by] == UNLABELED) {
                    int z = mate[y];
                    if (z != NIL) {
                        label[y] = ODD;
                        label[z] = EVEN;
                        parent[y] = x;
                        parent[z] = y;
                        for (int w : graph[z]) scan_edge(z, w);
                    }
                } else if (label[by] == EVEN) {
                    int lca = find_lca(x, y);
                    if (lca != NIL) { shrink_path(lca, x, y); shrink_path(lca, y, x); }
                    else { found = true; break; }
                }
            }
            if (!found) Delta++;
        }
        return found;
    }

    void phase_2() {
        for (int start = 0; start < n; start++) {
            if (mate[start] != NIL || label[start] != EVEN) continue;

            std::vector<int> queue(n);
            std::vector<int> pred(n, NIL);
            std::vector<bool> vis(n, false);
            int qh = 0, qt = 0;

            queue[qt++] = start;
            vis[find_base(start)] = true;
            int endpoint = NIL;

            while (qh < qt && endpoint == NIL) {
                int u = queue[qh++];
                for (int v : graph[u]) {
                    int bu = find_base(u), bv = find_base(v);
                    if (bu == bv || vis[bv]) continue;
                    if (mate[v] == NIL && v != start) { pred[v] = u; endpoint = v; break; }
                    if (label[bv] != ODD) {
                        pred[v] = u;
                        vis[bv] = true;
                        int mv = mate[v];
                        if (mv != NIL && !vis[find_base(mv)]) {
                            pred[mv] = v;
                            vis[find_base(mv)] = true;
                            queue[qt++] = mv;
                        }
                    }
                }
            }

            if (endpoint != NIL) {
                std::vector<int> path;
                for (int c = endpoint; c != NIL; c = pred[c]) path.push_back(c);
                std::reverse(path.begin(), path.end());
                for (size_t i = 0; i + 1 < path.size(); i += 2) {
                    mate[path[i]] = path[i + 1];
                    mate[path[i + 1]] = path[i];
                }
            }
        }
    }

    std::vector<std::pair<int,int>> maximum_matching() {
        while (phase_1()) phase_2();

        std::vector<std::pair<int,int>> matching;
        for (int u = 0; u < n; u++) {
            if (mate[u] != NIL && mate[u] > u)
                matching.push_back({u, mate[u]});
        }
        std::sort(matching.begin(), matching.end());
        return matching;
    }
};

void validate_matching(int n, const std::vector<std::vector<int>>& graph,
                       const std::vector<std::pair<int,int>>& matching) {
    std::vector<int> deg(n, 0);
    int errors = 0;

    for (auto& e : matching) {
        if (!std::binary_search(graph[e.first].begin(), graph[e.first].end(), e.second)) {
            fprintf(stderr, "ERROR: Edge (%d, %d) not in graph!\n", e.first, e.second);
            errors++;
        }
        deg[e.first]++;
        deg[e.second]++;
    }
    for (int i = 0; i < n; i++) {
        if (deg[i] > 1) { fprintf(stderr, "ERROR: Vertex %d in %d edges!\n", i, deg[i]); errors++; }
    }
    int matched = 0;
    for (int i = 0; i < n; i++) if (deg[i] > 0) matched++;

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %d\n", (int)matching.size());
    printf("Matched vertices: %d\n", matched);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

int main(int argc, char* argv[]) {
    printf("Gabow's Scaling Algorithm (Optimized) - C++ Implementation\n");
    printf("============================================================\n\n");

    if (argc < 2) { printf("Usage: %s <filename>\n", argv[0]); return 1; }

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open file: %s\n", argv[1]); return 1; }

    int n, m;
    if (fscanf(f, "%d %d", &n, &m) != 2) { fprintf(stderr, "Bad header\n"); fclose(f); return 1; }

    std::vector<std::pair<int,int>> edges;
    edges.reserve(m);
    for (int i = 0; i < m; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edges.push_back({u, v});
    }
    fclose(f);

    printf("Graph: %d vertices, %d edges\n", n, (int)edges.size());

    auto t0 = std::chrono::high_resolution_clock::now();
    GabowOptimized gabow(n, edges);
    auto matching = gabow.maximum_matching();
    auto t1 = std::chrono::high_resolution_clock::now();

    validate_matching(n, gabow.graph, matching);

    printf("Matching size: %d\n", (int)matching.size());
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

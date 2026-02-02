/*
 * Gabow's Algorithm (Simple) - O(V * E) Maximum Matching
 *
 * BFS-based augmenting path search with blossom contraction using
 * path-compressed union-find. One augmenting path per iteration.
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <chrono>

static const int NIL = -1;

struct GabowSimple {
    int n;
    std::vector<std::vector<int>> graph;
    std::vector<int> mate;
    std::vector<int> base;
    std::vector<int> parent;
    std::vector<bool> blossom;
    std::vector<bool> visited;

    GabowSimple(int n, const std::vector<std::pair<int,int>>& edges) : n(n) {
        graph.resize(n);
        mate.assign(n, NIL);
        base.resize(n);
        parent.resize(n);
        blossom.resize(n);
        visited.resize(n);

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
            if (mate[x] == NIL) break;
            if (parent[mate[x]] == NIL) break;
            x = find_base(parent[mate[x]]);
        }
        int y = find_base(v);
        while (y != NIL) {
            if (marked[y]) return y;
            if (mate[y] == NIL) break;
            if (parent[mate[y]] == NIL) break;
            y = find_base(parent[mate[y]]);
        }
        return NIL;
    }

    void mark_blossom(int u, int lca, std::vector<int>& queue, int& qtail) {
        while (find_base(u) != lca) {
            int bv = find_base(u);
            int bw = find_base(mate[u]);
            blossom[bv] = blossom[bw] = true;
            if (!visited[bw]) { visited[bw] = true; queue[qtail++] = bw; }
            if (mate[u] == NIL || parent[mate[u]] == NIL) break;
            u = parent[mate[u]];
        }
    }

    void contract_blossom(int u, int v, std::vector<int>& queue, int& qtail) {
        int lca = find_lca(u, v);
        blossom.assign(n, false);
        mark_blossom(u, lca, queue, qtail);
        mark_blossom(v, lca, queue, qtail);
        for (int i = 0; i < n; i++) {
            if (blossom[find_base(i)]) {
                base[i] = lca;
                if (!visited[i]) { visited[i] = true; queue[qtail++] = i; }
            }
        }
    }

    bool find_augmenting_path(int start) {
        for (int i = 0; i < n; i++) { base[i] = i; parent[i] = NIL; }
        visited.assign(n, false);

        std::vector<int> queue(n);
        int qhead = 0, qtail = 0;
        queue[qtail++] = start;
        visited[start] = true;

        while (qhead < qtail) {
            int u = queue[qhead++];
            for (int v : graph[u]) {
                int bu = find_base(u), bv = find_base(v);
                if (bu == bv) continue;

                if (mate[v] == NIL) {
                    parent[v] = u;
                    return true;
                }
                if (!visited[bv]) {
                    parent[v] = u;
                    visited[bv] = true;
                    int w = mate[v];
                    visited[find_base(w)] = true;
                    queue[qtail++] = w;
                } else {
                    /* potential blossom — check same tree */
                    int ru = bu, rv = bv;
                    while (mate[ru] != NIL && parent[mate[ru]] != NIL) ru = find_base(parent[mate[ru]]);
                    while (mate[rv] != NIL && parent[mate[rv]] != NIL) rv = find_base(parent[mate[rv]]);
                    if (ru == rv) contract_blossom(u, v, queue, qtail);
                }
            }
        }
        return false;
    }

    void augment_path(int v) {
        while (v != NIL) {
            int pv = parent[v];
            int ppv = mate[pv];
            mate[v] = pv;
            mate[pv] = v;
            v = ppv;
        }
    }

    std::vector<std::pair<int,int>> maximum_matching() {
        bool found = true;
        while (found) {
            found = false;
            for (int v = 0; v < n; v++) {
                if (mate[v] == NIL) {
                    if (find_augmenting_path(v)) {
                        for (int u = 0; u < n; u++) {
                            if (mate[u] == NIL && parent[u] != NIL) {
                                augment_path(u);
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
        }

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
    printf("Gabow's Algorithm (Simple) - C++ Implementation\n");
    printf("=================================================\n\n");

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
    GabowSimple gabow(n, edges);
    auto matching = gabow.maximum_matching();
    auto t1 = std::chrono::high_resolution_clock::now();

    validate_matching(n, gabow.graph, matching);

    printf("Matching size: %d\n", (int)matching.size());
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

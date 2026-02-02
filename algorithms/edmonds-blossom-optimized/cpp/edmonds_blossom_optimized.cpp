/*
 * Edmonds' Blossom Algorithm (Optimized) - O(V^2 * E) Maximum Matching
 *
 * Optimized implementation with BFS augmenting path search and blossom contraction.
 * One augmenting path per iteration.
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <chrono>

static const int NIL = -1;

struct EdmondsBlossomOptimized {
    int n;
    std::vector<std::vector<int>> graph;
    std::vector<int> mate;

    EdmondsBlossomOptimized(int n, const std::vector<std::pair<int,int>>& edges) : n(n) {
        graph.resize(n);
        mate.assign(n, NIL);
        for (auto& e : edges) {
            int u = e.first, v = e.second;
            if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
                graph[u].push_back(v);
                graph[v].push_back(u);
            }
        }
        for (int i = 0; i < n; i++) std::sort(graph[i].begin(), graph[i].end());
    }

    int find_base(int v, std::vector<int>& base) {
        while (base[v] != v) v = base[v] = base[base[v]];
        return v;
    }

    int find_blossom_base(int v, int w, std::vector<int>& parent, std::vector<int>& base) {
        std::vector<bool> on_path(n, false);
        int cur = v;
        while (cur != NIL) {
            on_path[find_base(cur, base)] = true;
            cur = parent[cur];
        }
        cur = w;
        while (cur != NIL) {
            int b = find_base(cur, base);
            if (on_path[b]) return b;
            cur = parent[cur];
        }
        return find_base(v, base);
    }

    void trace_and_update(int start, int blossom_base, std::vector<int>& base,
                          std::vector<int>& label, std::vector<int>& parent,
                          std::vector<int>& queue, int& qtail) {
        int cur = start;
        while (true) {
            int cb = find_base(cur, base);
            if (cb == blossom_base) break;
            base[cur] = blossom_base;
            if (label[cur] == 2) {
                label[cur] = 1;
                queue[qtail++] = cur;
            }
            if (mate[cur] != NIL) {
                base[mate[cur]] = blossom_base;
                if (parent[mate[cur]] != NIL) cur = parent[mate[cur]];
                else break;
            } else break;
        }
    }

    std::vector<int> find_augmenting_path(int start) {
        std::vector<int> parent(n, NIL);
        std::vector<int> base(n);
        for (int i = 0; i < n; i++) base[i] = i;
        std::vector<int> label(n, 0);

        std::vector<int> queue(n);
        int qhead = 0, qtail = 0;

        label[start] = 1;
        queue[qtail++] = start;

        while (qhead < qtail) {
            int v = queue[qhead++];
            int vb = find_base(v, base);

            for (int w : graph[v]) {
                int wb = find_base(w, base);
                if (vb == wb) continue;

                if (label[w] == 0) {
                    if (mate[w] != NIL) {
                        label[w] = 2;
                        label[mate[w]] = 1;
                        parent[w] = v;
                        parent[mate[w]] = w;
                        queue[qtail++] = mate[w];
                    } else {
                        std::vector<int> path = {w, v};
                        int cur = v;
                        while (parent[cur] != NIL) { path.push_back(parent[cur]); cur = parent[cur]; }
                        return path;
                    }
                }
                else if (label[w] == 1) {
                    int bb = find_blossom_base(v, w, parent, base);
                    trace_and_update(v, bb, base, label, parent, queue, qtail);
                    trace_and_update(w, bb, base, label, parent, queue, qtail);
                }
            }
        }
        return {};
    }

    void augment(const std::vector<int>& path) {
        for (size_t i = 0; i + 1 < path.size(); i += 2) {
            mate[path[i]] = path[i + 1];
            mate[path[i + 1]] = path[i];
        }
    }

    std::vector<std::pair<int,int>> maximum_matching() {
        bool improved = true;
        while (improved) {
            improved = false;
            for (int v = 0; v < n; v++) {
                if (mate[v] == NIL) {
                    auto path = find_augmenting_path(v);
                    if (!path.empty()) { augment(path); improved = true; break; }
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
    printf("Edmonds' Blossom Algorithm (Optimized) - C++ Implementation\n");
    printf("=============================================================\n\n");

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
    EdmondsBlossomOptimized eb(n, edges);
    auto matching = eb.maximum_matching();
    auto t1 = std::chrono::high_resolution_clock::now();

    validate_matching(n, eb.graph, matching);

    printf("Matching size: %d\n", (int)matching.size());
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

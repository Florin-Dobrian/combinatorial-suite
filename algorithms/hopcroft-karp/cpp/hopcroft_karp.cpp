/*
 * Hopcroft-Karp Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * BFS to find shortest augmenting path length, then DFS to find
 * all vertex-disjoint augmenting paths of that length.
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <chrono>
#include <climits>

static const int NIL = -1;

struct HopcroftKarp {
    int left_count;
    int right_count;
    std::vector<std::vector<int>> graph; /* graph[u] = list of right nodes */
    std::vector<int> pair_left;
    std::vector<int> pair_right;
    std::vector<int> dist;

    HopcroftKarp(int lc, int rc, const std::vector<std::pair<int,int>>& edges)
        : left_count(lc), right_count(rc) {
        graph.resize(lc);
        for (auto& e : edges) {
            int u = e.first, v = e.second;
            if (u >= 0 && u < lc && v >= 0 && v < rc)
                graph[u].push_back(v);
        }
        for (int i = 0; i < lc; i++) std::sort(graph[i].begin(), graph[i].end());

        pair_left.assign(lc, NIL);
        pair_right.assign(rc, NIL);
        dist.resize(lc + 1);
    }

    bool bfs() {
        std::vector<int> queue(left_count);
        int qh = 0, qt = 0;

        for (int u = 0; u < left_count; u++) {
            if (pair_left[u] == NIL) { dist[u] = 0; queue[qt++] = u; }
            else dist[u] = INT_MAX;
        }
        dist[left_count] = INT_MAX; /* NIL sentinel */

        while (qh < qt) {
            int u = queue[qh++];
            if (dist[u] < dist[left_count]) {
                for (int v : graph[u]) {
                    int pn = (pair_right[v] == NIL) ? left_count : pair_right[v];
                    if (dist[pn] == INT_MAX) {
                        dist[pn] = dist[u] + 1;
                        if (pair_right[v] != NIL) queue[qt++] = pair_right[v];
                    }
                }
            }
        }
        return dist[left_count] != INT_MAX;
    }

    bool dfs(int u) {
        if (u == NIL) return true;
        for (int v : graph[u]) {
            int pn = (pair_right[v] == NIL) ? left_count : pair_right[v];
            if (dist[pn] == dist[u] + 1) {
                if (dfs(pair_right[v])) {
                    pair_right[v] = u;
                    pair_left[u] = v;
                    return true;
                }
            }
        }
        dist[u] = INT_MAX;
        return false;
    }

    std::vector<std::pair<int,int>> maximum_matching() {
        while (bfs()) {
            for (int u = 0; u < left_count; u++) {
                if (pair_left[u] == NIL) dfs(u);
            }
        }

        std::vector<std::pair<int,int>> matching;
        for (int u = 0; u < left_count; u++) {
            if (pair_left[u] != NIL) matching.push_back({u, pair_left[u]});
        }
        std::sort(matching.begin(), matching.end());
        return matching;
    }
};

void validate_matching(int lc, int rc, const std::vector<std::vector<int>>& graph,
                       const std::vector<std::pair<int,int>>& matching) {
    std::vector<int> ldeg(lc, 0), rdeg(rc, 0);
    int errors = 0;

    for (auto& e : matching) {
        if (!std::binary_search(graph[e.first].begin(), graph[e.first].end(), e.second)) {
            fprintf(stderr, "ERROR: Edge (%d, %d) not in graph!\n", e.first, e.second);
            errors++;
        }
        ldeg[e.first]++;
        rdeg[e.second]++;
    }
    for (int i = 0; i < lc; i++) {
        if (ldeg[i] > 1) { fprintf(stderr, "ERROR: Left %d in %d edges!\n", i, ldeg[i]); errors++; }
    }
    for (int i = 0; i < rc; i++) {
        if (rdeg[i] > 1) { fprintf(stderr, "ERROR: Right %d in %d edges!\n", i, rdeg[i]); errors++; }
    }
    int matched_l = 0, matched_r = 0;
    for (int i = 0; i < lc; i++) if (ldeg[i] > 0) matched_l++;
    for (int i = 0; i < rc; i++) if (rdeg[i] > 0) matched_r++;

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %d\n", (int)matching.size());
    printf("Left matched: %d, Right matched: %d\n", matched_l, matched_r);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

int main(int argc, char* argv[]) {
    printf("Hopcroft-Karp Algorithm - C++ Implementation\n");
    printf("==============================================\n\n");

    if (argc < 2) { printf("Usage: %s <filename>\n", argv[0]); return 1; }

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open file: %s\n", argv[1]); return 1; }

    int lc, rc, m;
    if (fscanf(f, "%d %d %d", &lc, &rc, &m) != 3) { fprintf(stderr, "Bad header\n"); fclose(f); return 1; }

    std::vector<std::pair<int,int>> edges;
    edges.reserve(m);
    for (int i = 0; i < m; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edges.push_back({u, v});
    }
    fclose(f);

    printf("Graph: %d left, %d right, %d edges\n", lc, rc, (int)edges.size());

    auto t0 = std::chrono::high_resolution_clock::now();
    HopcroftKarp hk(lc, rc, edges);
    auto matching = hk.maximum_matching();
    auto t1 = std::chrono::high_resolution_clock::now();

    validate_matching(lc, rc, hk.graph, matching);

    printf("Matching size: %d\n", (int)matching.size());
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

/*
 * Hopcroft-Karp Hybrid Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * CSR variant: adjacency stored as flat (adj_off, adj_edges) arrays
 * instead of vector<vector<int>>. Algorithm unchanged.
 *
 * Old HK's lean BFS (single dist[] array, sentinel trick, no status enums)
 * + iterative stack-based DFS with edge index array (no recursion, no rescan).
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <climits>

static const int NIL = -1;

struct HopcroftKarpHybrid {
    int left_count;
    int greedy_size = 0;
    int right_count;
    /* CSR adjacency: adj_edges[adj_off[u] .. adj_off[u+1]) are the right
     * nodes neighboring left node u.  Replaces vector<vector<int>>. */
    std::vector<int> adj_off;
    std::vector<int> adj_edges;
    std::vector<int> pair_left;
    std::vector<int> pair_right;
    std::vector<int> dist;

    /* Edge index: per left vertex, offset within its adjacency range (persistent within a phase) */
    std::vector<int> edge_idx;

    /* DFS stack: see original source for protocol. */
    std::vector<int> stk_u;
    std::vector<int> stk_v;
    int stk_top;

    HopcroftKarpHybrid(int lc, int rc, const std::vector<std::pair<int,int>>& edges)
        : left_count(lc), right_count(rc), stk_top(0) {
        std::vector<std::vector<int>> tmp(lc);
        for (auto& e : edges) {
            int u = e.first, v = e.second;
            if (u >= 0 && u < lc && v >= 0 && v < rc)
                tmp[u].push_back(v);
        }
        for (int i = 0; i < lc; i++) {
            std::sort(tmp[i].begin(), tmp[i].end());
            tmp[i].erase(std::unique(tmp[i].begin(), tmp[i].end()), tmp[i].end());
        }

        adj_off.assign(lc + 1, 0);
        for (int i = 0; i < lc; i++) adj_off[i + 1] = adj_off[i] + (int)tmp[i].size();
        adj_edges.resize(adj_off[lc]);
        for (int i = 0; i < lc; i++)
            std::copy(tmp[i].begin(), tmp[i].end(), adj_edges.begin() + adj_off[i]);

        pair_left.assign(lc, NIL);
        pair_right.assign(rc, NIL);
        dist.resize(lc + 1);
        edge_idx.assign(lc, 0);
        stk_u.resize(lc);
        stk_v.resize(lc);
    }

    /* BFS: identical to original HK, adapted for CSR adjacency */
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
                int s = adj_off[u], e = adj_off[u + 1];
                for (int j = s; j < e; j++) {
                    int v = adj_edges[j];
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

    /*
     * DFS: iterative with edge index.
     *
     * edge_idx[u] is an offset WITHIN u's adjacency range [adj_off[u], adj_off[u+1]).
     * So the "current candidate edge" is adj_edges[adj_off[u] + edge_idx[u]].
     */
    bool dfs(int root) {
        stk_top = 0;
        stk_u[0] = root;
        stk_top = 1;

        while (stk_top > 0) {
            int u = stk_u[stk_top - 1];
            int s = adj_off[u], e = adj_off[u + 1];
            int sz = e - s;

            bool pushed = false;
            while (edge_idx[u] < sz) {
                int v = adj_edges[s + edge_idx[u]];
                int pn = (pair_right[v] == NIL) ? left_count : pair_right[v];
                if (dist[pn] != dist[u] + 1) {
                    edge_idx[u]++;
                    continue;
                }

                stk_v[stk_top - 1] = v;
                edge_idx[u]++;

                if (pair_right[v] == NIL) {
                    /* Found augmenting path — augment all the way back */
                    for (int d = stk_top - 1; d >= 0; d--) {
                        pair_right[stk_v[d]] = stk_u[d];
                        pair_left[stk_u[d]] = stk_v[d];
                    }
                    return true;
                }

                stk_u[stk_top] = pair_right[v];
                stk_top++;
                pushed = true;
                break;
            }

            if (!pushed) {
                dist[u] = INT_MAX;
                stk_top--;
            }
        }
        return false;
    }

    /* Greedy: simple sequential */
    int greedy_init() {
        int cnt = 0;
        for (int u = 0; u < left_count; u++) {
            if (pair_left[u] != NIL) continue;
            int s = adj_off[u], e = adj_off[u + 1];
            for (int j = s; j < e; j++) {
                int v = adj_edges[j];
                if (pair_right[v] == NIL) { pair_left[u] = v; pair_right[v] = u; cnt++; break; }
            }
        }
        return cnt;
    }

    /* Greedy: min-degree */
    int greedy_init_md() {
        int cnt = 0;
        std::vector<int> deg(right_count, 0);
        for (int u = 0; u < left_count; u++) {
            int s = adj_off[u], e = adj_off[u + 1];
            for (int j = s; j < e; j++) deg[adj_edges[j]]++;
        }
        std::vector<int> order(left_count);
        for (int i = 0; i < left_count; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b){
            int da = adj_off[a + 1] - adj_off[a];
            int db = adj_off[b + 1] - adj_off[b];
            return da < db || (da == db && a < b);
        });
        for (int u : order) {
            if (pair_left[u] != NIL) continue;
            int best = -1, best_deg = INT_MAX;
            int s = adj_off[u], e = adj_off[u + 1];
            for (int j = s; j < e; j++) {
                int v = adj_edges[j];
                if (pair_right[v] == NIL && deg[v] < best_deg) {
                    best = v; best_deg = deg[v];
                }
            }
            if (best >= 0) { pair_left[u] = best; pair_right[best] = u; cnt++; }
        }
        return cnt;
    }

    std::vector<std::pair<int,int>> maximum_matching(int greedy_mode = 0) {
        int greedy_count = 0;
        if (greedy_mode == 1) greedy_count = greedy_init();
        else if (greedy_mode == 2) greedy_count = greedy_init_md();
        greedy_size = greedy_count;
        int phases = 0;
        while (bfs()) {
            phases++;
            for (int u = 0; u < left_count; u++) edge_idx[u] = 0;
            for (int u = 0; u < left_count; u++) {
                if (pair_left[u] == NIL) dfs(u);
            }
        }

        printf("Phases: %d\n", phases);

        std::vector<std::pair<int,int>> matching;
        for (int u = 0; u < left_count; u++) {
            if (pair_left[u] != NIL) matching.push_back({u, pair_left[u]});
        }
        std::sort(matching.begin(), matching.end());
        return matching;
    }
};

void validate_matching(int lc, int rc,
                       const std::vector<int>& adj_off,
                       const std::vector<int>& adj_edges,
                       const std::vector<std::pair<int,int>>& matching) {
    std::vector<int> ldeg(lc, 0), rdeg(rc, 0);
    int errors = 0;

    for (auto& e : matching) {
        int u = e.first, v = e.second;
        int s = adj_off[u], en = adj_off[u + 1];
        if (!std::binary_search(adj_edges.begin() + s, adj_edges.begin() + en, v)) {
            fprintf(stderr, "ERROR: Edge (%d, %d) not in graph!\n", u, v);
            errors++;
        }
        ldeg[u]++;
        rdeg[v]++;
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
    printf("Hopcroft-Karp Hybrid CSR Algorithm - C++ Implementation\n");
    printf("==========================================================\n\n");

    if (argc < 2) { printf("Usage: %s <filename> [--greedy|--greedy-md]\n", argv[0]); return 1; }
    int greedy_mode = 0;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedy_mode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedy_mode = 2;
    }

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
    HopcroftKarpHybrid hk(lc, rc, edges);
    auto matching = hk.maximum_matching(greedy_mode);
    auto t1 = std::chrono::high_resolution_clock::now();

    validate_matching(lc, rc, hk.adj_off, hk.adj_edges, matching);

    printf("Matching size: %d\n", (int)matching.size());

    if (greedy_mode > 0) {
        int gs = hk.greedy_size;
        int fs = (int)matching.size();
        printf("Greedy init size: %d\n", gs);
        if (fs > 0) printf("Greedy/Final: %.2f%%\n", 100.0 * gs / fs);
        else printf("Greedy/Final: NA\n");
    }
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

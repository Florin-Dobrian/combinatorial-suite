/*
 * Hopcroft-Karp Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * VV (vector-of-vectors) adjacency.
 *
 * Refactored to separate input (BipartiteGraph), output (BipartiteMatching),
 * and algorithm state (HKState).
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <climits>

constexpr int32_t NIL = -1;

/* ---------- Input: BipartiteGraph ---------- */

struct BipartiteGraph {
    size_t sNumVtxs;
    size_t tNumVtxs;
    size_t numEdgs;
    std::vector<std::vector<int32_t>> sAdj;  // sAdj[u] = t-neighbors of s-vertex u
    std::vector<std::vector<int32_t>> tAdj;  // tAdj[v] = s-neighbors of t-vertex v
};

/* ---------- Output: BipartiteMatching ---------- */

struct BipartiteMatching {
    size_t sNumVtxs;
    size_t tNumVtxs;
    size_t numEdgs;
    std::vector<int32_t> sMate;   // sMate[u] = t-partner of s-vertex u, or NIL
    std::vector<int32_t> tMate;   // tMate[v] = s-partner of t-vertex v, or NIL
};

/* ---------- State: HKState ---------- */

struct HKState {
    std::vector<int32_t> dist;    // length sNumVtxs+1; dist[sNumVtxs] is the NIL sentinel
};

/* ---------- BipartiteGraph construction ---------- */

BipartiteGraph buildBipartiteGraph(size_t sNumVtxs, size_t tNumVtxs,
                                   const std::vector<std::pair<int32_t,int32_t>>& edges) {
    BipartiteGraph g;
    g.sNumVtxs = sNumVtxs;
    g.tNumVtxs = tNumVtxs;
    g.sAdj.resize(sNumVtxs);
    g.tAdj.resize(tNumVtxs);
    for (auto& e : edges) {
        int32_t u = e.first, v = e.second;
        if (u >= 0 && (size_t)u < sNumVtxs && v >= 0 && (size_t)v < tNumVtxs) {
            g.sAdj[u].push_back(v);
            g.tAdj[v].push_back(u);
        }
    }
    for (size_t u = 0; u < sNumVtxs; u++) {
        std::sort(g.sAdj[u].begin(), g.sAdj[u].end());
        g.sAdj[u].erase(std::unique(g.sAdj[u].begin(), g.sAdj[u].end()), g.sAdj[u].end());
    }
    for (size_t v = 0; v < tNumVtxs; v++) {
        std::sort(g.tAdj[v].begin(), g.tAdj[v].end());
        g.tAdj[v].erase(std::unique(g.tAdj[v].begin(), g.tAdj[v].end()), g.tAdj[v].end());
    }
    size_t total = 0;
    for (size_t u = 0; u < sNumVtxs; u++) total += g.sAdj[u].size();
    g.numEdgs = total;
    return g;
}

/* ---------- BipartiteMatching construction ---------- */

BipartiteMatching emptyBipartiteMatching(const BipartiteGraph& g) {
    BipartiteMatching m;
    m.sNumVtxs = g.sNumVtxs;
    m.tNumVtxs = g.tNumVtxs;
    m.numEdgs = 0;
    m.sMate.assign(g.sNumVtxs, NIL);
    m.tMate.assign(g.tNumVtxs, NIL);
    return m;
}

/* ---------- Greedy initial matching: simple ---------- */

int32_t greedyInit(const BipartiteGraph& g, BipartiteMatching& m) {
    int32_t cnt = 0;
    for (size_t u = 0; u < g.sNumVtxs; u++) {
        if (m.sMate[u] != NIL) continue;
        for (int32_t v : g.sAdj[u]) {
            if (m.tMate[v] == NIL) {
                m.sMate[u] = v;
                m.tMate[v] = (int32_t)u;
                cnt++;
                break;
            }
        }
    }
    m.numEdgs += cnt;
    return cnt;
}

/* ---------- Greedy initial matching: min-degree ---------- */

int32_t greedyInitMd(const BipartiteGraph& g, BipartiteMatching& m) {
    int32_t cnt = 0;
    std::vector<int32_t> deg(g.tNumVtxs, 0);
    for (size_t u = 0; u < g.sNumVtxs; u++)
        for (int32_t v : g.sAdj[u])
            deg[v]++;

    std::vector<int32_t> order(g.sNumVtxs);
    for (size_t i = 0; i < g.sNumVtxs; i++) order[i] = (int32_t)i;
    std::sort(order.begin(), order.end(), [&](int32_t a, int32_t b){
        return g.sAdj[a].size() < g.sAdj[b].size() ||
               (g.sAdj[a].size() == g.sAdj[b].size() && a < b);
    });

    for (int32_t u : order) {
        if (m.sMate[u] != NIL) continue;
        int32_t best = NIL, bestDeg = INT_MAX;
        for (int32_t v : g.sAdj[u]) {
            if (m.tMate[v] == NIL && deg[v] < bestDeg) {
                best = v;
                bestDeg = deg[v];
            }
        }
        if (best != NIL) {
            m.sMate[u] = best;
            m.tMate[best] = u;
            cnt++;
        }
    }
    m.numEdgs += cnt;
    return cnt;
}

/* ---------- HK BFS ---------- */

static bool bfs(const BipartiteGraph& g, const BipartiteMatching& m, HKState& s) {
    std::vector<int32_t> queue(g.sNumVtxs);
    int32_t qh = 0, qt = 0;

    for (size_t u = 0; u < g.sNumVtxs; u++) {
        if (m.sMate[u] == NIL) { s.dist[u] = 0; queue[qt++] = (int32_t)u; }
        else s.dist[u] = INT_MAX;
    }
    s.dist[g.sNumVtxs] = INT_MAX;  /* NIL sentinel */

    while (qh < qt) {
        int32_t u = queue[qh++];
        if (s.dist[u] < s.dist[g.sNumVtxs]) {
            for (int32_t v : g.sAdj[u]) {
                int32_t pn = (m.tMate[v] == NIL) ? (int32_t)g.sNumVtxs : m.tMate[v];
                if (s.dist[pn] == INT_MAX) {
                    s.dist[pn] = s.dist[u] + 1;
                    if (m.tMate[v] != NIL) queue[qt++] = m.tMate[v];
                }
            }
        }
    }
    return s.dist[g.sNumVtxs] != INT_MAX;
}

/* ---------- HK DFS ---------- */

static bool dfs(int32_t u, const BipartiteGraph& g, BipartiteMatching& m, HKState& s) {
    if (u == NIL) return true;
    for (int32_t v : g.sAdj[u]) {
        int32_t pn = (m.tMate[v] == NIL) ? (int32_t)g.sNumVtxs : m.tMate[v];
        if (s.dist[pn] == s.dist[u] + 1) {
            if (dfs(m.tMate[v], g, m, s)) {
                m.tMate[v] = u;
                m.sMate[u] = v;
                return true;
            }
        }
    }
    s.dist[u] = INT_MAX;
    return false;
}

/* ---------- Top-level Hopcroft-Karp ---------- */

void hopcroftKarp(const BipartiteGraph& g, BipartiteMatching& m) {
    HKState s;
    s.dist.assign(g.sNumVtxs + 1, 0);

    int32_t phases = 0;
    int32_t newEdgs = 0;
    while (bfs(g, m, s)) {
        phases++;
        for (size_t u = 0; u < g.sNumVtxs; u++) {
            if (m.sMate[u] == NIL && dfs((int32_t)u, g, m, s)) {
                newEdgs++;
            }
        }
    }
    m.numEdgs += newEdgs;
    printf("Phases: %d\n", phases);
}

/* ---------- Validation ---------- */

void validateBipartiteMatching(const BipartiteGraph& g, const BipartiteMatching& m) {
    int32_t errors = 0;
    int32_t matchedS = 0, matchedT = 0;

    for (size_t u = 0; u < g.sNumVtxs; u++) {
        if (m.sMate[u] != NIL) {
            matchedS++;
            int32_t v = m.sMate[u];
            if (v < 0 || (size_t)v >= g.tNumVtxs) {
                fprintf(stderr, "ERROR: sMate[%zu] = %d out of range\n", u, v);
                errors++;
            } else if (m.tMate[v] != (int32_t)u) {
                fprintf(stderr, "ERROR: sMate[%zu]=%d but tMate[%d]=%d\n", u, v, v, m.tMate[v]);
                errors++;
            } else if (!std::binary_search(g.sAdj[u].begin(), g.sAdj[u].end(), v)) {
                fprintf(stderr, "ERROR: edge (%zu,%d) not in graph\n", u, v);
                errors++;
            }
        }
    }
    for (size_t v = 0; v < g.tNumVtxs; v++) {
        if (m.tMate[v] != NIL) matchedT++;
    }

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %zu\n", m.numEdgs);
    printf("S matched: %d, T matched: %d\n", matchedS, matchedT);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

/* ---------- Main ---------- */

int main(int argc, char* argv[]) {
    printf("Hopcroft-Karp Algorithm - C++ Implementation (VV)\n");
    printf("===================================================\n\n");

    if (argc < 2) { printf("Usage: %s <filename> [--greedy|--greedy-md]\n", argv[0]); return 1; }
    int greedyMode = 0;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedyMode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedyMode = 2;
    }

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open file: %s\n", argv[1]); return 1; }

    int sNumVtxs, tNumVtxs, numEdgsIn;
    if (fscanf(f, "%d %d %d", &sNumVtxs, &tNumVtxs, &numEdgsIn) != 3) {
        fprintf(stderr, "Bad header\n"); fclose(f); return 1;
    }

    std::vector<std::pair<int32_t,int32_t>> edges;
    edges.reserve(numEdgsIn);
    for (int i = 0; i < numEdgsIn; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edges.push_back({(int32_t)u, (int32_t)v});
    }
    fclose(f);

    printf("Graph: %d s-vertices, %d t-vertices, %zu edges\n",
           sNumVtxs, tNumVtxs, edges.size());

    BipartiteGraph bipartiteGraph = buildBipartiteGraph(sNumVtxs, tNumVtxs, edges);
    BipartiteMatching bipartiteMatching = emptyBipartiteMatching(bipartiteGraph);

    auto t0 = std::chrono::high_resolution_clock::now();

    int32_t greedySize = 0;
    if (greedyMode == 1) greedySize = greedyInit(bipartiteGraph, bipartiteMatching);
    else if (greedyMode == 2) greedySize = greedyInitMd(bipartiteGraph, bipartiteMatching);

    hopcroftKarp(bipartiteGraph, bipartiteMatching);

    auto t1 = std::chrono::high_resolution_clock::now();

    validateBipartiteMatching(bipartiteGraph, bipartiteMatching);

    printf("Matching size: %zu\n", bipartiteMatching.numEdgs);
    if (greedyMode > 0) {
        printf("Greedy init size: %d\n", greedySize);
        if (bipartiteMatching.numEdgs > 0)
            printf("Greedy/Final: %.2f%%\n", 100.0 * greedySize / bipartiteMatching.numEdgs);
        else
            printf("Greedy/Final: NA\n");
    }
    printf("Time: %ld ms\n",
           (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
    return 0;
}

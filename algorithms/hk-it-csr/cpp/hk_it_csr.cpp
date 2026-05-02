/*
 * Hopcroft-Karp Iterative Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * CSR adjacency: contiguous flat arrays.
 * Old HK's lean BFS (single sLevel[] array, sentinel trick) + iterative
 * stack-based DFS with edge index array (no recursion, no rescan).
 *
 * Refactored to separate input (BipartiteGraph), output (BipartiteMatching),
 * and algorithm state (HKIState).
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
    std::vector<size_t>  sIdx;
    std::vector<size_t>  tIdx;
    std::vector<int32_t> sAdj;
    std::vector<int32_t> tAdj;
};

/* ---------- Output: BipartiteMatching ---------- */

struct BipartiteMatching {
    size_t sNumVtxs;
    size_t tNumVtxs;
    size_t numEdgs;
    std::vector<int32_t> sMate;
    std::vector<int32_t> tMate;
};

/* ---------- State: HKIState ---------- */

struct HKIState {
    std::vector<int32_t> sLevel;       // length sNumVtxs+1; sLevel[sNumVtxs] is the NIL sentinel
    std::vector<int32_t> sIdx;       // length sNumVtxs; relative offset within s's adjacency range, persistent within a phase
    std::vector<int32_t> sPrcbStk;       // length sNumVtxs; DFS stack (s-vertices)
    std::vector<int32_t> tPrcbStk;       // length sNumVtxs; t-vertex chosen at each depth
    int32_t stkTop;
};

/* ---------- BipartiteGraph construction ---------- */

BipartiteGraph buildBipartiteGraph(size_t sNumVtxs, size_t tNumVtxs,
                                   const std::vector<std::pair<int32_t,int32_t>>& edges) {
    std::vector<std::vector<int32_t>> sTmp(sNumVtxs);
    std::vector<std::vector<int32_t>> tTmp(tNumVtxs);
    for (auto& e : edges) {
        int32_t s = e.first, t = e.second;
        if (s >= 0 && static_cast<size_t>(s) < sNumVtxs && t >= 0 && static_cast<size_t>(t) < tNumVtxs) {
            sTmp[s].push_back(t);
            tTmp[t].push_back(s);
        }
    }
    for (size_t s = 0; s < sNumVtxs; s++) {
        std::sort(sTmp[s].begin(), sTmp[s].end());
        sTmp[s].erase(std::unique(sTmp[s].begin(), sTmp[s].end()), sTmp[s].end());
    }
    for (size_t t = 0; t < tNumVtxs; t++) {
        std::sort(tTmp[t].begin(), tTmp[t].end());
        tTmp[t].erase(std::unique(tTmp[t].begin(), tTmp[t].end()), tTmp[t].end());
    }

    BipartiteGraph graph;
    graph.sNumVtxs = sNumVtxs;
    graph.tNumVtxs = tNumVtxs;

    graph.sIdx.assign(sNumVtxs + 1, 0);
    for (size_t s = 0; s < sNumVtxs; s++)
        graph.sIdx[s + 1] = graph.sIdx[s] + sTmp[s].size();
    graph.sAdj.resize(graph.sIdx[sNumVtxs]);
    for (size_t s = 0; s < sNumVtxs; s++)
        std::copy(sTmp[s].begin(), sTmp[s].end(), graph.sAdj.begin() + graph.sIdx[s]);

    graph.tIdx.assign(tNumVtxs + 1, 0);
    for (size_t t = 0; t < tNumVtxs; t++)
        graph.tIdx[t + 1] = graph.tIdx[t] + tTmp[t].size();
    graph.tAdj.resize(graph.tIdx[tNumVtxs]);
    for (size_t t = 0; t < tNumVtxs; t++)
        std::copy(tTmp[t].begin(), tTmp[t].end(), graph.tAdj.begin() + graph.tIdx[t]);

    graph.numEdgs = graph.sIdx[sNumVtxs];
    return graph;
}

/* ---------- BipartiteMatching construction ---------- */

BipartiteMatching emptyBipartiteMatching(const BipartiteGraph& graph) {
    BipartiteMatching matching;
    matching.sNumVtxs = graph.sNumVtxs;
    matching.tNumVtxs = graph.tNumVtxs;
    matching.numEdgs = 0;
    matching.sMate.assign(graph.sNumVtxs, NIL);
    matching.tMate.assign(graph.tNumVtxs, NIL);
    return matching;
}

/* ---------- Greedy initial matching: simple ---------- */

int32_t greedyInit(const BipartiteGraph& graph, BipartiteMatching& matching) {
    int32_t numEdgs = 0;
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] != NIL) continue;
        size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
        for (size_t k = sBegin; k < sEnd; k++) {
            int32_t t = graph.sAdj[k];
            if (matching.tMate[t] == NIL) {
                matching.sMate[s] = t;
                matching.tMate[t] = static_cast<int32_t>(s);
                numEdgs++;
                break;
            }
        }
    }
    matching.numEdgs += numEdgs;
    return numEdgs;
}

/* ---------- Greedy initial matching: min-degree ---------- */

int32_t greedyInitMd(const BipartiteGraph& graph, BipartiteMatching& matching) {
    int32_t numEdgs = 0;
    std::vector<int32_t> deg(graph.tNumVtxs, 0);
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
        for (size_t k = sBegin; k < sEnd; k++) deg[graph.sAdj[k]]++;
    }
    std::vector<int32_t> sOrder(graph.sNumVtxs);
    for (size_t s = 0; s < graph.sNumVtxs; s++) sOrder[s] = static_cast<int32_t>(s);
    /* Sort s-vertices in increasing order of degree, breaking ties by vertex label. */
    std::sort(sOrder.begin(), sOrder.end(), [&](int32_t s1, int32_t s2){
        size_t s1Deg = graph.sIdx[s1 + 1] - graph.sIdx[s1];
        size_t s2Deg = graph.sIdx[s2 + 1] - graph.sIdx[s2];
        return s1Deg < s2Deg || (s1Deg == s2Deg && s1 < s2);
    });
    for (int32_t s : sOrder) {
        if (matching.sMate[s] != NIL) continue;
        int32_t best = NIL, bestDeg = INT_MAX;
        size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
        for (size_t k = sBegin; k < sEnd; k++) {
            int32_t t = graph.sAdj[k];
            if (matching.tMate[t] == NIL && deg[t] < bestDeg) {
                best = t;
                bestDeg = deg[t];
            }
        }
        if (best != NIL) {
            matching.sMate[s] = best;
            matching.tMate[best] = s;
            numEdgs++;
        }
    }
    matching.numEdgs += numEdgs;
    return numEdgs;
}

/* ---------- HK BFS ---------- */

static bool bfs(const BipartiteGraph& graph, const BipartiteMatching& matching, HKIState& state) {
    std::vector<int32_t> sPrcbQue(graph.sNumVtxs);
    int32_t queHead = 0, queTail = 0;

    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] == NIL) { state.sLevel[s] = 0; sPrcbQue[queTail++] = static_cast<int32_t>(s); }
        else state.sLevel[s] = INT_MAX;
    }
    state.sLevel[graph.sNumVtxs] = INT_MAX;  /* NIL sentinel */

    while (queHead < queTail) {
        int32_t s = sPrcbQue[queHead++];
        if (state.sLevel[s] < state.sLevel[graph.sNumVtxs]) {
            size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
            for (size_t k = sBegin; k < sEnd; k++) {
                int32_t t = graph.sAdj[k];
                int32_t ss = (matching.tMate[t] == NIL) ? static_cast<int32_t>(graph.sNumVtxs) : matching.tMate[t];
                if (state.sLevel[ss] == INT_MAX) {
                    state.sLevel[ss] = state.sLevel[s] + 1;
                    if (matching.tMate[t] != NIL) sPrcbQue[queTail++] = matching.tMate[t];
                }
            }
        }
    }
    return state.sLevel[graph.sNumVtxs] != INT_MAX;
}

/*
 * DFS: iterative with edge index.
 *
 * state.sIdx[s] is an offset WITHIN s's adjacency range [graph.sIdx[s], graph.sIdx[s+1]).
 * So the "current candidate edge" is graph.sAdj[graph.sIdx[s] + state.sIdx[s]].
 */
static bool dfs(int32_t sFirst, const BipartiteGraph& graph,
                BipartiteMatching& matching, HKIState& state) {
    state.stkTop = 0;
    state.sPrcbStk[state.stkTop++] = sFirst;

    while (state.stkTop > 0) {
        int32_t s = state.sPrcbStk[state.stkTop - 1];
        size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
        int32_t sNumEdgs = static_cast<int32_t>(sEnd - sBegin);

        bool pushed = false;
        while (state.sIdx[s] < sNumEdgs) {
            int32_t t = graph.sAdj[sBegin + state.sIdx[s]];
            int32_t ss = (matching.tMate[t] == NIL) ? static_cast<int32_t>(graph.sNumVtxs) : matching.tMate[t];
            if (state.sLevel[ss] != state.sLevel[s] + 1) {
                state.sIdx[s]++;
                continue;
            }

            state.tPrcbStk[state.stkTop - 1] = t;
            state.sIdx[s]++;

            if (matching.tMate[t] == NIL) {
                /* Found augmenting path — augment all the way back */
                for (int32_t k = state.stkTop - 1; k >= 0; k--) {
                    matching.tMate[state.tPrcbStk[k]] = state.sPrcbStk[k];
                    matching.sMate[state.sPrcbStk[k]] = state.tPrcbStk[k];
                }
                return true;
            }

            state.sPrcbStk[state.stkTop++] = matching.tMate[t];
            pushed = true;
            break;
        }

        if (!pushed) {
            state.sLevel[s] = INT_MAX;
            state.stkTop--;
        }
    }
    return false;
}

/* ---------- Top-level Hopcroft-Karp Iterative ---------- */

int32_t hkIterative(const BipartiteGraph& graph, BipartiteMatching& matching) {
    HKIState state;
    state.sLevel.assign(graph.sNumVtxs + 1, 0);
    state.sIdx.assign(graph.sNumVtxs, 0);
    state.sPrcbStk.assign(graph.sNumVtxs, 0);
    state.tPrcbStk.assign(graph.sNumVtxs, 0);
    state.stkTop = 0;

    int32_t numPhases = 0;
    int32_t newEdgs = 0;
    while (bfs(graph, matching, state)) {
        numPhases++;
        for (size_t s = 0; s < graph.sNumVtxs; s++) state.sIdx[s] = 0;
        for (size_t s = 0; s < graph.sNumVtxs; s++) {
            if (matching.sMate[s] == NIL && dfs(static_cast<int32_t>(s), graph, matching, state)) {
                newEdgs++;
            }
        }
    }
    matching.numEdgs += newEdgs;
    return numPhases;
}

/* ---------- Validation ---------- */

void validateBipartiteMatching(const BipartiteGraph& graph, const BipartiteMatching& matching) {
    int32_t errors = 0;
    int32_t matchedS = 0, matchedT = 0;

    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] != NIL) {
            matchedS++;
            int32_t t = matching.sMate[s];
            if (t < 0 || static_cast<size_t>(t) >= graph.tNumVtxs) {
                fprintf(stderr, "ERROR: sMate[%zu] = %d out of range\n", s, t);
                errors++;
            } else if (matching.tMate[t] != static_cast<int32_t>(s)) {
                fprintf(stderr, "ERROR: sMate[%zu]=%d but tMate[%d]=%d\n", s, t, t, matching.tMate[t]);
                errors++;
            } else {
                size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
                if (!std::binary_search(graph.sAdj.begin() + sBegin,
                                        graph.sAdj.begin() + sEnd, t)) {
                    fprintf(stderr, "ERROR: edge (%zu,%d) not in graph\n", s, t);
                    errors++;
                }
            }
        }
    }
    for (size_t t = 0; t < graph.tNumVtxs; t++) {
        if (matching.tMate[t] != NIL) matchedT++;
    }

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %zu\n", matching.numEdgs);
    printf("S matched: %d, T matched: %d\n", matchedS, matchedT);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

/* ---------- Main ---------- */

int main(int argc, char* argv[]) {
    printf("Hopcroft-Karp Iterative Algorithm - C++ Implementation (CSR)\n");
    printf("==============================================================\n\n");

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
        int s, t;
        if (fscanf(f, "%d %d", &s, &t) != 2) break;
        edges.push_back({static_cast<int32_t>(s), static_cast<int32_t>(t)});
    }
    fclose(f);

    printf("Graph: %d s-vertices, %d t-vertices, %zu edges\n",
           sNumVtxs, tNumVtxs, edges.size());

    BipartiteGraph graph = buildBipartiteGraph(sNumVtxs, tNumVtxs, edges);
    BipartiteMatching matching = emptyBipartiteMatching(graph);

    auto t0 = std::chrono::high_resolution_clock::now();

    int32_t greedySize = 0;
    if (greedyMode == 1) greedySize = greedyInit(graph, matching);
    else if (greedyMode == 2) greedySize = greedyInitMd(graph, matching);

    int32_t numPhases = hkIterative(graph, matching);

    auto t1 = std::chrono::high_resolution_clock::now();

    validateBipartiteMatching(graph, matching);

    printf("Phases: %d\n", numPhases);
    printf("Matching size: %zu\n", matching.numEdgs);
    if (greedyMode > 0) {
        printf("Greedy init size: %d\n", greedySize);
        if (matching.numEdgs > 0)
            printf("Greedy/Final: %.2f%%\n", 100.0 * greedySize / matching.numEdgs);
        else
            printf("Greedy/Final: NA\n");
    }
    printf("Time: %ld ms\n",
           static_cast<long>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()));
    return 0;
}

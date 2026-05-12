/*
 * Hopcroft-Karp Bidirectional, Iterative (HKB1) - O(E sqrt(V)) Maximum Bipartite Matching
 *
 * Variant 1 ("G2 style"): two parallel HK-style BFSes, one rooted at free S
 * vertices and one rooted at free T vertices, sharing a single FIFO. Each side
 * has its own level array with its own NIL sentinel. A side "fires" when its
 * BFS reaches a free vertex on the opposite partition; that pins its sentinel
 * and fixes the AP length on that side. DFS extraction runs from whichever
 * side(s) fired, exactly as in plain HK on the corresponding level array.
 *
 * This variant does strictly more BFS work than plain HK (twice the sources,
 * two level arrays maintained) and the same DFS work; it is not expected to
 * beat plain HK except under heavy free-vertex skew. It is a stepping stone
 * toward the MV-style variant (HKB2), where BFS and DDFS interleave per level
 * and the AP walk is halved by middle-out search at bridges.
 *
 * CSR adjacency, three-object architecture (BipartiteGraph / BipartiteMatching
 * / HKB1IState), top-level returns numPhases.
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
constexpr int8_t  SIDE_S = 0;
constexpr int8_t  SIDE_T = 1;

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

/* ---------- State: HKB1IState ----------
 *
 * sLevel, tLevel each have an extra slot at the end for the NIL sentinel.
 * sLevel[sNumVtxs] is pinned when a free T is reached from the S-side BFS;
 * tLevel[tNumVtxs] is pinned when a free S is reached from the T-side BFS.
 */

struct HKB1IState {
    std::vector<int32_t> sLevel;     // length sNumVtxs + 1; sLevel[sNumVtxs] = NIL sentinel
    std::vector<int32_t> tLevel;     // length tNumVtxs + 1; tLevel[tNumVtxs] = NIL sentinel
    std::vector<int32_t> sIdx;       // length sNumVtxs; relative offset within s's adjacency range
    std::vector<int32_t> tIdx;       // length tNumVtxs; relative offset within t's adjacency range
    std::vector<int32_t> sPrcbStk;   // length max(sNumVtxs, tNumVtxs); DFS stack (s-vertices)
    std::vector<int32_t> tPrcbStk;   // length max(sNumVtxs, tNumVtxs); t-vertex chosen at each depth
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

/* ---------- HKB1 BFS ----------
 *
 * Two parallel HK-style BFSes sharing a single FIFO. S-pops walk
 * S -> T (non-matching) -> S (matching follow); T-pops walk
 * T -> S (non-matching) -> T (matching follow). Each side has its own
 * NIL sentinel: sLevel[sNumVtxs] pinned when free T reached from S-side,
 * tLevel[tNumVtxs] pinned when free S reached from T-side.
 *
 * Returns true if either side fired (some AP exists).
 */

static bool bfs(const BipartiteGraph& graph, const BipartiteMatching& matching, HKB1IState& state) {
    /* FIFO holds (vtx, side) pairs via parallel arrays. */
    std::vector<int32_t> vtxQue(graph.sNumVtxs + graph.tNumVtxs);
    std::vector<int8_t>  sideQue(graph.sNumVtxs + graph.tNumVtxs);
    int32_t queHead = 0, queTail = 0;

    /* Init: free S at sLevel 0, free T at tLevel 0. */
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] == NIL) {
            state.sLevel[s] = 0;
            vtxQue[queTail] = static_cast<int32_t>(s);
            sideQue[queTail++] = SIDE_S;
        } else {
            state.sLevel[s] = INT_MAX;
        }
    }
    for (size_t t = 0; t < graph.tNumVtxs; t++) {
        if (matching.tMate[t] == NIL) {
            state.tLevel[t] = 0;
            vtxQue[queTail] = static_cast<int32_t>(t);
            sideQue[queTail++] = SIDE_T;
        } else {
            state.tLevel[t] = INT_MAX;
        }
    }
    state.sLevel[graph.sNumVtxs] = INT_MAX;  /* S-side NIL sentinel */
    state.tLevel[graph.tNumVtxs] = INT_MAX;  /* T-side NIL sentinel */

    while (queHead < queTail) {
        int32_t v = vtxQue[queHead];
        int8_t  side = sideQue[queHead];
        queHead++;

        if (side == SIDE_S) {
            int32_t s = v;
            /* Gate: don't expand past the S-side's pinned ceiling. */
            if (state.sLevel[s] >= state.sLevel[graph.sNumVtxs]) continue;

            size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
            for (size_t k = sBegin; k < sEnd; k++) {
                int32_t t = graph.sAdj[k];
                int32_t ss = (matching.tMate[t] == NIL) ? static_cast<int32_t>(graph.sNumVtxs) : matching.tMate[t];
                if (state.sLevel[ss] == INT_MAX) {
                    state.sLevel[ss] = state.sLevel[s] + 1;
                    if (matching.tMate[t] != NIL) {
                        vtxQue[queTail] = matching.tMate[t];
                        sideQue[queTail++] = SIDE_S;
                    }
                }
            }
        } else {
            int32_t t = v;
            if (state.tLevel[t] >= state.tLevel[graph.tNumVtxs]) continue;

            size_t tBegin = graph.tIdx[t], tEnd = graph.tIdx[t + 1];
            for (size_t k = tBegin; k < tEnd; k++) {
                int32_t s = graph.tAdj[k];
                int32_t tt = (matching.sMate[s] == NIL) ? static_cast<int32_t>(graph.tNumVtxs) : matching.sMate[s];
                if (state.tLevel[tt] == INT_MAX) {
                    state.tLevel[tt] = state.tLevel[t] + 1;
                    if (matching.sMate[s] != NIL) {
                        vtxQue[queTail] = matching.sMate[s];
                        sideQue[queTail++] = SIDE_T;
                    }
                }
            }
        }
    }

    bool sFired = state.sLevel[graph.sNumVtxs] != INT_MAX;
    bool tFired = state.tLevel[graph.tNumVtxs] != INT_MAX;
    return sFired || tFired;
}

/*
 * S-side DFS: identical to plain HK's DFS, walking S -> T via non-matching
 * edges where tMate[t] (mapped to NIL via the sNumVtxs slot) sits at
 * sLevel[s] + 1. Augments on free T discovery. state.sIdx[s] resumes
 * partial scans within a phase.
 */

static bool dfsFromS(int32_t sFirst, const BipartiteGraph& graph,
                     BipartiteMatching& matching, HKB1IState& state) {
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

/*
 * T-side DFS: mirror image of dfsFromS. Walks T -> S via non-matching edges
 * where sMate[s] (mapped to NIL via the tNumVtxs slot) sits at
 * tLevel[t] + 1. Augments on free S discovery. Stack convention reused:
 * sPrcbStk holds T-vertices, tPrcbStk holds S-vertices chosen at each depth.
 *
 * NB: variable names follow the partition role (s = T-vertex from this
 * function's perspective) rather than the partition tag, to keep the code
 * structurally identical to dfsFromS. Comments clarify.
 */

static bool dfsFromT(int32_t tFirst, const BipartiteGraph& graph,
                     BipartiteMatching& matching, HKB1IState& state) {
    state.stkTop = 0;
    state.sPrcbStk[state.stkTop++] = tFirst;  /* sPrcbStk now holds T-vertices */

    while (state.stkTop > 0) {
        int32_t t = state.sPrcbStk[state.stkTop - 1];
        size_t tBegin = graph.tIdx[t], tEnd = graph.tIdx[t + 1];
        int32_t tNumEdgs = static_cast<int32_t>(tEnd - tBegin);

        bool pushed = false;
        while (state.tIdx[t] < tNumEdgs) {
            int32_t s = graph.tAdj[tBegin + state.tIdx[t]];
            int32_t tt = (matching.sMate[s] == NIL) ? static_cast<int32_t>(graph.tNumVtxs) : matching.sMate[s];
            if (state.tLevel[tt] != state.tLevel[t] + 1) {
                state.tIdx[t]++;
                continue;
            }

            state.tPrcbStk[state.stkTop - 1] = s;  /* tPrcbStk now holds S-vertices */
            state.tIdx[t]++;

            if (matching.sMate[s] == NIL) {
                /* Augmenting path found. Stack stores the AP as alternating
                 * (t, s) pairs walking from tFirst to free s. Each (t, s)
                 * pair is a non-matching edge to install. */
                for (int32_t k = state.stkTop - 1; k >= 0; k--) {
                    matching.tMate[state.sPrcbStk[k]] = state.tPrcbStk[k];
                    matching.sMate[state.tPrcbStk[k]] = state.sPrcbStk[k];
                }
                return true;
            }

            state.sPrcbStk[state.stkTop++] = matching.sMate[s];
            pushed = true;
            break;
        }

        if (!pushed) {
            state.tLevel[t] = INT_MAX;
            state.stkTop--;
        }
    }
    return false;
}

/* ---------- Top-level HKB1 ---------- */

int32_t hkb1IterativeMcm(const BipartiteGraph& graph, BipartiteMatching& matching) {
    HKB1IState state;
    state.sLevel.assign(graph.sNumVtxs + 1, 0);
    state.tLevel.assign(graph.tNumVtxs + 1, 0);
    state.sIdx.assign(graph.sNumVtxs, 0);
    state.tIdx.assign(graph.tNumVtxs, 0);
    size_t stkCap = std::max(graph.sNumVtxs, graph.tNumVtxs);
    state.sPrcbStk.assign(stkCap, 0);
    state.tPrcbStk.assign(stkCap, 0);
    state.stkTop = 0;

    int32_t numPhases = 0;
    int32_t newEdgs = 0;
    while (bfs(graph, matching, state)) {
        numPhases++;

        /* Both sides may have fired; we only run one DFS per phase. Running
         * both back-to-back on the same BFS labeling is unsafe: the second
         * DFS sees stale levels for vertices consumed by the first DFS's
         * augmentations, and could augment through already-matched vertices.
         *
         * Choose the side with the shorter AP length when both fired (more
         * APs available at the smaller length per phase, faster sqrt(n)
         * convergence). Ties: prefer S-side for parity with plain HK. The
         * other side's APs, if any, are picked up next phase.
         */
        bool sFired = state.sLevel[graph.sNumVtxs] != INT_MAX;
        bool tFired = state.tLevel[graph.tNumVtxs] != INT_MAX;
        bool useS;
        if (sFired && tFired) {
            useS = (state.sLevel[graph.sNumVtxs] <= state.tLevel[graph.tNumVtxs]);
        } else {
            useS = sFired;
        }

        if (useS) {
            for (size_t s = 0; s < graph.sNumVtxs; s++) state.sIdx[s] = 0;
            for (size_t s = 0; s < graph.sNumVtxs; s++) {
                if (matching.sMate[s] == NIL && dfsFromS(static_cast<int32_t>(s), graph, matching, state)) {
                    newEdgs++;
                }
            }
        } else {
            for (size_t t = 0; t < graph.tNumVtxs; t++) state.tIdx[t] = 0;
            for (size_t t = 0; t < graph.tNumVtxs; t++) {
                if (matching.tMate[t] == NIL && dfsFromT(static_cast<int32_t>(t), graph, matching, state)) {
                    newEdgs++;
                }
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
    printf("Hopcroft-Karp Bidirectional Iterative (HKB1) - C++ Implementation (CSR)\n");
    printf("=========================================================================\n\n");

    if (argc < 2) { printf("Usage: %s <filename> [--greedy|--greedy-md]\n", argv[0]); return 1; }
    int32_t greedyMode = 0;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedyMode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedyMode = 2;
    }

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open file: %s\n", argv[1]); return 1; }

    int sNumVtxs, tNumVtxs, numEdgs;
    if (fscanf(f, "%d %d %d", &sNumVtxs, &tNumVtxs, &numEdgs) != 3) {
        fprintf(stderr, "Bad header\n"); fclose(f); return 1;
    }

    std::vector<std::pair<int32_t,int32_t>> edges;
    edges.reserve(numEdgs);
    for (int i = 0; i < numEdgs; i++) {
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

    int32_t numPhases = hkb1IterativeMcm(graph, matching);

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

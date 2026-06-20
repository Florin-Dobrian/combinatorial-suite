/*
 * Hopcroft-Karp Bidirectional, Iterative (HKB1) - O(E sqrt(V)) Maximum Bipartite Matching
 *
 * Variant 1 ("BG2"): bidirectional bucket-PQ keyed by BFS level. No FIFO.
 * Free vertices on BOTH partitions are seeded as EVEN tree roots at level 0,
 * with S-roots and T-roots growing distinct trees. Edges are scheduled into
 * buckets at the level where they become "ready":
 *
 *   EVEN-UNLABELED edge (will GROW):  bucket[level(EVEN) + 1].
 *   EVEN-EVEN cross-tree (a BRIDGE):  bucket[level(u) + level(v)].
 *
 * Bucket index for a length-(2k+1) AP equals k. Draining bucket[k] either
 * extends the trees by one level (GROW) or records bridges. The first bucket
 * with any bridge fixes the SAP length for the phase.
 *
 * Genuine meet-in-the-middle search: each tree grows to AP_length/2 before
 * bridges fire, vs HKB0 where each independent FIFO BFS walks the full AP
 * length. Augmentation is plain single-walker DFS, forward DFS in the G2
 * sense — no DDFS red/green (that's HKB2's specialization).
 *
 * This is "G2 minus blossoms minus duals" specialized to bipartite-MCM. The
 * Δ-clock and y-duals of G2 fall away under the 0/2-encoding (level = Δ/2,
 * y ≡ 1 by construction), leaving only the structural skeleton: bucket-PQ,
 * GROW, bridge detection.
 *
 * Bridge endpoints are walked back to free roots via per-ODD-vertex parent
 * pointers captured at GROW time. Augmentation is multi-AP per phase, with
 * a visited[] flag to keep APs vertex-disjoint.
 *
 * CSR adjacency, three-object architecture (BipartiteGraph / BipartiteMatching
 * / HKB1IState), top-level returns numPhases.
 */

#include <cstdio>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <chrono>
#include <climits>

constexpr int32_t NIL = -1;

constexpr int8_t  LBL_UNLABELED = 0;
constexpr int8_t  LBL_EVEN      = 1;
constexpr int8_t  LBL_ODD       = 2;

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
 * Per-vertex search state:
 *   sLabel/tLabel  - LBL_UNLABELED / LBL_EVEN / LBL_ODD per partition.
 *   sLevel/tLevel  - BFS level in the tree the vertex belongs to.
 *   sOddParent/tOddParent - for ODD vertices, the EVEN that GROWed them
 *                           (captured at GROW time; used to walk back to root).
 *   sVisited/tVisited - per-augment-phase flag for vertex-disjoint AP extraction.
 *
 * Bucket-PQ:
 *   buckets[d] - list of (s, t) pairs for edges that become ready at level d.
 *                One entry per scheduling event (so an edge may appear twice
 *                if both endpoints schedule it; duplicate-bridge detection
 *                is handled lazily in the augment phase via visited[]).
 *
 * bridges - list of (s, t) pairs found in the first non-empty bridge bucket.
 *           Multiple bridges may form a maximal vertex-disjoint set of SAPs.
 *
 * DFS scratch:
 *   sPath, tPath  - S-tree side walk: bridge S-endpoint → ... → free S root.
 *   sPath2, tPath2 - T-tree side walk: bridge T-endpoint → ... → free T root.
 */

struct HKB1IState {
    std::vector<int8_t>  sLabel;
    std::vector<int8_t>  tLabel;
    std::vector<int32_t> sLevel;
    std::vector<int32_t> tLevel;
    std::vector<int32_t> sOddParent;
    std::vector<int32_t> tOddParent;
    std::vector<bool>    sVisited;
    std::vector<bool>    tVisited;

    std::vector<std::vector<std::pair<int32_t, int32_t>>> buckets;
    std::vector<std::pair<int32_t, int32_t>> bridges;

    std::vector<int32_t> sPath, tPath;
    std::vector<int32_t> sPath2, tPath2;
    int32_t committedD;          /* bucket level the bridges fired at */
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

    BipartiteGraph g;
    g.sNumVtxs = sNumVtxs;
    g.tNumVtxs = tNumVtxs;
    g.sIdx.assign(sNumVtxs + 1, 0);
    g.tIdx.assign(tNumVtxs + 1, 0);
    for (size_t s = 0; s < sNumVtxs; s++) g.sIdx[s + 1] = g.sIdx[s] + sTmp[s].size();
    for (size_t t = 0; t < tNumVtxs; t++) g.tIdx[t + 1] = g.tIdx[t] + tTmp[t].size();
    g.sAdj.reserve(g.sIdx[sNumVtxs]);
    g.tAdj.reserve(g.tIdx[tNumVtxs]);
    for (size_t s = 0; s < sNumVtxs; s++) for (int32_t t : sTmp[s]) g.sAdj.push_back(t);
    for (size_t t = 0; t < tNumVtxs; t++) for (int32_t s : tTmp[t]) g.tAdj.push_back(s);
    g.numEdgs = g.sAdj.size();
    return g;
}

BipartiteMatching emptyBipartiteMatching(const BipartiteGraph& g) {
    BipartiteMatching m;
    m.sNumVtxs = g.sNumVtxs;
    m.tNumVtxs = g.tNumVtxs;
    m.numEdgs = 0;
    m.sMate.assign(g.sNumVtxs, NIL);
    m.tMate.assign(g.tNumVtxs, NIL);
    return m;
}

/* ---------- Greedy initialization ---------- */

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

/* ---------- HKB1 search: bidirectional bucket-PQ by level ---------- */

static inline void scheduleS(int32_t s,
                             const BipartiteGraph& graph,
                             const BipartiteMatching& matching,
                             HKB1IState& state) {
    /* Schedule outgoing unmatched edges of a newly-EVEN S-vertex `s`. */
    int32_t levelS = state.sLevel[s];
    int32_t matedT = matching.sMate[s];   /* NIL for free S */
    size_t sBegin = graph.sIdx[s], sEnd = graph.sIdx[s + 1];
    for (size_t k = sBegin; k < sEnd; k++) {
        int32_t t = graph.sAdj[k];
        if (t == matedT) continue;        /* skip matched edge */
        int8_t tLbl = state.tLabel[t];
        if (tLbl == LBL_ODD) continue;    /* already in tree on T-side */
        int32_t target = (tLbl == LBL_EVEN)
                       ? (levelS + state.tLevel[t])    /* EVEN-EVEN cross-tree bridge */
                       : (levelS + 1);                 /* EVEN-UNLABELED GROW */
        while (static_cast<int32_t>(state.buckets.size()) <= target) state.buckets.emplace_back();
        state.buckets[target].emplace_back(s, t);
    }
}

static inline void scheduleT(int32_t t,
                             const BipartiteGraph& graph,
                             const BipartiteMatching& matching,
                             HKB1IState& state) {
    /* Schedule outgoing unmatched edges of a newly-EVEN T-vertex `t`. */
    int32_t levelT = state.tLevel[t];
    int32_t matedS = matching.tMate[t];   /* NIL for free T */
    size_t tBegin = graph.tIdx[t], tEnd = graph.tIdx[t + 1];
    for (size_t k = tBegin; k < tEnd; k++) {
        int32_t s = graph.tAdj[k];
        if (s == matedS) continue;
        int8_t sLbl = state.sLabel[s];
        if (sLbl == LBL_ODD) continue;
        int32_t target = (sLbl == LBL_EVEN)
                       ? (state.sLevel[s] + levelT)
                       : (levelT + 1);
        while (static_cast<int32_t>(state.buckets.size()) <= target) state.buckets.emplace_back();
        state.buckets[target].emplace_back(s, t);
    }
}

static bool search(const BipartiteGraph& graph,
                   const BipartiteMatching& matching,
                   HKB1IState& state) {
    /* Reset per-iteration state. */
    std::fill(state.sLabel.begin(), state.sLabel.end(), LBL_UNLABELED);
    std::fill(state.tLabel.begin(), state.tLabel.end(), LBL_UNLABELED);
    std::fill(state.sLevel.begin(), state.sLevel.end(), NIL);
    std::fill(state.tLevel.begin(), state.tLevel.end(), NIL);
    state.buckets.clear();
    state.bridges.clear();

    /* SEED: free S as EVEN level 0 (will grow the S-tree); free T as EVEN
     * level 0 (will grow the T-tree). In bipartite the tree affiliation
     * follows from the label-and-partition: EVEN-S is always in the S-tree,
     * EVEN-T always in the T-tree — no explicit tree-side bookkeeping needed.
     * Label first, schedule second — so schedule sees consistent labels. */
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] == NIL) {
            state.sLabel[s] = LBL_EVEN;
            state.sLevel[s] = 0;
        }
    }
    for (size_t t = 0; t < graph.tNumVtxs; t++) {
        if (matching.tMate[t] == NIL) {
            state.tLabel[t] = LBL_EVEN;
            state.tLevel[t] = 0;
        }
    }
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] == NIL) scheduleS(static_cast<int32_t>(s), graph, matching, state);
    }
    for (size_t t = 0; t < graph.tNumVtxs; t++) {
        if (matching.tMate[t] == NIL) scheduleT(static_cast<int32_t>(t), graph, matching, state);
    }

    /* Drain buckets in level-increasing order; stop at first level with bridges.
     * FIFO within a bucket via qi++ walk: mirrors g2_csr.cpp's drain discipline
     * (g2 lines 386-388). Mid-drain pushes (a GROW spawning bridges to free
     * EVENs into the current bucket) sit at the tail and get processed after
     * older entries, preserving BFS-order within the bucket. */
    for (int32_t d = 0; d < static_cast<int32_t>(state.buckets.size()); d++) {
        bool foundThisLevel = false;
        size_t qi = 0;
        while (qi < state.buckets[d].size()) {
            auto edge = state.buckets[d][qi++];
            int32_t s = edge.first, t = edge.second;
            int8_t sLbl = state.sLabel[s], tLbl = state.tLabel[t];

            /* Doubly-booked or partner-already-in-tree: skip. */
            if (sLbl == LBL_ODD || tLbl == LBL_ODD) continue;

            if (sLbl == LBL_EVEN && tLbl == LBL_EVEN) {
                /* Both EVEN -> bridge candidate. Record only when the edge is
                 * tight at THIS level, i.e. sLvl + tLvl == d, giving a length
                 * 2d+1 shortest AP. An early bridge (scheduled stale, with
                 * sLvl + tLvl > d) is not the shortest yet; skip it here and
                 * let it fire at its true level d' = sLvl + tLvl later. */
                if (state.sLevel[s] + state.tLevel[t] == d) {
                    state.bridges.emplace_back(s, t);
                    foundThisLevel = true;
                }
                continue;
            }

            /* GROW: one EVEN, one UNLABELED. The UNLABELED side has a mate
             * (free vertices were seeded EVEN at init). */
            if (sLbl == LBL_EVEN) {
                int32_t mateT = matching.tMate[t];
                if (mateT == NIL) continue;   /* defensive — shouldn't fire */
                state.tLabel[t] = LBL_ODD;
                state.tLevel[t] = d;
                state.tOddParent[t] = s;
                state.sLabel[mateT] = LBL_EVEN;
                state.sLevel[mateT] = d;
                scheduleS(mateT, graph, matching, state);
            } else {
                int32_t mateS = matching.sMate[s];
                if (mateS == NIL) continue;
                state.sLabel[s] = LBL_ODD;
                state.sLevel[s] = d;
                state.sOddParent[s] = t;
                state.tLabel[mateS] = LBL_EVEN;
                state.tLevel[mateS] = d;
                scheduleT(mateS, graph, matching, state);
            }
        }
        if (foundThisLevel) { state.committedD = d; break; }
    }

    return !state.bridges.empty();
}

/* ---------- Augment: HK-style forward DFS from each free S ----------
 *
 * Phase 1 built a bidirectional BFS structure:
 *   S-tree (rooted at free S): EVEN-S at sLvl 0,1,2,...; ODD-T at tLvl 1,2,...
 *   T-tree (rooted at free T): EVEN-T at tLvl 0,1,2,...; ODD-S at sLvl 1,2,...
 *   Bridge: an EVEN-EVEN cross-tree edge.
 *
 * An AP from free_S to free_T crosses exactly one bridge edge (the unique
 * tree-switch). It has structure:
 *
 *   free_S (sLvl 0) -> ... -> bridge_S (sLvl a) - bridge - bridge_T (tLvl b)
 *                                                       -> ... -> free_T (tLvl 0)
 *
 * Note: a and b are independent -- they don't have to be equal. "Early
 * bridges" scheduled when one endpoint was UNLABELED can have a+b > bucket
 * level. The AP is still valid as long as both halves of the walk land on
 * level-monotone tight edges.
 *
 * Single forward DFS function. Augmentation starts from S side (free S)
 * only -- the walk goes forward through both trees in one pass. Dispatch
 * on sLabel of the current S-vertex:
 *
 *   EVEN-S (in S-tree at sLvl L):
 *     For each unmatched adjacent t:
 *       - tLbl ODD, tLvl == L+1   -> S-tree forward (matched edge -> recurse)
 *       - tLbl EVEN               -> bridge crossing into T-tree
 *
 *   ODD-S (in T-tree at sLvl m, post-bridge):
 *     For each unmatched adjacent t:
 *       - tLbl EVEN, tLvl == m-1  -> T-tree forward (matched edge -> recurse)
 *
 * Each successful step flips the unmatched edge to matched on return.
 * visited[] persists HK-style.
 *
 * MSAP completeness: adjacency enumeration at each step provides multiple
 * routes through any vertex -- distinct bridges can yield vertex-disjoint
 * APs through different intermediate paths, matching HK / MV / HKB2.
 */

static bool dfsAP(int32_t s,
                  const BipartiteGraph& graph,
                  BipartiteMatching& matching,
                  HKB1IState& state) {
    if (state.sVisited[s]) return false;
    state.sVisited[s] = true;

    int32_t sLvl = state.sLevel[s];
    int8_t  sLbl = state.sLabel[s];
    size_t a_begin = graph.sIdx[s], a_end = graph.sIdx[s + 1];

    for (size_t k = a_begin; k < a_end; k++) {
        int32_t t = graph.sAdj[k];
        if (t == matching.sMate[s]) continue;   /* skip matched edge */
        if (state.tVisited[t]) continue;

        int8_t  tLbl = state.tLabel[t];
        int32_t tLvl = state.tLevel[t];
        bool admit = false;

        if (sLbl == LBL_EVEN) {
            /* S-tree side. */
            if (tLbl == LBL_ODD && tLvl == sLvl + 1) admit = true;  /* forward in S-tree */
            else if (tLbl == LBL_EVEN && sLvl + tLvl == state.committedD) admit = true;  /* bridge at committed level only */
        } else if (sLbl == LBL_ODD) {
            /* T-tree side (post-bridge). */
            if (tLbl == LBL_EVEN && tLvl == sLvl - 1) admit = true;
        }
        if (!admit) continue;

        state.tVisited[t] = true;
        int32_t next_s = matching.tMate[t];

        if (next_s == NIL) {
            /* t is free T -- AP found. */
            matching.sMate[s] = t;
            matching.tMate[t] = s;
            return true;
        }
        if (dfsAP(next_s, graph, matching, state)) {
            matching.sMate[s] = t;
            matching.tMate[t] = s;
            return true;
        }
    }
    return false;
}

static int32_t augment(const BipartiteGraph& graph,
                       BipartiteMatching& matching,
                       HKB1IState& state) {
    std::fill(state.sVisited.begin(), state.sVisited.end(), false);
    std::fill(state.tVisited.begin(), state.tVisited.end(), false);

    int32_t newEdgs = 0;
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] != NIL) continue;       /* not free */
        if (state.sVisited[s]) continue;
        if (state.sLabel[s] != LBL_EVEN) continue;    /* must be S-tree EVEN at sLvl 0 */
        if (dfsAP(static_cast<int32_t>(s), graph, matching, state)) {
            newEdgs++;
        }
    }
    return newEdgs;
}

/* ---------- Top-level HKB1 ---------- */

int32_t hkb1IterativeMcm(const BipartiteGraph& graph, BipartiteMatching& matching) {
    HKB1IState state;
    state.sLabel.assign(graph.sNumVtxs, LBL_UNLABELED);
    state.tLabel.assign(graph.tNumVtxs, LBL_UNLABELED);
    state.sLevel.assign(graph.sNumVtxs, NIL);
    state.tLevel.assign(graph.tNumVtxs, NIL);
    state.sOddParent.assign(graph.sNumVtxs, NIL);
    state.tOddParent.assign(graph.tNumVtxs, NIL);
    state.sVisited.assign(graph.sNumVtxs, false);
    state.tVisited.assign(graph.tNumVtxs, false);

    int32_t numPhases = 0;
    int32_t newEdgs = 0;
    while (search(graph, matching, state)) {
        numPhases++;
        newEdgs += augment(graph, matching, state);
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

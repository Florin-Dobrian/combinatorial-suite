/*
 * Gabow MCM with proper duals and edge IDs. CSR variant.
 * O(E * sqrt(V)) Maximum Cardinality Matching.
 *
 * Edge model: each undirected edge gets a unique ID with fixed
 * (src, tgt) endpoints. The G2State carries sVtx/tVtx arrays and
 * a CSR-by-edge-id (idx + edgAdj) so iteration sites loop
 * over edge IDs incident to a vertex. This supports the exact
 * Phase 2 edge semantics from LEDA:
 *   G.source(e) -> sVtx[e]
 *   G.target(e) -> tVtx[e]
 *   G.opposite(v, e) -> (sVtx[e]==v ? tVtx[e] : sVtx[e])
 *
 * The H-side structure (contractedInto) remains a vector-of-vectors
 * because it grows dynamically during each phase.
 *
 * All integers. No floating point in algorithm. No dependencies.
 *
 * Three-object architecture:
 *   GeneralGraph    -- const input (CSR adjacency by neighbor vertex)
 *   GeneralMatching -- output (mate array + size)
 *   G2State         -- algorithm scratch (everything else)
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <climits>
#include <chrono>
#include <utility>

static const int NIL = -1;
enum { EVEN = 0, ODD = 1, UNLABELED = 2 };

/* ====================================================================
 *                          GeneralGraph (input)
 * ==================================================================== */

struct GeneralGraph {
    size_t numVtxs;
    size_t numEdgs;                      // number of undirected edges
    std::vector<size_t>  idx;            // row-pointer, length numVtxs+1
    std::vector<int32_t> adj;            // length 2*numEdgs, neighbor vertex IDs
};

/* Build from a list of (u, v) edges. Each undirected edge is
 * recorded twice (once at each endpoint); per-vertex adjacency
 * rows are sorted and deduplicated explicitly so adj is always
 * sorted within each row. Self-loops and out-of-range endpoints
 * are discarded. */
GeneralGraph buildGeneralGraph(size_t numVtxs,
                               const std::vector<std::pair<int32_t,int32_t>>& edges) {
    std::vector<std::vector<int32_t>> tmp(numVtxs);
    for (auto& e : edges) {
        int32_t u = e.first, v = e.second;
        if (u >= 0 && static_cast<size_t>(u) < numVtxs &&
            v >= 0 && static_cast<size_t>(v) < numVtxs && u != v) {
            tmp[u].push_back(v);
            tmp[v].push_back(u);
        }
    }
    for (size_t v = 0; v < numVtxs; v++) {
        std::sort(tmp[v].begin(), tmp[v].end());
        tmp[v].erase(std::unique(tmp[v].begin(), tmp[v].end()), tmp[v].end());
    }

    GeneralGraph graph;
    graph.numVtxs = numVtxs;

    graph.idx.assign(numVtxs + 1, 0);
    for (size_t v = 0; v < numVtxs; v++)
        graph.idx[v + 1] = graph.idx[v] + tmp[v].size();
    graph.adj.resize(graph.idx[numVtxs]);
    for (size_t v = 0; v < numVtxs; v++)
        std::copy(tmp[v].begin(), tmp[v].end(), graph.adj.begin() + graph.idx[v]);

    /* Each undirected edge is stored twice (once per endpoint),
     * so the total adj length is 2 * numEdgs. */
    graph.numEdgs = graph.idx[numVtxs] / 2;
    return graph;
}

/* ====================================================================
 *                       GeneralMatching (output)
 * ==================================================================== */

struct GeneralMatching {
    std::vector<int> mate;
    size_t numEdgs;
};

GeneralMatching emptyGeneralMatching(const GeneralGraph& graph) {
    GeneralMatching matching;
    matching.mate.assign(graph.numVtxs, NIL);
    matching.numEdgs = 0;
    return matching;
}

/* ====================================================================
 *                          G2State (scratch)
 * ==================================================================== */

struct G2State {
    /* Edge storage: sVtx[eid], tVtx[eid] are fixed endpoints */
    std::vector<int32_t> sVtx, tVtx;
    /* edgAdj: vertex -> edge IDs incident to it. Uses graph.idx as
     * the row pointer (graph.idx[v]..graph.idx[v+1] gives the same
     * row range, since each undirected edge contributes once to
     * each endpoint's adjacency in both views). */
    std::vector<int32_t> edgAdj;

    /* ---- Phase 1: alternating tree ---- */
    std::vector<int32_t> label;
    std::vector<int32_t> parent;          /* parent vertex in alternating tree */
    std::vector<int32_t> sourceBridge, targetBridge;

    /* ---- Phase 1: dual checkpoints ---- */
    std::vector<int32_t> bd, bDelta;
    int32_t Delta;

    /* ---- Phase 1: blossom DSU (current) ---- */
    std::vector<int32_t> basePar;

    /* ---- Phase 1/2: blossom DSU (snapshotted across H-construction) ---- */
    std::vector<int32_t> dbasePar, dbaseRank;

    /* ---- Phase 1: bucket priority queue ---- */
    int32_t maxPq;
    std::vector<std::vector<int32_t>> L;  /* L[d] = edge IDs becoming tight at Delta=d */
    int32_t maxDeltaUsed;

    /* ---- Phase 1: augmenting-path-vs-blossom walking trick ---- */
    std::vector<double> path1, path2;
    double strue;

    /* ---- Phase 1: visited tracking ---- */
    std::vector<int32_t> treeNodes;
    std::vector<int32_t> prevTreeNodes;

    /* ---- Free vertex tracking ---- */
    std::vector<int32_t> free_vertices;
    bool freeListBuilt;

    /* ---- H minor ---- */
    std::vector<int32_t> rep;
    std::vector<int32_t> mateH;
    std::vector<bool>    isH;             /* per edge ID: is this edge in H? */
    std::vector<std::vector<int32_t>> contractedInto;

    /* ---- Phase 2: search on H ---- */
    std::vector<int32_t> labelH;
    std::vector<int32_t> parentH;         /* parentH[uh] = edge ID */
    std::vector<int32_t> bridgeH;         /* bridgeH[vh] = edge ID */
    std::vector<int32_t> dirH;            /* dirH[vh] = 1 or -1 */
    std::vector<int32_t> evenTimeH;
    int32_t tH;
};

G2State emptyG2State(const GeneralGraph& graph) {
    G2State state;
    state.Delta = 0;
    state.maxDeltaUsed = 0;
    state.strue = 0.0;
    state.freeListBuilt = false;
    state.tH = 0;

    /* Build sVtx/tVtx from the graph's adjacency. We pick each
     * undirected edge exactly once by iterating u < v. */
    state.sVtx.reserve(graph.numEdgs);
    state.tVtx.reserve(graph.numEdgs);
    for (size_t u = 0; u < graph.numVtxs; u++) {
        for (size_t j = graph.idx[u]; j < graph.idx[u + 1]; j++) {
            int32_t v = graph.adj[j];
            if (static_cast<int32_t>(u) < v) {
                state.sVtx.push_back(static_cast<int32_t>(u));
                state.tVtx.push_back(v);
            }
        }
    }

    /* Build edgAdj (vertex -> edge IDs), reusing graph.idx as the row
     * pointer (graph.idx[v]..graph.idx[v+1] gives the same row range
     * since each undirected edge contributes once to each endpoint). */
    state.edgAdj.resize(graph.idx[graph.numVtxs]);
    std::vector<size_t> tmp_pos(graph.idx.begin(),
                                graph.idx.begin() + graph.numVtxs);
    for (size_t i = 0; i < graph.numEdgs; i++) {
        state.edgAdj[tmp_pos[state.sVtx[i]]++] = static_cast<int32_t>(i);
        state.edgAdj[tmp_pos[state.tVtx[i]]++] = static_cast<int32_t>(i);
    }
    state.isH.assign(graph.numEdgs, false);

    state.label.assign(graph.numVtxs, UNLABELED);
    state.parent.assign(graph.numVtxs, NIL);
    state.sourceBridge.assign(graph.numVtxs, NIL);
    state.targetBridge.assign(graph.numVtxs, NIL);
    state.bd.assign(graph.numVtxs, 1);
    state.bDelta.assign(graph.numVtxs, 0);

    state.basePar.resize(graph.numVtxs);
    for (size_t i = 0; i < graph.numVtxs; i++) state.basePar[i] = static_cast<int32_t>(i);
    state.dbasePar.resize(graph.numVtxs);
    for (size_t i = 0; i < graph.numVtxs; i++) state.dbasePar[i] = static_cast<int32_t>(i);
    state.dbaseRank.assign(graph.numVtxs, 0);

    state.maxPq = graph.numVtxs / 2 + 2;
    state.L.resize(state.maxPq);
    state.path1.assign(graph.numVtxs, 0.0);
    state.path2.assign(graph.numVtxs, 0.0);

    state.rep.resize(graph.numVtxs);
    state.mateH.assign(graph.numVtxs, NIL);
    state.labelH.assign(graph.numVtxs, UNLABELED);
    state.parentH.assign(graph.numVtxs, NIL);
    state.bridgeH.assign(graph.numVtxs, NIL);
    state.dirH.assign(graph.numVtxs, 0);
    state.evenTimeH.assign(graph.numVtxs, 0);
    state.contractedInto.resize(graph.numVtxs);
    return state;
}

/* ====================================================================
 *                          Helpers (free functions)
 * ==================================================================== */

inline int opposite(const G2State& state, int v, int eid) {
    return state.sVtx[eid] == v ? state.tVtx[eid] : state.sVtx[eid];
}

inline int wght(const GeneralMatching& matching, const G2State& state, int eid) {
    return (matching.mate[state.sVtx[eid]] == state.tVtx[eid]) ? 2 : 0;
}

inline int findBase(G2State& state, int v) {
    while (state.basePar[v] != v) {
        state.basePar[v] = state.basePar[state.basePar[v]];
        v = state.basePar[v];
    }
    return v;
}

inline int findDbase(G2State& state, int v) {
    while (state.dbasePar[v] != v) {
        state.dbasePar[v] = state.dbasePar[state.dbasePar[v]];
        v = state.dbasePar[v];
    }
    return v;
}

inline void unionDbase(G2State& state, int a, int b) {
    a = findDbase(state, a); b = findDbase(state, b);
    if (a == b) return;
    if (state.dbaseRank[a] < state.dbaseRank[b]) std::swap(a, b);
    state.dbasePar[b] = a;
    if (state.dbaseRank[a] == state.dbaseRank[b]) state.dbaseRank[a]++;
}

inline void makeRepDbase(G2State& state, int v) {
    int r = findDbase(state, v);
    if (r != v) { state.dbasePar[r] = v; state.dbasePar[v] = v; }
}

/* dual(v): dual variable, computed on demand from label and bd/bDelta checkpoints */
inline int dual(G2State& state, int v) {
    int bv = findBase(state, v);
    if (state.label[bv] == UNLABELED) return 1;
    if (state.label[bv] == EVEN) return state.bd[v] - (state.Delta - state.bDelta[v]);
    return state.bd[v] + (state.Delta - state.bDelta[v]);
}

/* scanEdge: schedule edge into PQ at its projected tightness Delta. */
void scanEdge(const GeneralMatching& matching, G2State& state, int eid, int z) {
    int u = opposite(state, z, eid);
    if (matching.mate[u] == z || state.label[findBase(state, u)] == ODD) return;
    int p = dual(state, z) + dual(state, u);
    int tight_at;
    if (state.label[findBase(state, u)] == UNLABELED)
        tight_at = state.Delta + p;
    else
        tight_at = state.Delta + p / 2;
    if (tight_at >= 0 && tight_at < state.maxPq) {
        state.L[tight_at].push_back(eid);
        if (tight_at > state.maxDeltaUsed) state.maxDeltaUsed = tight_at;
    }
}

/* shrinkPath: contract a half-cycle into the blossom rooted at b. */
void shrinkPath(const GeneralGraph& graph, const GeneralMatching& matching, G2State& state,
                 int b, int x, int y,
                 std::vector<std::pair<int,int>>& dunions) {
    int v = findBase(state, x);
    while (v != b) {
        state.basePar[v] = b;
        dunions.push_back({v, b});
        v = matching.mate[v];
        state.basePar[v] = b;
        dunions.push_back({v, b});
        state.basePar[b] = b;
        state.sourceBridge[v] = x;
        state.targetBridge[v] = y;
        state.bd[v] = state.bd[v] + (state.Delta - state.bDelta[v]);
        state.bDelta[v] = state.Delta;
        {
            int s = graph.idx[v], e = graph.idx[v + 1];
            for (int j = s; j < e; j++)
                scanEdge(matching, state, state.edgAdj[j], v);
        }
        v = findBase(state, state.parent[v]);
    }
    dunions.push_back({b, b});
}

/* ---- Free-vertex list maintenance ---- */
void buildFreeList(const GeneralGraph& graph, const GeneralMatching& matching, G2State& state) {
    state.free_vertices.clear();
    for (size_t v = 0; v < graph.numVtxs; v++)
        if (matching.mate[v] == NIL) state.free_vertices.push_back(static_cast<int32_t>(v));
    state.freeListBuilt = true;
}

void updateFreeList(const GeneralMatching& matching, G2State& state) {
    int j = 0;
    for (int i = 0; i < static_cast<int>(state.free_vertices.size()); i++)
        if (matching.mate[state.free_vertices[i]] == NIL)
            state.free_vertices[j++] = state.free_vertices[i];
    state.free_vertices.resize(j);
}

/* ====================================================================
 *                             PHASE 1
 * ==================================================================== */
bool phase1(const GeneralGraph& graph, const GeneralMatching& matching, G2State& state) {
    state.Delta = 0;
    state.treeNodes.clear();

    /* Only clear L entries up to maxDeltaUsed from previous phase */
    int clear_limit = std::min(state.maxDeltaUsed + 1, state.maxPq);
    for (int i = 0; i < clear_limit; i++) state.L[i].clear();
    state.maxDeltaUsed = 0;
    std::vector<std::pair<int,int>> dunions;

    /* Reset only previous tree nodes. */
    for (int v : state.prevTreeNodes) {
        state.basePar[v] = v;
        state.dbasePar[v] = v;
        state.dbaseRank[v] = 0;
        state.label[v] = UNLABELED;
        state.parent[v] = NIL;
        state.sourceBridge[v] = NIL;
        state.targetBridge[v] = NIL;
        state.bd[v] = 1;
        state.bDelta[v] = 0;
        {
            int s = graph.idx[v], e = graph.idx[v + 1];
            for (int j = s; j < e; j++)
                state.isH[state.edgAdj[j]] = false;
        }
    }

    /* Build or update free vertex list */
    if (!state.freeListBuilt)
        buildFreeList(graph, matching, state);
    else
        updateFreeList(matching, state);

    /* Label free vertices EVEN, then scan */
    for (int v : state.free_vertices) {
        state.label[v] = EVEN;
        state.treeNodes.push_back(v);
    }
    for (int v : state.free_vertices) {
        int s = graph.idx[v], e = graph.idx[v + 1];
        for (int j = s; j < e; j++)
            scanEdge(matching, state, state.edgAdj[j], v);
    }

    bool found_sap = false;

    while (state.Delta <= state.maxDeltaUsed) {
        /* Skip empty levels */
        while (state.Delta <= state.maxDeltaUsed && state.L[state.Delta].empty()) state.Delta++;
        if (state.Delta > state.maxDeltaUsed) break;

        int qi = 0;
        while (qi < static_cast<int>(state.L[state.Delta].size())) {
            int eid = state.L[state.Delta][qi++];
            int x = state.sVtx[eid], y = state.tVtx[eid];

            /* Stale-entry guard */
            if (dual(state, x) + dual(state, y) != wght(matching, state, eid)) continue;

            if (state.label[findBase(state, x)] != EVEN) std::swap(x, y);
            if (y == matching.mate[x] || findBase(state, x) == findBase(state, y) ||
                state.label[findBase(state, y)] == ODD) continue;

            if (state.label[findBase(state, y)] == UNLABELED) {
                int z = matching.mate[y];
                state.bd[y] = 1; state.bDelta[y] = state.Delta;
                state.bd[z] = 1; state.bDelta[z] = state.Delta;
                state.parent[z] = y;
                state.parent[y] = x;
                state.label[y] = ODD;
                state.label[z] = EVEN;
                state.treeNodes.push_back(y);
                state.treeNodes.push_back(z);
                {
                    int s = graph.idx[z], e = graph.idx[z + 1];
                    for (int j = s; j < e; j++)
                        scanEdge(matching, state, state.edgAdj[j], z);
                }

            } else if (state.label[findBase(state, y)] == EVEN) {
                state.strue += 1.0;
                int hx = findBase(state, x), hy = findBase(state, y);
                state.path1[hx] = state.strue; state.path2[hy] = state.strue;
                int lca = NIL;
                while (true) {
                    if (state.path1[hy] == state.strue) { lca = hy; break; }
                    if (state.path2[hx] == state.strue) { lca = hx; break; }
                    bool hxr = (matching.mate[hx] == NIL || state.parent[matching.mate[hx]] == NIL);
                    bool hyr = (matching.mate[hy] == NIL || state.parent[matching.mate[hy]] == NIL);
                    if (hxr && hyr) break;
                    if (!hxr) { hx = findBase(state, state.parent[matching.mate[hx]]); state.path1[hx] = state.strue; }
                    if (!hyr) { hy = findBase(state, state.parent[matching.mate[hy]]); state.path2[hy] = state.strue; }
                }
                if (lca != NIL) {
                    shrinkPath(graph, matching, state, lca, x, y, dunions);
                    shrinkPath(graph, matching, state, lca, y, x, dunions);
                } else {
                    found_sap = true;
                }
            }
        }
        state.L[state.Delta].clear();

        if (found_sap) {
            /* Build H */
            for (int v : state.treeNodes) {
                state.contractedInto[findDbase(state, v)].push_back(v);
                state.mateH[v] = NIL;
            }
            /* Mark tight edges */
            for (int u : state.treeNodes) {
                int uh = findDbase(state, u);
                int s = graph.idx[u], e = graph.idx[u + 1];
                for (int j = s; j < e; j++) {
                    int eid = state.edgAdj[j];
                    int v = opposite(state, u, eid);
                    int vh = findDbase(state, v);
                    if (uh != vh && dual(state, u) + dual(state, v) == wght(matching, state, eid)) {
                        state.isH[eid] = true;
                        if (wght(matching, state, eid) == 2) {
                            state.mateH[uh] = vh;
                            state.mateH[vh] = uh;
                        }
                    }
                }
            }
            state.prevTreeNodes = state.treeNodes;
            return true;
        }

        for (auto& [a, b] : dunions) {
            if (a == b) makeRepDbase(state, a);
            else unionDbase(state, a, b);
        }
        dunions.clear();
        state.Delta++;
    }
    state.prevTreeNodes = state.treeNodes;
    return false;
}

/* ====================================================================
 *                             PHASE 2
 * ==================================================================== */

/* findApHG: recursive DFS on H. */
int findApHG(const GeneralGraph& graph, const GeneralMatching& matching, G2State& state, int vh) {
    (void)matching;
    for (int v : state.contractedInto[vh]) {
        int a_beg = graph.idx[v], a_end = graph.idx[v + 1];
        for (int j = a_beg; j < a_end; j++) {
            int eid = state.edgAdj[j];
            if (!state.isH[eid]) continue;
            int w = opposite(state, v, eid);
            int uh = state.rep[w];
            if (state.mateH[vh] == uh) continue;

            if (state.labelH[uh] == UNLABELED) {
                int muh = state.mateH[uh];
                if (muh == NIL) {
                    state.labelH[uh] = ODD;
                    state.parentH[uh] = eid;
                    return uh;
                }
                state.labelH[uh] = ODD;
                state.labelH[muh] = EVEN;
                state.parentH[uh] = eid;
                state.evenTimeH[muh] = state.tH++;
                int s = findApHG(graph, matching, state, muh);
                if (s != NIL) return s;

            } else {
                int bh = findDbase(state, vh);
                int zh = findDbase(state, uh);
                if (state.evenTimeH[bh] < state.evenTimeH[zh]) {
                    std::vector<int> tmp;
                    std::vector<int> endpoints;
                    while (zh != bh) {
                        endpoints.push_back(zh);
                        zh = state.mateH[zh];
                        endpoints.push_back(zh);
                        tmp.insert(tmp.begin(), zh);
                        int pe = state.parentH[zh];
                        zh = findDbase(state, state.rep[state.rep[state.sVtx[pe]] == zh ? state.tVtx[pe] : state.sVtx[pe]]);
                    }
                    for (int nd : endpoints) unionDbase(state, nd, bh);
                    makeRepDbase(state, bh);
                    for (int odd_node : tmp) {
                        state.bridgeH[odd_node] = eid;
                        state.dirH[odd_node] = (state.tVtx[eid] == v) ? 1 : -1;
                    }
                    for (int odd_node : tmp) {
                        int s = findApHG(graph, matching, state, odd_node);
                        if (s != NIL) return s;
                    }
                }
            }
        }
    }
    return NIL;
}

/* findPathInHG: trace augmenting path in H from vh to uh. */
void findPathInHG(const GeneralMatching& matching, G2State& state,
                     std::vector<int>& path, int vh, int uh) {
    if (vh == uh) return;
    if (state.labelH[vh] == EVEN) {
        int mvh = state.mateH[vh];
        int pe = state.parentH[mvh];
        path.push_back(pe);
        int next = state.rep[state.rep[state.sVtx[pe]] == mvh ? state.tVtx[pe] : state.sVtx[pe]];
        findPathInHG(matching, state, path, next, uh);
    } else {
        int be = state.bridgeH[vh];
        int mate_side, uh_side;
        if (state.dirH[vh] == 1) {
            mate_side = state.rep[state.sVtx[be]];
            uh_side   = state.rep[state.tVtx[be]];
        } else {
            mate_side = state.rep[state.tVtx[be]];
            uh_side   = state.rep[state.sVtx[be]];
        }
        int mt = (state.mateH[vh] != NIL) ? state.rep[state.mateH[vh]] : vh;
        findPathInHG(matching, state, path, mate_side, mt);
        path.push_back(be);
        findPathInHG(matching, state, path, uh_side, uh);
    }
}

/* findPathInG: unfold within Phase 1 blossom. */
void findPathInG(const GeneralMatching& matching, G2State& state,
                    std::vector<std::pair<int,int>>& pairs, int v, int u) {
    if (v == u) return;
    if (state.label[v] == EVEN) {
        pairs.push_back({matching.mate[v], state.parent[matching.mate[v]]});
        findPathInG(matching, state, pairs, state.parent[matching.mate[v]], u);
    } else {
        findPathInG(matching, state, pairs, state.sourceBridge[v], matching.mate[v]);
        pairs.push_back({state.sourceBridge[v], state.targetBridge[v]});
        findPathInG(matching, state, pairs, state.targetBridge[v], u);
    }
}

/* augmentG: unfold H-path edges to G, augment matching. */
void augmentG(GeneralMatching& matching, G2State& state,
              const std::vector<int>& h_edge_ids) {
    std::vector<std::pair<int,int>> pairs;
    for (int eid : h_edge_ids) {
        int u = state.sVtx[eid], v = state.tVtx[eid];
        pairs.push_back({u, v});
        findPathInG(matching, state, pairs, u, state.rep[u]);
        findPathInG(matching, state, pairs, v, state.rep[v]);
    }
    for (auto& [a, b] : pairs) {
        matching.mate[a] = b; matching.mate[b] = a;
    }
    matching.numEdgs++;
}

void phase2(const GeneralGraph& graph, GeneralMatching& matching, G2State& state) {
    state.tH = 0;
    for (int v : state.treeNodes) {
        state.rep[v] = findDbase(state, v);
        state.labelH[v] = UNLABELED;
        state.parentH[v] = NIL;
        state.bridgeH[v] = NIL;
        state.dirH[v] = 0;
        state.evenTimeH[v] = 0;
    }

    std::vector<std::vector<int>> all_paths;

    for (int vh : state.treeNodes) {
        if (vh != state.rep[vh]) continue;
        if (state.labelH[vh] == UNLABELED && state.mateH[vh] == NIL) {
            state.labelH[vh] = EVEN;
            state.evenTimeH[vh] = state.tH++;
            int found = findApHG(graph, matching, state, vh);
            if (found != NIL) {
                std::vector<int> path;
                int pe = state.parentH[found];
                path.push_back(pe);
                int next = state.rep[state.rep[state.sVtx[pe]] == found ? state.tVtx[pe] : state.sVtx[pe]];
                findPathInHG(matching, state, path, next, vh);
                all_paths.push_back(std::move(path));
            }
        }
    }

    for (auto& p : all_paths)
        augmentG(matching, state, p);

    for (int v : state.treeNodes) {
        state.contractedInto[v].clear();
        state.mateH[v] = NIL;
    }
}

/* ====================================================================
 *                          Greedy initialization
 * ==================================================================== */

/* ---------- Greedy initial matching: simple ---------- */
int32_t greedyInit(const GeneralGraph& graph, GeneralMatching& matching) {
    int32_t numEdgs = 0;
    for (size_t u = 0; u < graph.numVtxs; u++) {
        if (matching.mate[u] != NIL) continue;
        size_t uBegin = graph.idx[u], uEnd = graph.idx[u + 1];
        for (size_t k = uBegin; k < uEnd; k++) {
            int32_t v = graph.adj[k];
            if (matching.mate[v] == NIL) {
                matching.mate[u] = v;
                matching.mate[v] = static_cast<int32_t>(u);
                numEdgs++;
                break;
            }
        }
    }
    matching.numEdgs += numEdgs;
    return numEdgs;
}

/* ---------- Greedy initial matching: min-degree ---------- */
int32_t greedyInitMd(const GeneralGraph& graph, GeneralMatching& matching) {
    int32_t numEdgs = 0;
    std::vector<int32_t> deg(graph.numVtxs, 0);
    for (size_t u = 0; u < graph.numVtxs; u++) {
        deg[u] = static_cast<int32_t>(graph.idx[u + 1] - graph.idx[u]);
    }
    std::vector<int32_t> order(graph.numVtxs);
    for (size_t u = 0; u < graph.numVtxs; u++) order[u] = static_cast<int32_t>(u);
    /* Sort vertices in increasing order of degree, breaking ties by vertex label. */
    std::sort(order.begin(), order.end(), [&](int32_t u1, int32_t u2){
        return deg[u1] < deg[u2] || (deg[u1] == deg[u2] && u1 < u2);
    });
    for (int32_t u : order) {
        if (matching.mate[u] != NIL) continue;
        int32_t best = NIL, bestDeg = INT_MAX;
        size_t uBegin = graph.idx[u], uEnd = graph.idx[u + 1];
        for (size_t k = uBegin; k < uEnd; k++) {
            int32_t v = graph.adj[k];
            if (matching.mate[v] == NIL && deg[v] < bestDeg) {
                best = v;
                bestDeg = deg[v];
            }
        }
        if (best != NIL) {
            matching.mate[u] = best;
            matching.mate[best] = u;
            numEdgs++;
        }
    }
    matching.numEdgs += numEdgs;
    return numEdgs;
}

/* ====================================================================
 *                            Top-level driver
 * ==================================================================== */
int32_t g2Mcm(const GeneralGraph& graph, GeneralMatching& matching) {
    G2State state = emptyG2State(graph);
    int32_t numPhases = 0;
    while (true) {
        bool hasSap = phase1(graph, matching, state);
        if (!hasSap) break;
        phase2(graph, matching, state);
        numPhases++;
    }
    return numPhases;
}

/* ====================================================================
 *                            Validation + main
 * ==================================================================== */
void validateGeneralMatching(const GeneralGraph& graph, const GeneralMatching& matching) {
    int32_t errors = 0;
    int32_t numMatched = 0;

    for (size_t u = 0; u < graph.numVtxs; u++) {
        int32_t v = matching.mate[u];
        if (v == NIL) continue;
        numMatched++;
        if (v < 0 || static_cast<size_t>(v) >= graph.numVtxs) {
            fprintf(stderr, "ERROR: mate[%zu] = %d out of range\n", u, v);
            errors++;
        } else if (matching.mate[v] != static_cast<int32_t>(u)) {
            fprintf(stderr, "ERROR: mate[%zu]=%d but mate[%d]=%d\n", u, v, v, matching.mate[v]);
            errors++;
        } else if (static_cast<int32_t>(u) < v) {
            /* Check edge presence in graph (each undirected edge once, when u < v). */
            size_t uBegin = graph.idx[u], uEnd = graph.idx[u + 1];
            if (!std::binary_search(graph.adj.begin() + uBegin,
                                    graph.adj.begin() + uEnd, v)) {
                fprintf(stderr, "ERROR: edge (%zu,%d) not in graph\n", u, v);
                errors++;
            }
        }
    }

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %zu\n", matching.numEdgs);
    printf("Vertices matched: %d\n", numMatched);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

int main(int argc, char* argv[]) {
    printf("Gabow MCM (duals + edge IDs) - C++\n");
    printf("===================================\n\n");

    if (argc < 2) { printf("Usage: %s <filename> [--greedy|--greedy-md]\n", argv[0]); return 1; }
    int32_t greedyMode = 0;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedyMode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedyMode = 2;
    }

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    int numVtxs, numEdgs;
    if (fscanf(f, "%d %d", &numVtxs, &numEdgs) != 2) { fprintf(stderr, "Bad header\n"); fclose(f); return 1; }

    std::vector<std::pair<int32_t,int32_t>> edges;
    edges.reserve(numEdgs);
    for (int i = 0; i < numEdgs; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edges.push_back({static_cast<int32_t>(u), static_cast<int32_t>(v)});
    }
    fclose(f);

    printf("Graph: %d vertices, %zu edges\n", numVtxs, edges.size());

    GeneralGraph graph = buildGeneralGraph(numVtxs, edges);
    GeneralMatching matching = emptyGeneralMatching(graph);

    auto t0 = std::chrono::high_resolution_clock::now();

    int32_t greedySize = 0;
    if (greedyMode == 1) greedySize = greedyInit(graph, matching);
    else if (greedyMode == 2) greedySize = greedyInitMd(graph, matching);

    int32_t numPhases = g2Mcm(graph, matching);

    auto t1 = std::chrono::high_resolution_clock::now();

    validateGeneralMatching(graph, matching);

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

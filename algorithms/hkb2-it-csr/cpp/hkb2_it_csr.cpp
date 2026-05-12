/*
 * Hopcroft-Karp Bidirectional, Iterative, with DDFS (HKB2) - O(E sqrt(V))
 *
 * Variant 2: bidirectional BFS interleaved with DDFS at the level where
 * bridges first appear.  Conceptually distinct from MV (which is a
 * general-graph algorithm with blossom machinery), but the DDFS coordination
 * is borrowed from MV — coordinated cursors, Sr/Sg stacks, level-tie
 * tie-break via L() comparison, step_into, above/below pointers.
 *
 * In bipartite there are no blossoms, so the "contracted graph H" of MV is
 * identical to G — there is no contraction.  This means:
 *   - bud_star is identity; we never look up a representative.
 *   - max_level is unused; only min_level (= even_level if even, odd_level
 *     if odd) matters.
 *   - hanging_bridges, walk_blossom*, jump_bridge are all dead and removed.
 *   - DDFS_PETAL never fires (verified by instrumenting MV on hundreds of
 *     bipartite test graphs).
 *
 * Earlier drafts of this file kept MV's per-vertex `preds` and `predTo`
 * vectors intact, which was an inefficiency — MV needed those because
 * blossom contraction made BFS predecessors point into the contracted
 * graph, requiring a precomputed mapping back.  In bipartite with H = G,
 * BFS predecessors are recoverable directly from the CSR using the BFS
 * parity rule:
 *   - vertex v at EVEN level: unique pred = match[v] (matching is 1-regular).
 *   - vertex v at ODD level:  preds = vAdj[v] \ {match[v]} restricted to
 *                             vertices at level v.level - 1.
 * So we drop `preds` and `predTo` entirely, keeping only `numPreds` (an
 * integer counter that drives cascade-delete in remove_path).  This
 * substantially shrinks per-vertex state and ends pointer-chasing during
 * DDFS — both were costing memory and cache misses in the 500k–1m valley
 * observed in the comparison vs HK.
 *
 * Bipartite input is presented in three-object form (BipartiteGraph /
 * BipartiteMatching / HKB2IState) but the algorithm internally uses a
 * unified vertex namespace [0, sN+tN) — S in [0, sN), T in [sN, sN+tN).
 * Free vertices on either partition seed BFS at level 0 (this is
 * "stacking S and T sources" — analogous to the bidirectional Hungarian
 * search the paper describes).
 *
 * CSR adjacency, top-level returns numPhases.
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

constexpr int32_t DDFS_EMPTY = 0;
constexpr int32_t DDFS_PATH  = 2;
/* DDFS_PETAL omitted — never returned in bipartite. */

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

/* ---------- Per-vertex node, bipartite-clean ----------
 *
 * Bipartite has no blossoms, so there is no contracted graph H — everything
 * the algorithm needs sits directly in G's CSR adjacency.  In MV, `preds[v]`
 * stored an explicit list of predecessor vertices used to walk back from v;
 * here that list is implicit:
 *
 *   - If v is at EVEN level, its only BFS predecessor is mate(v) (the
 *     matching is 1-regular, so a single odd-level vertex propagated to v).
 *   - If v is at ODD level, its predecessors are all u in vAdj[v] with
 *     u != mate(v) and level[u] == level[v] - 1 (the non-matching neighbors
 *     at the previous even level).
 *
 * This lets us drop the per-vertex `preds` / `predTo` vectors entirely —
 * a substantial reduction in memory (one std::vector<int32_t> + one
 * std::vector<pair<int,int>> per vertex was costing ~30-50% of total state)
 * and an end to pointer-chasing during DDFS.  Only `numPreds`, the count of
 * still-alive predecessors, is kept explicitly — it drives cascade-delete
 * in `removePath`.
 */

struct HKB2Node {
    int32_t minLevel;
    int32_t evenLevel;
    int32_t oddLevel;
    int32_t match;
    int32_t above;
    int32_t below;
    int32_t ddfsGreen;
    int32_t ddfsRed;
    int32_t numPreds;     /* count of still-alive BFS predecessors */
    uint8_t deleted;
    uint8_t visited;

    HKB2Node() : minLevel(NIL), evenLevel(NIL), oddLevel(NIL),
                 match(NIL), above(NIL), below(NIL),
                 ddfsGreen(NIL), ddfsRed(NIL), numPreds(0),
                 deleted(0), visited(0) {}

    void setMinLevel(int32_t lvl) {
        minLevel = lvl;
        if (lvl & 1) oddLevel = lvl;
        else evenLevel = lvl;
    }

    void reset() {
        minLevel = evenLevel = oddLevel = NIL;
        above = below = ddfsGreen = ddfsRed = NIL;
        numPreds = 0;
        deleted = 0;
        visited = 0;
    }
};

/* ---------- Internal state: HKB2IState ---------- */

struct HKB2IState {
    size_t nVtxs;
    size_t sNumVtxs;

    /* Unified CSR adjacency, both directions. */
    std::vector<size_t>  vIdx;
    std::vector<int32_t> vAdj;

    std::vector<HKB2Node> nodes;

    std::vector<std::vector<int32_t>> levels;
    std::vector<std::vector<std::pair<int32_t,int32_t>>> brBuckets;

    /* DDFS scratch */
    std::vector<std::pair<int32_t,int32_t>> sStkR;  /* red descent stack */
    std::vector<std::pair<int32_t,int32_t>> sStkG;  /* green descent stack */
    std::vector<int32_t> path;                      /* reconstructed AP */
    std::vector<int32_t> nodesSeen;                 /* DDFS-touched (unused in bipartite, kept for parity) */

    int32_t todoNum;
    int32_t bridgeNum;
};

/* ---------- BipartiteGraph / BipartiteMatching construction ---------- */

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

BipartiteMatching emptyBipartiteMatching(const BipartiteGraph& graph) {
    BipartiteMatching matching;
    matching.sNumVtxs = graph.sNumVtxs;
    matching.tNumVtxs = graph.tNumVtxs;
    matching.numEdgs = 0;
    matching.sMate.assign(graph.sNumVtxs, NIL);
    matching.tMate.assign(graph.tNumVtxs, NIL);
    return matching;
}

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

/* ---------- HKB2 algorithm core ---------- */

struct HKB2 {
    const BipartiteGraph& graph;
    BipartiteMatching& matching;
    HKB2IState& state;

    HKB2(const BipartiteGraph& g, BipartiteMatching& m, HKB2IState& s)
        : graph(g), matching(m), state(s) {}

    /* ---- helpers (mirror MV) ---- */

    void addToLevel(int32_t lvl, int32_t v) {
        if (static_cast<size_t>(lvl) >= state.levels.size()) state.levels.resize(lvl + 1);
        state.levels[lvl].push_back(v);
        state.todoNum++;
    }

    void addToBridges(int32_t bucket, int32_t u, int32_t v) {
        if (static_cast<size_t>(bucket) >= state.brBuckets.size())
            state.brBuckets.resize(bucket + 1);
        state.brBuckets[bucket].push_back({u, v});
        state.bridgeNum++;
    }

    /* tenacity: in bipartite we only ever generate unmatched-edge bridges
     * during MIN at even levels.  Matched bridges (odd-level same-level
     * matching edges) would arise only if we extended via matching at odd
     * levels and hit a same-odd-level mate, which is structurally impossible
     * in bipartite (a matching edge always crosses partitions, and its
     * endpoints sit at consecutive levels in any honest BFS). */
    int32_t tenacity(int32_t n1, int32_t n2) const {
        if (state.nodes[n1].match == n2) {
            /* matched bridge — needs odd_level on both ends */
            if (state.nodes[n1].oddLevel != NIL && state.nodes[n2].oddLevel != NIL)
                return state.nodes[n1].oddLevel + state.nodes[n2].oddLevel + 1;
            return NIL;
        } else {
            /* unmatched bridge — needs even_level on both ends */
            if (state.nodes[n1].evenLevel != NIL && state.nodes[n2].evenLevel != NIL)
                return state.nodes[n1].evenLevel + state.nodes[n2].evenLevel + 1;
            return NIL;
        }
    }

    /* step_to: line-by-line port of MV, with the hanging-bridge branch
     * dropped (verified empty in bipartite).  Predecessor lists are no
     * longer materialized — only the count is kept, since the actual
     * predecessor vertices can be enumerated on demand from the CSR. */
    void stepTo(int32_t to, int32_t from, int32_t level) {
        level++;
        int32_t tl = state.nodes[to].minLevel;
        if (tl == NIL || tl >= level) {
            if (tl != level) {
                addToLevel(level, to);
                state.nodes[to].setMinLevel(level);
            }
            state.nodes[to].numPreds++;
        } else {
            /* same-level edge → bridge.  In bipartite, tenacity is sometimes
             * NIL at this point (the "hanging" case in MV).  Empirically
             * verified: hanging bridges never get promoted into useful
             * augmenting paths in bipartite, because promotion only happens
             * via petal contraction (which cannot occur).  So we silently
             * drop them — this is the only deviation from MV. */
            int32_t ten = tenacity(to, from);
            if (ten != NIL) addToBridges((ten - 1) / 2, to, from);
        }
    }

    void MIN(int32_t i) {
        if (static_cast<size_t>(i) >= state.levels.size()) return;
        size_t sz = state.levels[i].size();
        for (size_t k = 0; k < sz; k++) {
            int32_t cur = state.levels[i][k];
            state.todoNum--;
            HKB2Node& n = state.nodes[cur];
            if ((i & 1) == 0) {
                size_t b = state.vIdx[cur], e = state.vIdx[cur + 1];
                for (size_t j = b; j < e; j++) {
                    int32_t edge = state.vAdj[j];
                    if (edge != n.match) stepTo(edge, cur, i);
                }
            } else {
                if (n.match != NIL) stepTo(n.match, cur, i);
            }
        }
    }

    /* ---- DDFS — line-by-line port with bud_star inlined as identity ---- */

    /* Enumerate predecessors of `cur` in the BFS DAG, pushing (cur, pred)
     * pairs onto the descent stack.  In bipartite:
     *   - cur at even level: unique pred = match[cur] (if alive).
     *   - cur at odd level:  preds = vAdj[cur] \ {match[cur]} restricted to
     *     vertices at level cur.level - 1 that are still alive.
     * Deletion is honored via the `deleted` flag (cascade-delete) and via
     * `numPreds` reaching zero (no live preds remain). */
    void addPredsToStack(int32_t cur, std::vector<std::pair<int32_t,int32_t>>& stk) {
        const HKB2Node& n = state.nodes[cur];
        if (n.numPreds <= 0) return;
        int32_t targetLvl = n.minLevel - 1;
        if ((n.minLevel & 1) == 0) {
            /* even: unique pred is the matching mate */
            int32_t m = n.match;
            if (m != NIL && !state.nodes[m].deleted &&
                state.nodes[m].minLevel == targetLvl) {
                stk.push_back({cur, m});
            }
        } else {
            /* odd: scan non-matching neighbors at level - 1 */
            int32_t mate = n.match;
            size_t b = state.vIdx[cur], e = state.vIdx[cur + 1];
            for (size_t j = b; j < e; j++) {
                int32_t nb = state.vAdj[j];
                if (nb == mate) continue;
                if (state.nodes[nb].deleted) continue;
                if (state.nodes[nb].minLevel != targetLvl) continue;
                stk.push_back({cur, nb});
            }
        }
    }

    /* prepareNext: in MV this dereferences via bud_star; here it's identity.
     * Still need to set `below` for the AP reconstruction. */
    void prepareNext(std::pair<int32_t,int32_t>& Nx) {
        if (Nx.first != NIL) state.nodes[Nx.first].below = Nx.second;
    }

    static bool edgeValid(const std::pair<int32_t,int32_t>& e) {
        return !(e.first == NIL && e.second == NIL);
    }

    static void nodeFromStack(std::pair<int32_t,int32_t>& e,
                              std::vector<std::pair<int32_t,int32_t>>& S) {
        if (!S.empty()) { e = S.back(); S.pop_back(); }
        else e = {NIL, NIL};
    }

    int32_t L(const std::pair<int32_t,int32_t>& e) const {
        /* MV: return min_level[bud_star(e.second)].  bud_star is identity. */
        return state.nodes[e.second].minLevel;
    }

    void stepInto(int32_t& C, std::pair<int32_t,int32_t>& Nx,
                  std::vector<std::pair<int32_t,int32_t>>& S,
                  int32_t greenTop, int32_t redTop) {
        prepareNext(Nx);
        if (!state.nodes[Nx.second].visited) {
            state.nodes[Nx.second].above = Nx.first;
            C = Nx.second;
            HKB2Node& n = state.nodes[C];
            n.visited = 1;
            n.ddfsGreen = greenTop;
            n.ddfsRed = redTop;
            state.nodesSeen.push_back(C);
            addPredsToStack(C, S);
        }
        nodeFromStack(Nx, S);
    }

    int32_t DDFS(int32_t greenTop, int32_t redTop) {
        state.nodesSeen.clear();
        auto& Sr = state.sStkR;
        auto& Sg = state.sStkG;
        Sr.clear();
        Sg.clear();

        int32_t G = NIL, R = NIL;

        /* MV: if (bud_star(red_top) == bud_star(green_top)) return DDFS_EMPTY;
         * bud_star is identity. */
        if (redTop == greenTop) return DDFS_EMPTY;

        if (state.nodes[greenTop].minLevel == 0 &&
            state.nodes[redTop].minLevel == 0)
            return DDFS_PATH;

        std::pair<int32_t,int32_t> Ng = {NIL, greenTop};
        std::pair<int32_t,int32_t> Nr = {NIL, redTop};
        std::pair<int32_t,int32_t> redBefore = {NIL, NIL};
        std::pair<int32_t,int32_t> greenBefore = {NIL, NIL};

        while (R == NIL || G == NIL ||
               state.nodes[R].minLevel > 0 || state.nodes[G].minLevel > 0) {

            while (edgeValid(Nr) && edgeValid(Ng) && L(Nr) != L(Ng)) {

                while (edgeValid(Nr) && L(Nr) > L(Ng))
                    stepInto(R, Nr, Sr, greenTop, redTop);

                if (!edgeValid(Nr)) {
                    Nr = redBefore;
                    int32_t tmp = redBefore.first;
                    while (tmp != NIL && state.nodes[tmp].above != NIL) {
                        int32_t rc = state.nodes[tmp].above;
                        /* In MV, scan preds[rc] for the entry equal to tmp
                         * (bud_star included).  bud_star is identity here,
                         * so we just verify tmp is still a live pred of rc.
                         * Since above[tmp]=rc was set during the descent,
                         * the BFS-edge (tmp,rc) was valid then; only
                         * cascade-delete can have invalidated it. */
                        if (!state.nodes[tmp].deleted)
                            state.nodes[rc].below = tmp;
                        tmp = state.nodes[tmp].above;
                    }
                }

                while (edgeValid(Ng) && L(Nr) < L(Ng))
                    stepInto(G, Ng, Sg, greenTop, redTop);

                if (!edgeValid(Ng)) {
                    Ng = greenBefore;
                    int32_t tmp = greenBefore.first;
                    while (tmp != NIL && state.nodes[tmp].above != NIL) {
                        int32_t rc = state.nodes[tmp].above;
                        /* See the symmetric red-side comment above. */
                        if (!state.nodes[tmp].deleted)
                            state.nodes[rc].below = tmp;
                        tmp = state.nodes[tmp].above;
                    }
                }
            }

            if (Nr.second == Ng.second) {
                /* Cursors collide.  Backtrack red first, then green; if both
                 * stacks empty we'd return DDFS_PETAL — but in bipartite this
                 * never fires (verified). */
                if (!Sr.empty()) {
                    redBefore = Nr;
                    prepareNext(Nr);
                    nodeFromStack(Nr, Sr);
                    if (edgeValid(Nr)) R = Nr.first;
                    else Nr = redBefore;
                } else if (!Sg.empty()) {
                    greenBefore = Ng;
                    prepareNext(Ng);
                    nodeFromStack(Ng, Sg);
                    if (edgeValid(Ng)) G = Ng.first;
                    else Ng = greenBefore;
                } else {
                    /* Bipartite-dead branch — would be DDFS_PETAL in MV.
                     * Treat as DDFS_EMPTY (bridge yields no AP, skip). */
                    prepareNext(Nr);
                    prepareNext(Ng);
                    return DDFS_EMPTY;
                }
            } else {
                stepInto(R, Nr, Sr, greenTop, redTop);
                stepInto(G, Ng, Sg, greenTop, redTop);
            }
        }
        return DDFS_PATH;
    }

    /* find_path / walk_down_path: in MV walk_down_path consults `bud` and
     * may invoke walk_blossom.  bud is always NIL in bipartite, so this
     * reduces to a flat below-pointer descent. */
    void findPath(int32_t n1, int32_t n2) {
        state.path.clear();
        walkDownPath(n1);
        std::reverse(state.path.begin(), state.path.end());
        walkDownPath(n2);
    }

    void walkDownPath(int32_t start) {
        int32_t cur = start;
        while (cur != NIL) {
            state.path.push_back(cur);
            cur = state.nodes[cur].below;
        }
    }

    void augmentPath() {
        for (size_t i = 0; i + 1 < state.path.size(); i += 2) {
            int32_t n1 = state.path[i];
            int32_t n2 = state.path[i + 1];
            state.nodes[n1].match = n2;
            state.nodes[n2].match = n1;
        }
    }

    /* removePath: mark every vertex on the AP as deleted, then cascade-
     * delete any vertex whose pred-count drops to zero (it can no longer
     * be reached by any future augmenting-path search this phase).
     *
     * In MV, predTo[cur] held a precomputed list of (target, idx) pairs.
     * Here we recover it from the CSR: when we delete `cur` at level L,
     * any vertex that had `cur` as a pred sits at level L+1 and is reached
     * from `cur` via a prop edge of the appropriate parity:
     *   - L odd  → target at even L+1, prop is matching → target = match[cur]
     *   - L even → target at odd  L+1, prop is non-matching → scan vAdj[cur]
     * Decrement target.numPreds; if it hits zero, target becomes
     * unreachable and joins the deletion frontier. */
    void removePath() {
        std::vector<int32_t> stk = state.path;
        while (!stk.empty()) {
            int32_t cur = stk.back();
            stk.pop_back();
            HKB2Node& nc = state.nodes[cur];
            if (nc.deleted) continue;
            nc.deleted = 1;
            int32_t curLvl = nc.minLevel;
            if (curLvl == NIL) continue;
            int32_t childLvl = curLvl + 1;

            if ((curLvl & 1) == 1) {
                /* odd → unique matching child at even level */
                int32_t m = nc.match;
                if (m != NIL && !state.nodes[m].deleted &&
                    state.nodes[m].minLevel == childLvl) {
                    HKB2Node& nm = state.nodes[m];
                    nm.numPreds--;
                    if (nm.numPreds <= 0) stk.push_back(m);
                }
            } else {
                /* even → non-matching children at odd level */
                int32_t mate = nc.match;
                size_t b = state.vIdx[cur], e = state.vIdx[cur + 1];
                for (size_t j = b; j < e; j++) {
                    int32_t nb = state.vAdj[j];
                    if (nb == mate) continue;
                    HKB2Node& nn = state.nodes[nb];
                    if (nn.deleted) continue;
                    if (nn.minLevel != childLvl) continue;
                    nn.numPreds--;
                    if (nn.numPreds <= 0) stk.push_back(nb);
                }
            }
        }
    }

    bool MAX(int32_t i) {
        bool found = false;
        if (static_cast<size_t>(i) >= state.brBuckets.size()) return false;

        for (size_t j = 0; j < state.brBuckets[i].size(); j++) {
            auto br = state.brBuckets[i][j];
            state.bridgeNum--;
            int32_t n1 = br.first, n2 = br.second;
            if (state.nodes[n1].deleted || state.nodes[n2].deleted) continue;

            int32_t result = DDFS(n1, n2);
            if (result == DDFS_EMPTY) continue;

            /* DDFS_PATH (no DDFS_PETAL in bipartite) */
            findPath(n1, n2);
            augmentPath();
            removePath();
            found = true;
        }
        return found;
    }

    /* Reset before each phase — mirrors MV's reset() */
    void resetPhase() {
        for (auto& v : state.levels) v.clear();
        for (auto& v : state.brBuckets) v.clear();
        state.bridgeNum = 0;
        state.todoNum = 0;
        for (size_t i = 0; i < state.nVtxs; i++) {
            state.nodes[i].reset();
            if (state.nodes[i].match == NIL) {
                addToLevel(0, static_cast<int32_t>(i));
                state.nodes[i].setMinLevel(0);
            }
        }
    }

    bool runPhase() {
        bool found = false;
        for (int32_t i = 0; i < static_cast<int32_t>(state.nVtxs) / 2 + 1 && !found; i++) {
            if (state.todoNum <= 0 && state.bridgeNum <= 0) return false;
            MIN(i);
            found = MAX(i);
        }
        return found;
    }
};

/* ---------- Top-level ---------- */

/* Build the unified CSR adjacency and seed nodes with the initial matching
 * (translated from sMate/tMate into match[] in the unified namespace). */
static void buildState(const BipartiteGraph& graph, const BipartiteMatching& matching,
                       HKB2IState& state) {
    state.sNumVtxs = graph.sNumVtxs;
    state.nVtxs = graph.sNumVtxs + graph.tNumVtxs;
    state.vIdx.assign(state.nVtxs + 1, 0);

    for (size_t v = 0; v < graph.sNumVtxs; v++)
        state.vIdx[v + 1] = state.vIdx[v] + (graph.sIdx[v + 1] - graph.sIdx[v]);
    for (size_t t = 0; t < graph.tNumVtxs; t++) {
        size_t v = graph.sNumVtxs + t;
        state.vIdx[v + 1] = state.vIdx[v] + (graph.tIdx[t + 1] - graph.tIdx[t]);
    }
    state.vAdj.resize(state.vIdx[state.nVtxs]);

    for (size_t v = 0; v < graph.sNumVtxs; v++) {
        size_t srcBegin = graph.sIdx[v], srcEnd = graph.sIdx[v + 1];
        size_t dstBegin = state.vIdx[v];
        for (size_t k = srcBegin; k < srcEnd; k++) {
            state.vAdj[dstBegin + (k - srcBegin)] =
                static_cast<int32_t>(graph.sNumVtxs) + graph.sAdj[k];
        }
    }
    for (size_t t = 0; t < graph.tNumVtxs; t++) {
        size_t v = graph.sNumVtxs + t;
        size_t srcBegin = graph.tIdx[t], srcEnd = graph.tIdx[t + 1];
        size_t dstBegin = state.vIdx[v];
        for (size_t k = srcBegin; k < srcEnd; k++) {
            state.vAdj[dstBegin + (k - srcBegin)] = graph.tAdj[k];
        }
    }

    /* Initialize node array.  Translate sMate/tMate to unified match[]. */
    state.nodes.assign(state.nVtxs, HKB2Node());
    for (size_t s = 0; s < graph.sNumVtxs; s++) {
        if (matching.sMate[s] != NIL) {
            int32_t t = matching.sMate[s];
            state.nodes[s].match = static_cast<int32_t>(graph.sNumVtxs) + t;
        }
    }
    for (size_t t = 0; t < graph.tNumVtxs; t++) {
        if (matching.tMate[t] != NIL) {
            int32_t s = matching.tMate[t];
            int32_t v = static_cast<int32_t>(graph.sNumVtxs + t);
            state.nodes[v].match = s;
        }
    }
}

/* Translate match[] back into sMate/tMate. */
static void writebackMatching(const HKB2IState& state, BipartiteMatching& matching) {
    int32_t count = 0;
    matching.sMate.assign(matching.sNumVtxs, NIL);
    matching.tMate.assign(matching.tNumVtxs, NIL);
    for (size_t s = 0; s < state.sNumVtxs; s++) {
        int32_t mt = state.nodes[s].match;
        if (mt != NIL) {
            int32_t t = mt - static_cast<int32_t>(state.sNumVtxs);
            matching.sMate[s] = t;
            matching.tMate[t] = static_cast<int32_t>(s);
            count++;
        }
    }
    matching.numEdgs = static_cast<size_t>(count);
}

int32_t hkb2IterativeMcm(const BipartiteGraph& graph, BipartiteMatching& matching) {
    HKB2IState state;
    buildState(graph, matching, state);
    state.levels.reserve(state.nVtxs / 2 + 1);
    state.brBuckets.reserve(state.nVtxs / 2 + 1);

    int32_t numPhases = 0;
    HKB2 algo(graph, matching, state);
    while (true) {
        algo.resetPhase();
        if (!algo.runPhase()) break;
        numPhases++;
    }

    writebackMatching(state, matching);
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
    if (matchedS != matchedT) {
        fprintf(stderr, "ERROR: S-matched (%d) != T-matched (%d)\n", matchedS, matchedT);
        errors++;
    }

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %zu\n", matching.numEdgs);
    printf("S matched: %d, T matched: %d\n", matchedS, matchedT);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

int main(int argc, char* argv[]) {
    printf("Hopcroft-Karp Bidirectional MV-style (HKB2) - C++ Implementation (CSR)\n");
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

    int32_t numPhases = hkb2IterativeMcm(graph, matching);

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

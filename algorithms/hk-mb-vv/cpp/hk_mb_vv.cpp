/*
 * Hopcroft-Karp Pure Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * VV (vector-of-vectors) adjacency.
 * Matchbox-style implementation:
 *   - Iterative stack-based DFS (no recursion)
 *   - Edge index array for O(E) per-phase DFS
 *   - Selective cleanup (reset only visited vertices)
 *   - Circular queue containers
 *   - Bidirectional search (from smaller partition)
 *   - Lookahead variant (optional, separate top-level function)
 *   - Greedy and greedy-md initialization
 *   - Comprehensive per-phase statistics (#ifdef STATS)
 *
 * Refactored to separate input (BipartiteGraph), output (BipartiteMatching),
 * and algorithm state (HKPState).
 *
 * State-side note on direction:
 *   When rvrs == true, HKPState's "s" fields refer to the input graph's
 *   t-vertices (we search starting from t). The output BipartiteMatching
 *   is always in original (s, t) coordinates regardless of rvrs.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <climits>
#include <cassert>

constexpr int32_t NIL = -1;
constexpr int32_t INF_LEVEL = INT_MAX;

/* ---------- Vertex status during BFS/DFS ---------- */

enum Stt : unsigned char {
    STT_IDLE = 0,
    STT_BFS_QUEUED,
    STT_BFS_DONE,
    STT_DFS_ACTIVE,
    STT_DFS_DONE,
    STT_LAST
};

/* ---------- Helper container types ---------- */

struct CircQueue {
    std::vector<int32_t> buf;
    int32_t cap, sz, head, tail;

    CircQueue() : cap(0), sz(0), head(0), tail(0) {}
    void init(int32_t c) { cap = c; sz = head = tail = 0; buf.assign(c, 0); }
    bool empty() const { return sz == 0; }
    void push(int32_t v) { assert(sz < cap); buf[tail] = v; tail = (tail + 1 < cap) ? tail + 1 : 0; ++sz; }
    int32_t front() const { assert(sz > 0); return buf[head]; }
    void pop() { assert(sz > 0); head = (head + 1 < cap) ? head + 1 : 0; --sz; }
    void clear() { sz = head = tail = 0; }
};

struct Stack {
    std::vector<int32_t> buf;
    int32_t sz;

    Stack() : sz(0) {}
    void init(int32_t c) { buf.resize(c); sz = 0; }
    bool empty() const { return sz == 0; }
    void push(int32_t v) { buf[sz++] = v; }
    int32_t top() const { return buf[sz - 1]; }
    void pop() { --sz; }
    void clear() { sz = 0; }
};

struct IdxQueue {
    std::vector<int32_t> nxt, prv;
    int32_t head, sz, cap;

    IdxQueue() : head(NIL), sz(0), cap(0) {}
    void init(int32_t c) { cap = c; sz = 0; head = NIL; nxt.assign(c, NIL); prv.assign(c, NIL); }
    bool empty() const { return sz == 0; }
    void push(int32_t v) {
        nxt[v] = head;
        prv[v] = NIL;
        if (head != NIL) prv[head] = v;
        head = v;
        ++sz;
    }
    void erase(int32_t v) {
        if (prv[v] != NIL) nxt[prv[v]] = nxt[v]; else head = nxt[v];
        if (nxt[v] != NIL) prv[nxt[v]] = prv[v];
        nxt[v] = prv[v] = NIL;
        --sz;
    }
    int32_t first() const { return head; }
    int32_t next(int32_t v) const { return nxt[v]; }
};

/* ---------- Statistics (compiled only with -DSTATS) ---------- */

#ifdef STATS
struct PhaseStats {
    long long bfs_vtx, bfs_edg, dfs_vtx, dfs_edg;
    int num_augmentations;
    int shortest_path_len;
    long long agg_aug_path_len;
    int min_aug_path_len, max_aug_path_len;
};

struct AlgoStats {
    std::vector<PhaseStats> phases;
    int greedy_card;
    bool reversed;

    AlgoStats() { reset(); }
    void reset() { phases.clear(); greedy_card = 0; reversed = false; }
    void print() const {
        long long tb_v = 0, tb_e = 0, td_v = 0, td_e = 0;
        int t_aug = 0; long long t_apl = 0;
        for (auto& p : phases) {
            tb_v += p.bfs_vtx; tb_e += p.bfs_edg;
            td_v += p.dfs_vtx; td_e += p.dfs_edg;
            t_aug += p.num_augmentations;
            t_apl += p.agg_aug_path_len;
        }
        printf("\n=== Algorithm Statistics ===\n");
        printf("Search direction: %s\n", reversed ? "T->S (reversed)" : "S->T (normal)");
        printf("Greedy init card: %d\n", greedy_card);
        printf("Phases: %d\n", (int)phases.size());
        printf("Total BFS vertices: %lld, edges: %lld\n", tb_v, tb_e);
        printf("Total DFS vertices: %lld, edges: %lld\n", td_v, td_e);
        printf("Total augmentations: %d\n", t_aug);
        printf("Total aug path length: %lld\n", t_apl);
        if (t_aug > 0) printf("Avg aug path length: %.2f\n", (double)t_apl / t_aug);
        printf("\nPer-phase breakdown:\n");
        for (int i = 0; i < (int)phases.size(); i++) {
            const auto& p = phases[i];
            printf("  Phase %d: shortest=%d augs=%d bfs_v=%lld bfs_e=%lld "
                   "dfs_v=%lld dfs_e=%lld avg_len=%.1f\n",
                   i + 1, p.shortest_path_len, p.num_augmentations,
                   p.bfs_vtx, p.bfs_edg, p.dfs_vtx, p.dfs_edg,
                   p.num_augmentations > 0 ?
                   (double)p.agg_aug_path_len / p.num_augmentations : 0.0);
        }
        printf("===========================\n\n");
    }
};
#endif

/* ---------- Input: BipartiteGraph ---------- */

struct BipartiteGraph {
    size_t sNumVtxs;
    size_t tNumVtxs;
    size_t numEdgs;
    std::vector<std::vector<int32_t>> sAdj;
    std::vector<std::vector<int32_t>> tAdj;
};

/* ---------- Output: BipartiteMatching ---------- */

struct BipartiteMatching {
    size_t sNumVtxs;
    size_t tNumVtxs;
    size_t numEdgs;
    std::vector<int32_t> sMate;
    std::vector<int32_t> tMate;
};

/* ---------- State: HKPState ----------
 *
 * "src" semantics: the side we are searching from. When rvrs=false, src=s
 * (input s-side). When rvrs=true, src=t (input t-side).
 *
 * The fields named sX correspond to the algorithm's source side; tX to the
 * algorithm's target side. They are NOT the same as the input graph's
 * sAdj/tAdj when rvrs=true.
 */

struct HKPState {
    bool rvrs;
    size_t sCount;       // size of the search-source side (depends on rvrs)
    size_t tCount;       // size of the search-target side (depends on rvrs)

    std::vector<int32_t> sLevel, tLevel;
    std::vector<Stt>     sStt,   tStt;
    std::vector<int32_t> sPtr;
    std::vector<int32_t> sEdgeIdx;
    std::vector<int32_t> sLkhdIdx;

    CircQueue sBfsQueue;
    CircQueue sDoneQueue;
    CircQueue tDoneQueue;
    Stack     sDfsStack;
    CircQueue sLastQueue;
    CircQueue tLastQueue;
    IdxQueue  sExposed;

#ifdef STATS
    AlgoStats stats;
#endif
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

/* ---------- HKPState construction ---------- */

void initHKPState(HKPState& s, const BipartiteGraph& g, bool rvrs) {
    s.rvrs = rvrs;
    s.sCount = rvrs ? g.tNumVtxs : g.sNumVtxs;
    s.tCount = rvrs ? g.sNumVtxs : g.tNumVtxs;

    s.sLevel.assign(s.sCount, INF_LEVEL);
    s.tLevel.assign(s.tCount, INF_LEVEL);
    s.sStt.assign(s.sCount, STT_IDLE);
    s.tStt.assign(s.tCount, STT_IDLE);
    s.sPtr.assign(s.sCount, NIL);
    s.sEdgeIdx.assign(s.sCount, 0);
    s.sLkhdIdx.assign(s.sCount, 0);

    int32_t sc = std::max((int32_t)s.sCount, 1);
    int32_t tc = std::max((int32_t)s.tCount, 1);
    s.sBfsQueue.init(sc);
    s.sDoneQueue.init(sc);
    s.tDoneQueue.init(tc);
    s.sDfsStack.init(sc);
    int32_t mp = std::max(std::min((int32_t)s.sCount, (int32_t)s.tCount), 1);
    s.sLastQueue.init(mp);
    s.tLastQueue.init(mp);
    s.sExposed.init((int32_t)s.sCount);

#ifdef STATS
    s.stats.reset();
    s.stats.reversed = rvrs;
#endif
}

/* ---------- Adjacency accessors honoring rvrs ---------- *
 *
 * These return the adjacency arrays in the algorithm's "search source" /
 * "search target" frame.
 */

static inline const std::vector<int32_t>& srcAdj(const BipartiteGraph& g, const HKPState& s, int32_t u) {
    return s.rvrs ? g.tAdj[u] : g.sAdj[u];
}

/* ---------- Mate accessors honoring rvrs ---------- */

static inline int32_t srcMate(const BipartiteMatching& m, const HKPState& s, int32_t u) {
    return s.rvrs ? m.tMate[u] : m.sMate[u];
}
static inline int32_t tgtMate(const BipartiteMatching& m, const HKPState& s, int32_t v) {
    return s.rvrs ? m.sMate[v] : m.tMate[v];
}
static inline void setMate(BipartiteMatching& m, const HKPState& s, int32_t u, int32_t v) {
    if (s.rvrs) { m.tMate[u] = v; m.sMate[v] = u; }
    else        { m.sMate[u] = v; m.tMate[v] = u; }
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

/* ---------- BFS: compute shortest augmenting path length ---------- */

static int32_t bfsShortestLevel(const BipartiteGraph& g, const BipartiteMatching& m, HKPState& s
#ifdef STATS
    , PhaseStats& ps
#endif
) {
    int32_t shortest = INF_LEVEL;

    int32_t sf = s.sExposed.first();
    while (sf != NIL) {
    #ifdef STATS
        ps.bfs_vtx++;
    #endif
        s.sLevel[sf] = 0;
        s.sBfsQueue.push(sf);
        s.sStt[sf] = STT_BFS_QUEUED;
        sf = s.sExposed.next(sf);
    }

    int32_t curLayer = 0;
    while (!s.sBfsQueue.empty()) {
        int32_t src = s.sBfsQueue.front();
        s.sBfsQueue.pop();
        s.sDoneQueue.push(src);
        s.sStt[src] = STT_BFS_DONE;

        if (curLayer < s.sLevel[src] / 2) {
            if (shortest != INF_LEVEL) break;
            curLayer++;
        }

        const std::vector<int32_t>& adj = srcAdj(g, s, src);
        int32_t numEdg = (int32_t)adj.size();
        for (int32_t i = 0; i < numEdg; i++) {
            int32_t t = adj[i];
        #ifdef STATS
            ps.bfs_edg++;
        #endif
            if (s.tStt[t] == STT_BFS_DONE) continue;
        #ifdef STATS
            ps.bfs_vtx++;
        #endif
            s.tLevel[t] = s.sLevel[src] + 1;
            s.tDoneQueue.push(t);
            s.tStt[t] = STT_BFS_DONE;

            int32_t ss = tgtMate(m, s, t);
            if (ss == NIL) {
                shortest = s.tLevel[t];
            } else {
            #ifdef STATS
                ps.bfs_edg++;
                ps.bfs_vtx++;
            #endif
                s.sLevel[ss] = s.tLevel[t] + 1;
                s.sBfsQueue.push(ss);
                s.sStt[ss] = STT_BFS_QUEUED;
            }
        }
    }
    return shortest;
}

/* ---------- DFS: HK-style (find maximal set of shortest augmenting paths) ---------- */

static void dfsFindPathsHk(const BipartiteGraph& g, const BipartiteMatching& m, HKPState& s,
                           int32_t shortest
#ifdef STATS
    , PhaseStats& ps
#endif
) {
    int32_t sf = s.sExposed.first();
    while (sf != NIL) {
        int32_t sfNext = s.sExposed.next(sf);
    #ifdef STATS
        ps.dfs_vtx++;
    #endif
        s.sDfsStack.push(sf);
        s.sStt[sf] = STT_DFS_ACTIVE;

        while (!s.sDfsStack.empty()) {
            int32_t src = s.sDfsStack.top();
            const std::vector<int32_t>& adj = srcAdj(g, s, src);
            int32_t numEdg = (int32_t)adj.size();

            while (s.sEdgeIdx[src] < numEdg) {
                int32_t t = adj[s.sEdgeIdx[src]];
            #ifdef STATS
                ps.dfs_edg++;
            #endif
                int32_t ss = tgtMate(m, s, t);
                if ((s.tLevel[t] != s.sLevel[src] + 1) ||
                    ((ss == NIL) && (s.tStt[t] == STT_LAST)) ||
                    ((ss != NIL) && ((s.sLevel[src] + 1 == shortest) ||
                                     (s.sStt[ss] == STT_DFS_ACTIVE) ||
                                     (s.sStt[ss] == STT_DFS_DONE)))) {
                    s.sEdgeIdx[src]++;
                } else {
                    break;
                }
            }

            if (s.sEdgeIdx[src] < numEdg) {
                int32_t t = adj[s.sEdgeIdx[src]];
                int32_t ss = tgtMate(m, s, t);
            #ifdef STATS
                ps.dfs_vtx++;
            #endif
                if (ss == NIL) {
                    s.sLastQueue.push(src);
                    s.tLastQueue.push(t);
                    s.tStt[t] = STT_LAST;
                    while (!s.sDfsStack.empty()) {
                        s.sStt[s.sDfsStack.top()] = STT_DFS_DONE;
                        s.sDfsStack.pop();
                    }
                    break;
                }
            #ifdef STATS
                ps.dfs_edg++;
                ps.dfs_vtx++;
            #endif
                s.sPtr[ss] = src;
                s.sDfsStack.push(ss);
                s.sStt[ss] = STT_DFS_ACTIVE;
                s.sEdgeIdx[src]++;
            } else {
                s.sDfsStack.pop();
                s.sStt[src] = STT_DFS_DONE;
            }
        }
        sf = sfNext;
    }
}

/* ---------- DFS: lookahead-style (no BFS layering, matchbox-style) ---------- */

static void dfsFindPathsLkhd(const BipartiteGraph& g, const BipartiteMatching& m, HKPState& s
#ifdef STATS
    , PhaseStats& ps
#endif
) {
    int32_t sf = s.sExposed.first();
    while (sf != NIL) {
        int32_t sfNext = s.sExposed.next(sf);
    #ifdef STATS
        ps.dfs_vtx++;
    #endif
        s.sDfsStack.push(sf);
        s.sStt[sf] = STT_DFS_ACTIVE;

        while (!s.sDfsStack.empty()) {
            int32_t src = s.sDfsStack.top();
            const std::vector<int32_t>& adj = srcAdj(g, s, src);
            int32_t numEdg = (int32_t)adj.size();

            /* Lookahead: check for direct free t-neighbor */
            while (s.sLkhdIdx[src] < numEdg) {
                int32_t t = adj[s.sLkhdIdx[src]];
            #ifdef STATS
                ps.bfs_edg++;
            #endif
                int32_t ss = tgtMate(m, s, t);
                if ((ss == NIL && s.tStt[t] != STT_IDLE) || ss != NIL) {
                    s.sLkhdIdx[src]++;
                } else {
                    break;
                }
            }

            if (s.sLkhdIdx[src] < numEdg) {
                int32_t t = adj[s.sLkhdIdx[src]];
            #ifdef STATS
                ps.bfs_vtx++;
            #endif
                s.sLastQueue.push(src);
                s.tLastQueue.push(t);
                s.tStt[t] = STT_LAST;
                while (!s.sDfsStack.empty()) {
                    int32_t sss = s.sDfsStack.top();
                    s.sDfsStack.pop();
                    s.sDoneQueue.push(sss);
                    s.sStt[sss] = STT_DFS_DONE;
                }
                s.sLkhdIdx[src]++;
                break;
            }

            /* Regular DFS deeper */
            while (s.sEdgeIdx[src] < numEdg) {
                int32_t t = adj[s.sEdgeIdx[src]];
            #ifdef STATS
                ps.dfs_edg++;
            #endif
                int32_t ss = tgtMate(m, s, t);
                if (ss == NIL || (ss != NIL && s.sStt[ss] != STT_IDLE)) {
                    s.sEdgeIdx[src]++;
                } else {
                    break;
                }
            }

            if (s.sEdgeIdx[src] < numEdg) {
                int32_t t = adj[s.sEdgeIdx[src]];
                int32_t ss = tgtMate(m, s, t);
            #ifdef STATS
                ps.dfs_vtx++;
                ps.dfs_edg++;
                ps.dfs_vtx++;
            #endif
                s.sPtr[ss] = src;
                s.sDfsStack.push(ss);
                s.sStt[ss] = STT_DFS_ACTIVE;
                s.sEdgeIdx[src]++;
            } else {
                s.sDfsStack.pop();
                s.sDoneQueue.push(src);
                s.sStt[src] = STT_DFS_DONE;
            }
        }
        sf = sfNext;
    }
}

/* ---------- Augment along path from sLast to root via sPtr ---------- */

static int32_t augment(BipartiteMatching& m, HKPState& s, int32_t sLast, int32_t tLast) {
    int32_t src = sLast, t = tLast, k = 0;
    while (src != NIL) {
        int32_t tt = srcMate(m, s, src);
        setMate(m, s, src, t);
        if (s.sPtr[src] == NIL) s.sExposed.erase(src);
        src = s.sPtr[src]; t = tt; k++;
    }
    return 2 * k - 1;
}

/* ---------- Cleanup after HK phase ---------- */

static void cleanupHk(HKPState& s) {
    while (!s.sBfsQueue.empty()) {
        int32_t u = s.sBfsQueue.front(); s.sBfsQueue.pop();
        s.sPtr[u] = NIL; s.sLevel[u] = INF_LEVEL;
        s.sStt[u] = STT_IDLE; s.sEdgeIdx[u] = 0; s.sLkhdIdx[u] = 0;
    }
    while (!s.sDoneQueue.empty()) {
        int32_t u = s.sDoneQueue.front(); s.sDoneQueue.pop();
        s.sPtr[u] = NIL; s.sLevel[u] = INF_LEVEL;
        s.sStt[u] = STT_IDLE; s.sEdgeIdx[u] = 0; s.sLkhdIdx[u] = 0;
    }
    while (!s.tDoneQueue.empty()) {
        int32_t t = s.tDoneQueue.front(); s.tDoneQueue.pop();
        s.tLevel[t] = INF_LEVEL; s.tStt[t] = STT_IDLE;
    }
}

/* ---------- Cleanup after lookahead phase ---------- */

static void cleanupLkhd(HKPState& s) {
    while (!s.sDoneQueue.empty()) {
        int32_t u = s.sDoneQueue.front(); s.sDoneQueue.pop();
        s.sPtr[u] = NIL; s.sStt[u] = STT_IDLE;
        s.sEdgeIdx[u] = 0; s.sLkhdIdx[u] = 0;
    }
}

/* ---------- Top-level Hopcroft-Karp Pure (HK mode) ---------- */

void hopcroftKarpPure(const BipartiteGraph& g, BipartiteMatching& m, bool rvrs) {
    HKPState s;
    initHKPState(s, g, rvrs);

#ifdef STATS
    s.stats.greedy_card = (int)m.numEdgs;
#endif

    /* Push exposed src vertices into the IdxQueue */
    for (size_t u = 0; u < s.sCount; u++) {
        if (srcMate(m, s, (int32_t)u) == NIL) s.sExposed.push((int32_t)u);
    }

    int32_t phaseCount = 0;
    int32_t newEdgs = 0;
    while (true) {
    #ifdef STATS
        PhaseStats ps = {0, 0, 0, 0, 0, 0, 0, INF_LEVEL, 0};
    #endif

        int32_t shortest = bfsShortestLevel(g, m, s
        #ifdef STATS
            , ps
        #endif
        );

        if (shortest == INF_LEVEL) {
            cleanupHk(s);
            break;
        }

        phaseCount++;
    #ifdef STATS
        ps.shortest_path_len = shortest;
    #endif

        dfsFindPathsHk(g, m, s, shortest
        #ifdef STATS
            , ps
        #endif
        );

        while (!s.tLastQueue.empty()) {
            int32_t sl = s.sLastQueue.front(); s.sLastQueue.pop();
            int32_t tl = s.tLastQueue.front(); s.tLastQueue.pop();
            int32_t pathLen = augment(m, s, sl, tl);
            newEdgs++;
        #ifdef STATS
            ps.num_augmentations++;
            ps.agg_aug_path_len += pathLen;
            if (pathLen < ps.min_aug_path_len) ps.min_aug_path_len = pathLen;
            if (pathLen > ps.max_aug_path_len) ps.max_aug_path_len = pathLen;
        #else
            (void)pathLen;
        #endif
        }

        cleanupHk(s);

    #ifdef STATS
        s.stats.phases.push_back(ps);
    #endif
    }

    m.numEdgs += newEdgs;

#ifdef STATS
    s.stats.print();
#endif
    printf("Phases: %d\n", phaseCount);
}

/* ---------- Top-level Hopcroft-Karp Pure (lookahead mode) ---------- */

void hopcroftKarpPureLkhd(const BipartiteGraph& g, BipartiteMatching& m, bool rvrs) {
    HKPState s;
    initHKPState(s, g, rvrs);

#ifdef STATS
    s.stats.greedy_card = (int)m.numEdgs;
#endif

    for (size_t u = 0; u < s.sCount; u++) {
        if (srcMate(m, s, (int32_t)u) == NIL) s.sExposed.push((int32_t)u);
    }

    int32_t passCount = 0;
    int32_t newEdgs = 0;
    while (true) {
    #ifdef STATS
        PhaseStats ps = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    #endif

        dfsFindPathsLkhd(g, m, s
        #ifdef STATS
            , ps
        #endif
        );

        if (s.tLastQueue.empty()) break;

        passCount++;
        while (!s.tLastQueue.empty()) {
            int32_t sl = s.sLastQueue.front(); s.sLastQueue.pop();
            int32_t tl = s.tLastQueue.front(); s.tLastQueue.pop();
            s.tStt[tl] = STT_IDLE;
            int32_t pathLen = augment(m, s, sl, tl);
            newEdgs++;
        #ifdef STATS
            ps.num_augmentations++;
            ps.agg_aug_path_len += pathLen;
        #else
            (void)pathLen;
        #endif
        }

        cleanupLkhd(s);

    #ifdef STATS
        s.stats.phases.push_back(ps);
    #endif
    }

    m.numEdgs += newEdgs;

#ifdef STATS
    s.stats.print();
#endif
    printf("Phases: %d\n", passCount);
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
    printf("Hopcroft-Karp Pure Algorithm - C++ Implementation (VV)\n");
    printf("========================================================\n\n");

    if (argc < 2) {
        printf("Usage: %s <filename> [--greedy|--greedy-md] [--lkhd] [--rvrs|--fwd]\n", argv[0]);
        return 1;
    }

    int greedyMode = 0;
    bool useLkhd = false;
    bool forceRvrs = false;
    bool autoDir = true;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedyMode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedyMode = 2;
        else if (std::string(argv[i]) == "--lkhd") useLkhd = true;
        else if (std::string(argv[i]) == "--rvrs") { forceRvrs = true; autoDir = false; }
        else if (std::string(argv[i]) == "--fwd")  { forceRvrs = false; autoDir = false; }
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

    bool useRvrs = forceRvrs;
    if (autoDir) useRvrs = (tNumVtxs < sNumVtxs);

    printf("Mode: %s\n", useLkhd ? "lookahead DFS" : "HK (BFS+DFS)");
    printf("Direction: %s%s\n", useRvrs ? "T->S (reversed)" : "S->T (normal)",
           autoDir ? " (auto)" : " (forced)");

    BipartiteGraph bipartiteGraph = buildBipartiteGraph(sNumVtxs, tNumVtxs, edges);
    BipartiteMatching bipartiteMatching = emptyBipartiteMatching(bipartiteGraph);

    auto t0 = std::chrono::high_resolution_clock::now();

    int32_t greedySize = 0;
    if (greedyMode == 1) greedySize = greedyInit(bipartiteGraph, bipartiteMatching);
    else if (greedyMode == 2) greedySize = greedyInitMd(bipartiteGraph, bipartiteMatching);

    if (useLkhd) hopcroftKarpPureLkhd(bipartiteGraph, bipartiteMatching, useRvrs);
    else         hopcroftKarpPure(bipartiteGraph, bipartiteMatching, useRvrs);

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

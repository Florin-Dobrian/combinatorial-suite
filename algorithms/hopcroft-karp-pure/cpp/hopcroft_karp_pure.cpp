/*
 * Hopcroft-Karp Pure Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * Incorporates all optimizations from the Matchbox library:
 * - Iterative stack-based DFS (no recursion)
 * - Edge index array (sIdxArr) for O(E) per-phase DFS
 * - Selective cleanup (reset only visited vertices)
 * - Circular queue containers (no overflow on reuse)
 * - Bidirectional search (from smaller partition)
 * - Lookahead variant (eMplPathDfsLkhd style, optional)
 * - Greedy and greedy-md initialization
 * - Comprehensive per-phase statistics (#ifdef STATS)
 *
 * Input format: "L R M" header then M lines of "u v" edges (0-indexed).
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <climits>
#include <cassert>

/* =========================================================================
 * Constants
 * ========================================================================= */
static const int NIL = -1;
static const int INF_LEVEL = INT_MAX;

/* Vertex status during BFS/DFS */
enum Stt : unsigned char {
    STT_IDLE = 0,
    STT_BFS_QUEUED,   /* in BFS queue, not yet processed */
    STT_BFS_DONE,     /* BFS processed */
    STT_DFS_ACTIVE,   /* on DFS stack */
    STT_DFS_DONE,     /* DFS finished (found path or dead end) */
    STT_LAST          /* T-vertex: endpoint of found augmenting path */
};

/* =========================================================================
 * Circular Queue (matches Matchbox VecItmQue)
 *
 * Fixed-capacity circular buffer. Capacity = max simultaneous items.
 * Push/Pop wrap around, so total pushes can exceed capacity safely.
 * ========================================================================= */
struct CircQueue {
    std::vector<int> buf;
    int cap, sz, head, tail;

    CircQueue() : cap(0), sz(0), head(0), tail(0) {}
    void init(int c) {
        cap = c; sz = head = tail = 0;
        buf.assign(c, 0);
    }
    bool empty() const { return sz == 0; }
    void push(int v) {
        assert(sz < cap);
        buf[tail] = v;
        tail = (tail + 1 < cap) ? tail + 1 : 0;
        ++sz;
    }
    int front() const { assert(sz > 0); return buf[head]; }
    void pop() {
        assert(sz > 0);
        head = (head + 1 < cap) ? head + 1 : 0;
        --sz;
    }
    void clear() { sz = head = tail = 0; }
};

/* =========================================================================
 * Stack (matches Matchbox VecItmStk)
 * ========================================================================= */
struct Stack {
    std::vector<int> buf;
    int sz;

    Stack() : sz(0) {}
    void init(int c) { buf.resize(c); sz = 0; }
    bool empty() const { return sz == 0; }
    void push(int v) { buf[sz++] = v; }
    int top() const { return buf[sz - 1]; }
    void pop() { --sz; }
    void clear() { sz = 0; }
};

/* =========================================================================
 * Indexed Queue (matches Matchbox LstItmIdxdQue but array-based)
 *
 * Supports O(1) Push, Erase, First, Next for the exposed-vertex set.
 * Doubly-linked list embedded in arrays.
 * ========================================================================= */
struct IdxQueue {
    std::vector<int> nxt, prv;
    int head, sz, cap;

    IdxQueue() : head(NIL), sz(0), cap(0) {}
    void init(int c) {
        cap = c; sz = 0; head = NIL;
        nxt.assign(c, NIL);
        prv.assign(c, NIL);
    }
    bool empty() const { return sz == 0; }
    void push(int v) {
        nxt[v] = head;
        prv[v] = NIL;
        if (head != NIL) prv[head] = v;
        head = v;
        ++sz;
    }
    void erase(int v) {
        if (prv[v] != NIL) nxt[prv[v]] = nxt[v];
        else head = nxt[v];
        if (nxt[v] != NIL) prv[nxt[v]] = prv[v];
        nxt[v] = prv[v] = NIL;
        --sz;
    }
    int first() const { return head; }
    int next(int v) const { return nxt[v]; }
};

/* =========================================================================
 * Statistics (compiled only with -DSTATS)
 * ========================================================================= */
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
        if (t_aug > 0)
            printf("Avg aug path length: %.2f\n", (double)t_apl / t_aug);
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

/* =========================================================================
 * HopcroftKarpPure
 * ========================================================================= */
struct HopcroftKarpPure {
    int s_count, t_count;
    std::vector<std::vector<int>> src_adj;

    std::vector<int> s_mate, t_mate;
    std::vector<int> s_level, t_level;
    std::vector<Stt> s_stt, t_stt;
    std::vector<int> s_ptr;
    std::vector<int> s_idx;
    std::vector<int> s_lkhd;

    CircQueue s_bfs_queue;
    CircQueue s_done_queue;
    CircQueue t_done_queue;
    Stack     s_dfs_stack;
    CircQueue s_last_queue;
    CircQueue t_last_queue;
    IdxQueue  s_exposed;

    int greedy_size;
    bool rvrs;

#ifdef STATS
    AlgoStats stats;
#endif

    HopcroftKarpPure() : s_count(0), t_count(0), greedy_size(0), rvrs(false) {}

    void build(int left_count, int right_count,
               const std::vector<std::pair<int,int>>& edges, bool force_rvrs) {
        rvrs = force_rvrs;
        if (rvrs) { s_count = right_count; t_count = left_count; }
        else      { s_count = left_count;  t_count = right_count; }

        src_adj.resize(s_count);
        for (auto& e : edges) {
            int u = e.first, v = e.second;
            if (rvrs) {
                if (v >= 0 && v < s_count && u >= 0 && u < t_count)
                    src_adj[v].push_back(u);
            } else {
                if (u >= 0 && u < s_count && v >= 0 && v < t_count)
                    src_adj[u].push_back(v);
            }
        }
        for (int i = 0; i < s_count; i++) {
            std::sort(src_adj[i].begin(), src_adj[i].end());
            src_adj[i].erase(std::unique(src_adj[i].begin(), src_adj[i].end()),
                             src_adj[i].end());
        }

        s_mate.assign(s_count, NIL);
        t_mate.assign(t_count, NIL);
        s_level.assign(s_count, INF_LEVEL);
        t_level.assign(t_count, INF_LEVEL);
        s_stt.assign(s_count, STT_IDLE);
        t_stt.assign(t_count, STT_IDLE);
        s_ptr.assign(s_count, NIL);
        s_idx.assign(s_count, 0);
        s_lkhd.assign(s_count, 0);

        int sc = std::max(s_count, 1);
        int tc = std::max(t_count, 1);
        s_bfs_queue.init(sc);
        s_done_queue.init(sc);
        t_done_queue.init(tc);
        s_dfs_stack.init(sc);
        int mp = std::max(std::min(s_count, t_count), 1);
        s_last_queue.init(mp);
        t_last_queue.init(mp);
        s_exposed.init(s_count);
        greedy_size = 0;
    }

    /* Greedy: simple sequential */
    int greedy_init() {
        int cnt = 0;
        for (int s = 0; s < s_count; s++) {
            if (s_mate[s] != NIL) continue;
            for (int t : src_adj[s]) {
                if (t_mate[t] == NIL) {
                    s_mate[s] = t; t_mate[t] = s; cnt++; break;
                }
            }
        }
        return cnt;
    }

    /* Greedy: min-degree heuristic */
    int greedy_init_md() {
        int cnt = 0;
        std::vector<int> t_deg(t_count, 0);
        for (int s = 0; s < s_count; s++)
            for (int t : src_adj[s]) t_deg[t]++;
        std::vector<int> order(s_count);
        for (int i = 0; i < s_count; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            int da = (int)src_adj[a].size(), db = (int)src_adj[b].size();
            return da < db || (da == db && a < b);
        });
        for (int s : order) {
            if (s_mate[s] != NIL) continue;
            int best = NIL, best_deg = INT_MAX;
            for (int t : src_adj[s]) {
                if (t_mate[t] == NIL && t_deg[t] < best_deg) {
                    best = t; best_deg = t_deg[t];
                }
            }
            if (best != NIL) { s_mate[s] = best; t_mate[best] = s; cnt++; }
        }
        return cnt;
    }

    /* BFS: compute shortest augmenting path length */
    int bfs_shortest_level(
    #ifdef STATS
        PhaseStats& ps
    #endif
    ) {
        int shortest = INF_LEVEL;

        int sf = s_exposed.first();
        while (sf != NIL) {
        #ifdef STATS
            ps.bfs_vtx++;
        #endif
            s_level[sf] = 0;
            s_bfs_queue.push(sf);
            s_stt[sf] = STT_BFS_QUEUED;
            sf = s_exposed.next(sf);
        }

        int cur_layer = 0;
        while (!s_bfs_queue.empty()) {
            int s = s_bfs_queue.front();
            s_bfs_queue.pop();
            s_done_queue.push(s);
            s_stt[s] = STT_BFS_DONE;

            if (cur_layer < s_level[s] / 2) {
                if (shortest != INF_LEVEL) break;
                cur_layer++;
            }

            int num_edg = (int)src_adj[s].size();
            const int* adj = num_edg > 0 ? src_adj[s].data() : nullptr;
            for (int i = 0; i < num_edg; i++) {
                int t = adj[i];
            #ifdef STATS
                ps.bfs_edg++;
            #endif
                if (t_stt[t] == STT_BFS_DONE) continue;
            #ifdef STATS
                ps.bfs_vtx++;
            #endif
                t_level[t] = s_level[s] + 1;
                t_done_queue.push(t);
                t_stt[t] = STT_BFS_DONE;

                int ss = t_mate[t];
                if (ss == NIL) {
                    shortest = t_level[t];
                } else {
                #ifdef STATS
                    ps.bfs_edg++;
                    ps.bfs_vtx++;
                #endif
                    s_level[ss] = t_level[t] + 1;
                    s_bfs_queue.push(ss);
                    s_stt[ss] = STT_BFS_QUEUED;
                }
            }
        }
        return shortest;
    }

    /* DFS: find maximal set of shortest augmenting paths (HK) */
    void dfs_find_paths_hk(int shortest
    #ifdef STATS
        , PhaseStats& ps
    #endif
    ) {
        int sf = s_exposed.first();
        while (sf != NIL) {
            int sf_next = s_exposed.next(sf);
        #ifdef STATS
            ps.dfs_vtx++;
        #endif
            s_dfs_stack.push(sf);
            s_stt[sf] = STT_DFS_ACTIVE;

            while (!s_dfs_stack.empty()) {
                int s = s_dfs_stack.top();
                int num_edg = (int)src_adj[s].size();
                const int* adj = num_edg > 0 ? src_adj[s].data() : nullptr;

                while (s_idx[s] < num_edg) {
                    int t = adj[s_idx[s]];
                #ifdef STATS
                    ps.dfs_edg++;
                #endif
                    int ss = t_mate[t];
                    if ((t_level[t] != s_level[s] + 1) ||
                        ((ss == NIL) && (t_stt[t] == STT_LAST)) ||
                        ((ss != NIL) && ((s_level[s] + 1 == shortest) ||
                                         (s_stt[ss] == STT_DFS_ACTIVE) ||
                                         (s_stt[ss] == STT_DFS_DONE)))) {
                        s_idx[s]++;
                    } else {
                        break;
                    }
                }

                if (s_idx[s] < num_edg) {
                    int t = adj[s_idx[s]];
                    int ss = t_mate[t];
                #ifdef STATS
                    ps.dfs_vtx++;
                #endif
                    if (ss == NIL) {
                        s_last_queue.push(s);
                        t_last_queue.push(t);
                        t_stt[t] = STT_LAST;
                        while (!s_dfs_stack.empty()) {
                            s_stt[s_dfs_stack.top()] = STT_DFS_DONE;
                            s_dfs_stack.pop();
                        }
                        break;
                    }
                #ifdef STATS
                    ps.dfs_edg++;
                    ps.dfs_vtx++;
                #endif
                    s_ptr[ss] = s;
                    s_dfs_stack.push(ss);
                    s_stt[ss] = STT_DFS_ACTIVE;
                    s_idx[s]++;
                } else {
                    s_dfs_stack.pop();
                    s_stt[s] = STT_DFS_DONE;
                }
            }
            sf = sf_next;
        }
    }

    /* DFS: find augmenting paths with lookahead (no BFS layering, matchbox-style) */
    void dfs_find_paths_lkhd(
    #ifdef STATS
        PhaseStats& ps
    #endif
    ) {
        int sf = s_exposed.first();
        while (sf != NIL) {
            int sf_next = s_exposed.next(sf);
        #ifdef STATS
            ps.dfs_vtx++;
        #endif
            s_dfs_stack.push(sf);
            s_stt[sf] = STT_DFS_ACTIVE;

            while (!s_dfs_stack.empty()) {
                int s = s_dfs_stack.top();
                int num_edg = (int)src_adj[s].size();
                const int* adj = num_edg > 0 ? src_adj[s].data() : nullptr;

                /* Lookahead: check for direct free t-neighbor */
                while (s_lkhd[s] < num_edg) {
                    int t = adj[s_lkhd[s]];
                #ifdef STATS
                    ps.bfs_edg++;
                #endif
                    int ss = t_mate[t];
                    if ((ss == NIL && t_stt[t] != STT_IDLE) ||
                        ss != NIL) {
                        s_lkhd[s]++;
                    } else {
                        break;
                    }
                }

                if (s_lkhd[s] < num_edg) {
                    int t = adj[s_lkhd[s]];
                #ifdef STATS
                    ps.bfs_vtx++;
                #endif
                    s_last_queue.push(s);
                    t_last_queue.push(t);
                    t_stt[t] = STT_LAST;
                    while (!s_dfs_stack.empty()) {
                        int sss = s_dfs_stack.top();
                        s_dfs_stack.pop();
                        s_done_queue.push(sss);
                        s_stt[sss] = STT_DFS_DONE;
                    }
                    s_lkhd[s]++;
                    break;
                }

                /* Regular DFS deeper */
                while (s_idx[s] < num_edg) {
                    int t = adj[s_idx[s]];
                #ifdef STATS
                    ps.dfs_edg++;
                #endif
                    int ss = t_mate[t];
                    if (ss == NIL ||
                        (ss != NIL && s_stt[ss] != STT_IDLE)) {
                        s_idx[s]++;
                    } else {
                        break;
                    }
                }

                if (s_idx[s] < num_edg) {
                    int t = adj[s_idx[s]];
                    int ss = t_mate[t];
                #ifdef STATS
                    ps.dfs_vtx++;
                    ps.dfs_edg++;
                    ps.dfs_vtx++;
                #endif
                    s_ptr[ss] = s;
                    s_dfs_stack.push(ss);
                    s_stt[ss] = STT_DFS_ACTIVE;
                    s_idx[s]++;
                } else {
                    s_dfs_stack.pop();
                    s_done_queue.push(s);
                    s_stt[s] = STT_DFS_DONE;
                }
            }
            sf = sf_next;
        }
    }

    /* Augment along path from s_last to root via s_ptr */
    int augment(int s_last, int t_last) {
        int s = s_last, t = t_last, k = 0;
        while (s != NIL) {
            int tt = s_mate[s];
            s_mate[s] = t; t_mate[t] = s;
            if (s_ptr[s] == NIL) s_exposed.erase(s);
            s = s_ptr[s]; t = tt; k++;
        }
        return 2 * k - 1;
    }

    /* Cleanup after HK phase (also resets lookahead indices) */
    void cleanup_hk() {
        while (!s_bfs_queue.empty()) {
            int s = s_bfs_queue.front(); s_bfs_queue.pop();
            s_ptr[s] = NIL; s_level[s] = INF_LEVEL;
            s_stt[s] = STT_IDLE; s_idx[s] = 0; s_lkhd[s] = 0;
        }
        while (!s_done_queue.empty()) {
            int s = s_done_queue.front(); s_done_queue.pop();
            s_ptr[s] = NIL; s_level[s] = INF_LEVEL;
            s_stt[s] = STT_IDLE; s_idx[s] = 0; s_lkhd[s] = 0;
        }
        while (!t_done_queue.empty()) {
            int t = t_done_queue.front(); t_done_queue.pop();
            t_level[t] = INF_LEVEL; t_stt[t] = STT_IDLE;
        }
    }

    /* Cleanup after lookahead phase */
    void cleanup_lkhd() {
        while (!s_done_queue.empty()) {
            int s = s_done_queue.front(); s_done_queue.pop();
            s_ptr[s] = NIL; s_stt[s] = STT_IDLE;
            s_idx[s] = 0; s_lkhd[s] = 0;
        }
    }

    /* Main solver: HK mode */
    int solve_hk(int greedy_mode) {
        if (greedy_mode == 1) greedy_size = greedy_init();
        else if (greedy_mode == 2) greedy_size = greedy_init_md();
        int card = greedy_size;

    #ifdef STATS
        stats.reset();
        stats.greedy_card = greedy_size;
        stats.reversed = rvrs;
    #endif

        for (int s = 0; s < s_count; s++)
            if (s_mate[s] == NIL) s_exposed.push(s);

        int phase_count = 0;
        while (true) {
        #ifdef STATS
            PhaseStats ps = {0, 0, 0, 0, 0, 0, 0, INF_LEVEL, 0};
        #endif

            int shortest = bfs_shortest_level(
            #ifdef STATS
                ps
            #endif
            );

            if (shortest == INF_LEVEL) {
                cleanup_hk();
                break;
            }

            phase_count++;
        #ifdef STATS
            ps.shortest_path_len = shortest;
        #endif

            dfs_find_paths_hk(shortest
            #ifdef STATS
                , ps
            #endif
            );

            while (!t_last_queue.empty()) {
                int sl = s_last_queue.front(); s_last_queue.pop();
                int tl = t_last_queue.front(); t_last_queue.pop();
                int path_len = augment(sl, tl);
                card++;
            #ifdef STATS
                ps.num_augmentations++;
                ps.agg_aug_path_len += path_len;
                if (path_len < ps.min_aug_path_len) ps.min_aug_path_len = path_len;
                if (path_len > ps.max_aug_path_len) ps.max_aug_path_len = path_len;
            #endif
            }

            cleanup_hk();

        #ifdef STATS
            stats.phases.push_back(ps);
        #endif
        }

    #ifdef STATS
        stats.print();
    #endif
        printf("Phases: %d\n", phase_count);
        return card;
    }

    /* Main solver: lookahead mode (DFS-only, matchbox-style) */
    int solve_lkhd(int greedy_mode) {
        if (greedy_mode == 1) greedy_size = greedy_init();
        else if (greedy_mode == 2) greedy_size = greedy_init_md();
        int card = greedy_size;

    #ifdef STATS
        stats.reset();
        stats.greedy_card = greedy_size;
        stats.reversed = rvrs;
    #endif

        for (int s = 0; s < s_count; s++)
            if (s_mate[s] == NIL) s_exposed.push(s);

        int pass_count = 0;
        while (true) {
        #ifdef STATS
            PhaseStats ps = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        #endif

            dfs_find_paths_lkhd(
            #ifdef STATS
                ps
            #endif
            );

            if (t_last_queue.empty()) break;

            pass_count++;
            while (!t_last_queue.empty()) {
                int sl = s_last_queue.front(); s_last_queue.pop();
                int tl = t_last_queue.front(); t_last_queue.pop();
                t_stt[tl] = STT_IDLE;
                int path_len = augment(sl, tl);
                card++;
            #ifdef STATS
                ps.num_augmentations++;
                ps.agg_aug_path_len += path_len;
            #endif
            }

            /* Cleanup: reset only visited vertices */
            while (!s_done_queue.empty()) {
                int s = s_done_queue.front(); s_done_queue.pop();
                s_ptr[s] = NIL; s_stt[s] = STT_IDLE;
                s_idx[s] = 0; s_lkhd[s] = 0;
            }

        #ifdef STATS
            stats.phases.push_back(ps);
        #endif
        }

    #ifdef STATS
        stats.print();
    #endif
        printf("Phases: %d\n", pass_count);
        return card;
    }

    /* Get matching in original (left, right) coordinates */
    std::vector<std::pair<int,int>> get_matching() const {
        std::vector<std::pair<int,int>> result;
        for (int s = 0; s < s_count; s++) {
            if (s_mate[s] != NIL) {
                if (rvrs) result.push_back({s_mate[s], s});
                else      result.push_back({s, s_mate[s]});
            }
        }
        std::sort(result.begin(), result.end());
        return result;
    }
};

/* =========================================================================
 * Validation
 * ========================================================================= */
void validate_matching(int lc, int rc,
                       const std::vector<std::vector<int>>& graph,
                       const std::vector<std::pair<int,int>>& matching,
                       bool rvrs) {
    std::vector<int> ldeg(lc, 0), rdeg(rc, 0);
    int errors = 0;
    for (auto& e : matching) {
        const auto& adj = rvrs ? graph[e.second] : graph[e.first];
        int target = rvrs ? e.first : e.second;
        if (!std::binary_search(adj.begin(), adj.end(), target)) {
            fprintf(stderr, "ERROR: Edge (%d, %d) not in graph!\n", e.first, e.second);
            errors++;
        }
        ldeg[e.first]++; rdeg[e.second]++;
    }
    for (int i = 0; i < lc; i++)
        if (ldeg[i] > 1) { fprintf(stderr, "ERROR: Left %d in %d edges!\n", i, ldeg[i]); errors++; }
    for (int i = 0; i < rc; i++)
        if (rdeg[i] > 1) { fprintf(stderr, "ERROR: Right %d in %d edges!\n", i, rdeg[i]); errors++; }
    int ml = 0, mr = 0;
    for (int i = 0; i < lc; i++) if (ldeg[i] > 0) ml++;
    for (int i = 0; i < rc; i++) if (rdeg[i] > 0) mr++;

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %d\n", (int)matching.size());
    printf("Left matched: %d, Right matched: %d\n", ml, mr);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

/* =========================================================================
 * Main
 * ========================================================================= */
int main(int argc, char* argv[]) {
    printf("Hopcroft-Karp Pure Algorithm - C++ Implementation\n");
    printf("==================================================\n\n");

    if (argc < 2) {
        printf("Usage: %s <filename> [--greedy|--greedy-md] [--lkhd] [--rvrs|--fwd]\n", argv[0]);
        return 1;
    }

    int greedy_mode = 0;
    bool use_lkhd = false;
    bool force_rvrs = false;
    bool auto_dir = true;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedy_mode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedy_mode = 2;
        else if (std::string(argv[i]) == "--lkhd") use_lkhd = true;
        else if (std::string(argv[i]) == "--rvrs") { force_rvrs = true; auto_dir = false; }
        else if (std::string(argv[i]) == "--fwd") { force_rvrs = false; auto_dir = false; }
    }

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open file: %s\n", argv[1]); return 1; }

    int lc, rc, m;
    if (fscanf(f, "%d %d %d", &lc, &rc, &m) != 3) {
        fprintf(stderr, "Bad header (expected: left_count right_count edge_count)\n");
        fclose(f); return 1;
    }

    std::vector<std::pair<int,int>> edges;
    edges.reserve(m);
    for (int i = 0; i < m; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edges.push_back({u, v});
    }
    fclose(f);

    printf("Graph: %d left, %d right, %d edges\n", lc, rc, (int)edges.size());

    bool use_rvrs = force_rvrs;
    if (auto_dir) use_rvrs = (rc < lc);

    printf("Mode: %s\n", use_lkhd ? "lookahead DFS" : "HK (BFS+DFS)");
    printf("Direction: %s%s\n", use_rvrs ? "T->S (reversed)" : "S->T (normal)",
           auto_dir ? " (auto)" : " (forced)");

    auto t0 = std::chrono::high_resolution_clock::now();
    HopcroftKarpPure hk;
    hk.build(lc, rc, edges, use_rvrs);
    int card;
    if (use_lkhd) card = hk.solve_lkhd(greedy_mode);
    else          card = hk.solve_hk(greedy_mode);
    auto t1 = std::chrono::high_resolution_clock::now();

    auto matching = hk.get_matching();
    validate_matching(lc, rc, hk.src_adj, matching, use_rvrs);

    printf("Matching size: %d\n", card);
    if (greedy_mode > 0) {
        printf("Greedy init size: %d\n", hk.greedy_size);
        if (card > 0) printf("Greedy/Final: %.2f%%\n", 100.0 * hk.greedy_size / card);
    }
    printf("Time: %ld ms\n",
           (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

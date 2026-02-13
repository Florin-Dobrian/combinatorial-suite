/*
 * Gabow's Scaling Algorithm (Optimized) - O(E√V) Maximum Matching
 *
 * Pure cardinality (unweighted) version — integer weights conceptually all 1.
 *
 * Phase 1: BFS by levels (Delta), detect blossoms.
 *          Build contracted graph H: edges connecting different dbase
 *          components that were processed during BFS.
 * Phase 2: Find all shortest augmenting paths in H (iterative DFS
 *          with blossom contraction), unfold to G via bridges.
 *
 * Based on LEDA-7's mc_matching_gabow architecture, stripped of weighted
 * dual machinery (dval, bd, bDelta, priority queue) which is incorrect
 * for pure cardinality at large Delta.
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <chrono>

static const int NIL = -1;
static const int UNLABELED = 0;
static const int EVEN = 1;
static const int ODD = 2;

struct GabowOptimized {
    int n;
    std::vector<std::vector<int>> graph;
    std::vector<int> mate;

    /* phase_1 BFS tree */
    std::vector<int> label;
    std::vector<int> parent;
    std::vector<int> source_bridge;
    std::vector<int> target_bridge;

    /* base union-find (immediate unions during shrink_path) */
    std::vector<int> base_par;

    /* dbase union-find (deferred unions at Delta boundaries) */
    std::vector<int> dbase_par;

    /* BFS level queue: edges to process at each Delta */
    std::vector<std::vector<std::pair<int,int>>> level_queue;

    /* interleaved LCA with epoch (size_t to avoid overflow) */
    std::vector<size_t> lca_tag1, lca_tag2;
    size_t lca_epoch;

    /* tree membership */
    std::vector<bool> in_tree;
    std::vector<int> tree_nodes;

    int Delta;

    /* H construction: mark which edges are in H */
    /* For each tree vertex u, we mark edges to other dbase components */
    /* Using epoch-based adjacency: h_edge_epoch[v] == h_epoch means v's
       edges to other dbase components are H-edges */
    /* Actually simpler: just build H-adjacency during phase_1 */

    /* phase_2 / H state */
    std::vector<int> rep;         /* rep[v] = dbase(v) at start of phase_2 */
    std::vector<int> mateH;
    std::vector<int> labelH;
    std::vector<int> parentH_src, parentH_tgt;
    std::vector<int> bridgeH_src, bridgeH_tgt;
    std::vector<int> dirH;
    std::vector<int> even_timeH;
    int tH;
    std::vector<int> dbase2_par;  /* blossoms in H */
    std::vector<std::vector<int>> contracted_into;

    GabowOptimized(int n_, const std::vector<std::pair<int,int>>& edges) : n(n_) {
        graph.resize(n);
        mate.assign(n, NIL);
        label.assign(n, UNLABELED);
        parent.assign(n, NIL);
        source_bridge.assign(n, NIL);
        target_bridge.assign(n, NIL);
        base_par.resize(n);
        dbase_par.resize(n);
        lca_tag1.assign(n, 0);
        lca_tag2.assign(n, 0);
        lca_epoch = 0;
        in_tree.assign(n, false);

        rep.resize(n);
        mateH.assign(n, NIL);
        labelH.assign(n, UNLABELED);
        parentH_src.assign(n, NIL);
        parentH_tgt.assign(n, NIL);
        bridgeH_src.assign(n, NIL);
        bridgeH_tgt.assign(n, NIL);
        dirH.assign(n, 0);
        even_timeH.assign(n, 0);
        tH = 0;
        dbase2_par.resize(n);
        contracted_into.resize(n);

        for (auto& e : edges) {
            int u = e.first, v = e.second;
            if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
                graph[u].push_back(v);
                graph[v].push_back(u);
            }
        }
        for (int i = 0; i < n; i++) {
            std::sort(graph[i].begin(), graph[i].end());
            graph[i].erase(std::unique(graph[i].begin(), graph[i].end()), graph[i].end());
        }
        level_queue.resize(n + 2);
    }

    /* ---- union-find: base ---- */
    int find_base(int v) {
        while (base_par[v] != v) { base_par[v] = base_par[base_par[v]]; v = base_par[v]; }
        return v;
    }
    void union_base(int a, int b, int r) {
        a = find_base(a); b = find_base(b);
        base_par[a] = r; base_par[b] = r;
    }

    /* ---- union-find: dbase ---- */
    int find_dbase(int v) {
        while (dbase_par[v] != v) { dbase_par[v] = dbase_par[dbase_par[v]]; v = dbase_par[v]; }
        return v;
    }
    void union_dbase(int a, int b) {
        a = find_dbase(a); b = find_dbase(b);
        if (a != b) dbase_par[a] = b;
    }
    void make_rep_dbase(int v) {
        int r = find_dbase(v);
        if (r != v) { dbase_par[r] = v; dbase_par[v] = v; }
    }

    /* ---- union-find: dbase2 (blossoms in H) ---- */
    int find_db2(int v) {
        while (dbase2_par[v] != v) { dbase2_par[v] = dbase2_par[dbase2_par[v]]; v = dbase2_par[v]; }
        return v;
    }
    void union_db2(int a, int b) {
        a = find_db2(a); b = find_db2(b);
        if (a != b) dbase2_par[a] = b;
    }
    void make_rep_db2(int v) {
        int r = find_db2(v);
        if (r != v) { dbase2_par[r] = v; dbase2_par[v] = v; }
    }

    /* ---- interleaved LCA ---- */
    int find_lca(int u, int v) {
        ++lca_epoch;
        size_t ep = lca_epoch;
        int hx = find_base(u), hy = find_base(v);
        lca_tag1[hx] = ep;
        lca_tag2[hy] = ep;
        while (true) {
            if (lca_tag1[hy] == ep) return hy;
            if (lca_tag2[hx] == ep) return hx;
            bool hxr = (mate[hx] == NIL || parent[mate[hx]] == NIL);
            bool hyr = (mate[hy] == NIL || parent[mate[hy]] == NIL);
            if (hxr && hyr) return NIL;
            if (!hxr) { hx = find_base(parent[mate[hx]]); lca_tag1[hx] = ep; }
            if (!hyr) { hy = find_base(parent[mate[hy]]); lca_tag2[hy] = ep; }
        }
    }

    /* Called during phase_1 when processing a non-matching edge (z, u) 
       that connects different dbase components */

    /* ---- shrink_path ---- */
    void shrink_path(int b, int x, int y,
                     std::vector<std::pair<int,int>>& dunions) {
        int v = find_base(x);
        while (v != b) {
            union_base(v, b, b);
            dunions.push_back({v, b});
            int mv = mate[v];
            union_base(mv, b, b);
            dunions.push_back({mv, b});
            base_par[b] = b;
            source_bridge[mv] = x;
            target_bridge[mv] = y;
            /* Scan newly-EVEN vertex for edges at next Delta level */
            for (int w : graph[mv]) {
                if (w == mate[mv]) continue;
                int bw = find_base(w);
                if (label[bw] == ODD) continue;
                if (label[bw] == UNLABELED) {
                    level_queue[Delta + 1].push_back({mv, w});
                } else if (label[bw] == EVEN) {
                    level_queue[Delta].push_back({mv, w});
                }
            }
            v = find_base(parent[mv]);
        }
        dunions.push_back({b, b});
    }

    /* ================================================================ */
    /*                          PHASE 1                                 */
    /* ================================================================ */
    bool phase_1() {
        Delta = 0;
        tree_nodes.clear();
        for (auto& q : level_queue) q.clear();
        std::vector<std::pair<int,int>> dunions;

        for (int i = 0; i < n; i++) {
            base_par[i] = i;
            dbase_par[i] = i;
            label[i] = UNLABELED;
            parent[i] = NIL;
            source_bridge[i] = NIL;
            target_bridge[i] = NIL;
            in_tree[i] = false;
        }

        /* Initialize: free vertices are EVEN roots at Delta=0 */
        for (int v = 0; v < n; v++) {
            if (mate[v] == NIL) {
                label[v] = EVEN;
                in_tree[v] = true;
                tree_nodes.push_back(v);
                for (int u : graph[v]) {
                    if (u == mate[v]) continue;
                    int bu = find_base(u);
                    if (label[bu] == ODD) continue;
                    if (label[bu] == UNLABELED)
                        level_queue[1].push_back({v, u});  /* grows at Delta=1 */
                    else if (label[bu] == EVEN)
                        level_queue[0].push_back({v, u});  /* EVEN-EVEN at Delta=0 */
                }
            }
        }

        bool found_sap = false;

        while (Delta <= n) {
            while (!level_queue[Delta].empty()) {
                auto [z, u] = level_queue[Delta].back();
                level_queue[Delta].pop_back();

                int bz = find_base(z), bu = find_base(u);
                if (label[bz] != EVEN) { std::swap(z, u); std::swap(bz, bu); }
                if (bz == bu || label[bz] != EVEN) continue;
                if (u == mate[z] || label[bu] == ODD) continue;

                if (label[bu] == UNLABELED) {
                    int mv = mate[u];
                    if (mv == NIL) continue;
                    parent[u] = z;
                    parent[mv] = u;
                    label[u] = ODD;
                    label[mv] = EVEN;
                    in_tree[u] = true;
                    in_tree[mv] = true;
                    tree_nodes.push_back(u);
                    tree_nodes.push_back(mv);
                    /* Record the grow edge as H-edge */
                    /* Scan from newly EVEN vertex mv */
                    for (int w : graph[mv]) {
                        if (w == mate[mv]) continue;
                        int bw = find_base(w);
                        if (label[bw] == ODD) continue;
                        if (label[bw] == UNLABELED)
                            level_queue[Delta + 1].push_back({mv, w});
                        else if (label[bw] == EVEN)
                            level_queue[Delta].push_back({mv, w});
                    }

                } else if (label[bu] == EVEN) {
                    int lca = find_lca(z, u);
                    if (lca != NIL) {
                        /* Blossom — record the shrink edge as H-edge */
                        shrink_path(lca, z, u, dunions);
                        shrink_path(lca, u, z, dunions);
                    } else {
                        /* Augmenting path found */
                        found_sap = true;
                        /* DON'T break — continue all edges at this Delta */
                    }
                }
            }

            if (found_sap) {
                /* Build H: contracted_into and mateH */
                for (int v : tree_nodes) {
                    int db = find_dbase(v);
                    contracted_into[db].push_back(v);
                    mateH[v] = NIL;
                }
                /* Set mateH for matching edges between different dbase components */
                for (int u : tree_nodes) {
                    int uh = find_dbase(u);
                    int mv = mate[u];
                    if (mv != NIL && in_tree[mv]) {
                        int vh = find_dbase(mv);
                        if (uh != vh) {
                            mateH[uh] = vh;
                            mateH[vh] = uh;
                        }
                    }
                }
                return true;
            }

            /* Execute deferred dbase unions for this Delta */
            for (auto& [a, b] : dunions) {
                if (a == b) make_rep_dbase(a);
                else union_dbase(a, b);
            }
            dunions.clear();
            Delta++;
        }
        return false;
    }

    /* ================================================================ */
    /*                          PHASE 2                                 */
    /* ================================================================ */

    /* find_apHG: ITERATIVE DFS in H to find augmenting path.
     * Returns the free H-node found, or NIL.
     * Uses explicit stack to handle 500k+ vertices. */
    int find_apHG(int root_vh) {
        /* Stack frame for iterative DFS */
        struct Frame {
            int vh;           /* current H-node being scanned */
            int ci_idx;       /* index into contracted_into[vh] */
            int adj_idx;      /* index into graph[v] */
            int v;            /* current G-vertex being scanned */
        };
        std::vector<Frame> stk;
        stk.push_back({root_vh, 0, 0, -1});

        while (!stk.empty()) {
            auto& f = stk.back();
            int vh = f.vh;

            /* Iterate through G-vertices in this H-node */
            while (f.ci_idx < (int)contracted_into[vh].size()) {
                f.v = contracted_into[vh][f.ci_idx];
                int v = f.v;

                /* Iterate through edges from this G-vertex */
                while (f.adj_idx < (int)graph[v].size()) {
                    int w = graph[v][f.adj_idx];
                    f.adj_idx++;

                    if (!in_tree[w]) continue;
                    if (mate[v] == w) continue;  /* skip matching edges */
                    int wh = find_dbase(w);
                    if (wh == find_dbase(v)) continue;  /* same H-node */
                    int uh = find_db2(rep[w]);
                    if (mateH[vh] == uh) continue;  /* skip mate (raw, not through db2) */

                    if (labelH[uh] == ODD) continue;  /* already ODD */

                    if (labelH[uh] == UNLABELED) {
                        int muh = mateH[uh];
                        if (muh == NIL) {
                            /* Free node — augmenting path found! */
                            labelH[uh] = ODD;
                            parentH_src[uh] = w;
                            parentH_tgt[uh] = v;
                            return uh;
                        }
                        /* Grow step: extend by two edges */
                        labelH[uh] = ODD;
                        parentH_src[uh] = w;
                        parentH_tgt[uh] = v;
                        labelH[muh] = EVEN;
                        even_timeH[muh] = tH++;
                        /* Push current state and recurse into muh */
                        stk.push_back({muh, 0, 0, -1});
                        goto next_frame;

                    } else if (labelH[uh] == EVEN) {
                        /* Blossom step */
                        int bh = find_db2(vh);
                        int zh = find_db2(uh);
                        if (even_timeH[bh] < even_timeH[zh]) {
                            std::vector<int> tmp, endpoints;
                            int cur = zh;
                            while (cur != bh) {
                                endpoints.push_back(cur);
                                int mc = mateH[cur];
                                endpoints.push_back(mc);
                                tmp.push_back(mc);
                                int ps = parentH_src[mc], pt = parentH_tgt[mc];
                                int next = rep[rep[ps] == mc ? pt : ps];
                                cur = find_db2(next);
                            }
                            for (int nd : endpoints) union_db2(nd, bh);
                            make_rep_db2(bh);

                            for (int mc : tmp) {
                                bridgeH_src[mc] = v;
                                bridgeH_tgt[mc] = w;
                                dirH[mc] = -1;
                            }
                            /* Push each new ODD node for scanning */
                            for (int i = (int)tmp.size() - 1; i >= 0; i--) {
                                stk.push_back({tmp[i], 0, 0, -1});
                            }
                            goto next_frame;
                        }
                    }
                }
                /* Done with this G-vertex, move to next */
                f.ci_idx++;
                f.adj_idx = 0;
            }
            /* Done with all G-vertices in vh — backtrack */
            stk.pop_back();
            continue;

            next_frame:;
        }
        return NIL;
    }

    /* trace_H_path: iterative, collects non-matching G-edges along H-path */
    void trace_H_path(int vh, int uh, std::vector<std::pair<int,int>>& edges_out) {
        struct Frame { int vh, uh, phase, bs, bt, side_a, side_b; };
        std::vector<Frame> stk;
        stk.push_back({vh, uh, 0, 0, 0, 0, 0});

        while (!stk.empty()) {
            auto& f = stk.back();
            if (f.vh == f.uh) { stk.pop_back(); continue; }

            if (labelH[f.vh] == EVEN) {
                int mvh = mateH[f.vh];
                int ps = parentH_src[mvh], pt = parentH_tgt[mvh];
                edges_out.push_back({ps, pt});
                f.vh = rep[rep[ps] == mvh ? pt : ps];
                continue;
            }
            if (f.phase == 0) {
                f.bs = bridgeH_src[f.vh];
                f.bt = bridgeH_tgt[f.vh];
                if (dirH[f.vh] == 1) {
                    f.side_a = rep[f.bs];
                    f.side_b = rep[f.bt];
                } else {
                    f.side_a = rep[f.bt];
                    f.side_b = rep[f.bs];
                }
                f.phase = 1;
                int mt = (mateH[f.vh] != NIL) ? rep[mateH[f.vh]] : f.vh;
                stk.push_back({f.side_a, mt, 0, 0, 0, 0, 0});
                continue;
            }
            if (f.phase == 1) {
                edges_out.push_back({f.bs, f.bt});
                f.phase = 2;
                stk.push_back({f.side_b, f.uh, 0, 0, 0, 0, 0});
                continue;
            }
            stk.pop_back();
        }
    }

    /* find_path_in_G: iterative unfold within single H-node */
    void find_path_in_G(int v, int u, std::vector<std::pair<int,int>>& pairs) {
        struct Frame { int v, u, phase, sb, tb; };
        std::vector<Frame> stk;
        stk.push_back({v, u, 0, 0, 0});

        while (!stk.empty()) {
            auto& f = stk.back();
            if (f.v == f.u) { stk.pop_back(); continue; }
            if (f.phase == 0) {
                if (label[f.v] == EVEN) {
                    int mv = mate[f.v], pmv = parent[mv];
                    pairs.push_back({mv, pmv});
                    f.v = pmv;
                    continue;
                }
                f.sb = source_bridge[f.v];
                f.tb = target_bridge[f.v];
                f.phase = 1;
                stk.push_back({f.sb, mate[f.v], 0, 0, 0});
                continue;
            }
            if (f.phase == 1) {
                pairs.push_back({f.sb, f.tb});
                f.phase = 2;
                stk.push_back({f.tb, f.u, 0, 0, 0});
                continue;
            }
            stk.pop_back();
        }
    }

    /* augmentG: unfold H-edges to G and augment */
    void augmentG(const std::vector<std::pair<int,int>>& h_edges) {
        std::vector<std::pair<int,int>> pairs;
        for (auto& [u, v] : h_edges) {
            pairs.push_back({u, v});
            find_path_in_G(u, rep[u], pairs);
            find_path_in_G(v, rep[v], pairs);
        }
        for (auto& [a, b] : pairs) {
            mate[a] = b;
            mate[b] = a;
        }
    }

    /* phase_2: find all SAPs in H, unfold and augment */
    void phase_2() {
        for (int v : tree_nodes) {
            rep[v] = find_dbase(v);
            labelH[v] = UNLABELED;
            parentH_src[v] = parentH_tgt[v] = NIL;
            bridgeH_src[v] = bridgeH_tgt[v] = NIL;
            dirH[v] = 0;
            even_timeH[v] = 0;
            dbase2_par[v] = v;
        }
        tH = 0;

        std::vector<std::vector<std::pair<int,int>>> all_paths;

        for (int vh : tree_nodes) {
            if (vh != rep[vh]) continue;
            if (labelH[vh] != UNLABELED || mateH[vh] != NIL) continue;

            labelH[vh] = EVEN;
            even_timeH[vh] = tH++;

            int free_node = find_apHG(vh);
            if (free_node != NIL) {
                std::vector<std::pair<int,int>> h_nm;
                int ps = parentH_src[free_node], pt = parentH_tgt[free_node];
                h_nm.push_back({ps, pt});
                int next = rep[rep[ps] == free_node ? pt : ps];
                trace_H_path(next, vh, h_nm);
                all_paths.push_back(std::move(h_nm));
            }
        }

        for (auto& he : all_paths) augmentG(he);

        /* Clean up */
        for (int v : tree_nodes) {
            int db = find_dbase(v);
            contracted_into[db].clear();
            contracted_into[v].clear();
            mateH[v] = NIL;
        }
    }

    /* ================================================================ */
    /*                      MAIN ENTRY POINT                            */
    /* ================================================================ */
    std::vector<std::pair<int,int>> maximum_matching() {
        /* greedy init */
        for (int u = 0; u < n; u++) {
            if (mate[u] != NIL) continue;
            for (int v : graph[u]) {
                if (mate[v] == NIL) { mate[u] = v; mate[v] = u; break; }
            }
        }
        while (phase_1()) phase_2();

        std::vector<std::pair<int,int>> result;
        for (int u = 0; u < n; u++)
            if (mate[u] != NIL && mate[u] > u)
                result.push_back({u, mate[u]});
        std::sort(result.begin(), result.end());
        return result;
    }
};

/* ================================================================ */
/*                    VALIDATION AND MAIN                            */
/* ================================================================ */

void validate_matching(int n, const std::vector<std::vector<int>>& graph,
                       const std::vector<std::pair<int,int>>& matching) {
    std::vector<int> deg(n, 0);
    int errors = 0;
    for (auto& e : matching) {
        if (!std::binary_search(graph[e.first].begin(), graph[e.first].end(), e.second)) {
            fprintf(stderr, "ERROR: Edge (%d, %d) not in graph!\n", e.first, e.second);
            errors++;
        }
        deg[e.first]++;
        deg[e.second]++;
    }
    for (int i = 0; i < n; i++)
        if (deg[i] > 1) { fprintf(stderr, "ERROR: Vertex %d in %d edges!\n", i, deg[i]); errors++; }
    int matched = 0;
    for (int i = 0; i < n; i++) if (deg[i] > 0) matched++;

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %d\n", (int)matching.size());
    printf("Matched vertices: %d\n", matched);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

int main(int argc, char* argv[]) {
    printf("Gabow's Scaling Algorithm (Optimized) - C++ Implementation\n");
    printf("============================================================\n\n");

    if (argc < 2) { printf("Usage: %s <filename>\n", argv[0]); return 1; }

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open file: %s\n", argv[1]); return 1; }

    int n, m;
    if (fscanf(f, "%d %d", &n, &m) != 2) { fprintf(stderr, "Bad header\n"); fclose(f); return 1; }

    std::vector<std::pair<int,int>> edges;
    edges.reserve(m);
    for (int i = 0; i < m; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edges.push_back({u, v});
    }
    fclose(f);

    printf("Graph: %d vertices, %d edges\n", n, (int)edges.size());

    auto t0 = std::chrono::high_resolution_clock::now();
    GabowOptimized gabow(n, edges);
    auto matching = gabow.maximum_matching();
    auto t1 = std::chrono::high_resolution_clock::now();

    validate_matching(n, gabow.graph, matching);

    printf("Matching size: %d\n", (int)matching.size());
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

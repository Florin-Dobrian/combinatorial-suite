/*
 * Gabow MCM with proper duals and edge IDs.
 * O(E * sqrt(V)) Maximum Cardinality Matching.
 *
 * Edge model: each undirected edge gets a unique ID with fixed
 * (src, tgt) endpoints. adj[v] stores edge IDs, not neighbors.
 * This enables exact translation of LEDA's Phase 2 edge semantics:
 *   G.source(e) -> esrc[e]
 *   G.target(e) -> etgt[e]
 *   G.opposite(v, e) -> (esrc[e]==v ? etgt[e] : esrc[e])
 *
 * All integers. No floating point in algorithm. No dependencies.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <climits>
#include <chrono>

static const int NIL = -1;
enum { EVEN = 0, ODD = 1, UNLABELED = 2 };

struct GabowMCM {
    int n, m;
    int greedy_size = 0;

    /* Edge storage: esrc[eid], etgt[eid] are fixed endpoints */
    std::vector<int> esrc, etgt;
    /* adj[v] = list of edge IDs incident to v */
    std::vector<std::vector<int>> adj;
    std::vector<int> mate;

    /* Helpers */
    int opposite(int v, int eid) { return esrc[eid] == v ? etgt[eid] : esrc[eid]; }
    int w_edge(int eid) { return (mate[esrc[eid]] == etgt[eid]) ? 2 : 0; }

    /* ---- Phase 1 ---- */
    std::vector<int> label;
    std::vector<int> parent;             /* parent vertex in alternating tree */
    std::vector<int> source_bridge, target_bridge;
    std::vector<int> bd, bDelta;

    std::vector<int> base_par;
    int find_base(int v) {
        while (base_par[v] != v) { base_par[v] = base_par[base_par[v]]; v = base_par[v]; }
        return v;
    }

    std::vector<int> dbase_par, dbase_rank;
    int find_dbase(int v) {
        while (dbase_par[v] != v) { dbase_par[v] = dbase_par[dbase_par[v]]; v = dbase_par[v]; }
        return v;
    }
    void union_dbase(int a, int b) {
        a = find_dbase(a); b = find_dbase(b);
        if (a == b) return;
        if (dbase_rank[a] < dbase_rank[b]) std::swap(a, b);
        dbase_par[b] = a;
        if (dbase_rank[a] == dbase_rank[b]) dbase_rank[a]++;
    }
    void make_rep_dbase(int v) {
        int r = find_dbase(v);
        if (r != v) { dbase_par[r] = v; dbase_par[v] = v; }
    }

    int max_pq;
    std::vector<std::vector<int>> L;  /* L[d] = edge IDs becoming tight at Delta=d */

    std::vector<double> path1, path2;
    double strue;

    std::vector<int> tree_nodes;
    int Delta;

    /* ---- H state ---- */
    std::vector<int> rep;
    std::vector<int> mateH;
    std::vector<bool> is_H;           /* per edge ID: is this edge in H? */

    /* ---- Phase 2 ---- */
    std::vector<int> labelH;
    std::vector<int> parentH;         /* parentH[uh] = edge ID */
    std::vector<int> bridgeH;         /* bridgeH[vh] = edge ID */
    std::vector<int> dirH;            /* dirH[vh] = 1 or -1 */
    std::vector<int> even_timeH;
    int tH;
    std::vector<std::vector<int>> contracted_into;
    int size_of_M;

    GabowMCM(int n_, const std::vector<std::pair<int,int>>& edges) : n(n_), m(0), size_of_M(0) {
        adj.resize(n);
        mate.assign(n, NIL);
        label.assign(n, UNLABELED);
        parent.assign(n, NIL);
        source_bridge.assign(n, NIL);
        target_bridge.assign(n, NIL);
        bd.assign(n, 1);
        bDelta.assign(n, 0);
        base_par.resize(n);
        for (int i = 0; i < n; i++) base_par[i] = i;
        dbase_par.resize(n);
        for (int i = 0; i < n; i++) dbase_par[i] = i;
        dbase_rank.assign(n, 0);
        max_pq = n / 2 + 2;
        L.resize(max_pq);
        path1.assign(n, 0.0);
        path2.assign(n, 0.0);
        strue = 0.0;
        Delta = 0;
        rep.resize(n);
        mateH.assign(n, NIL);
        labelH.assign(n, UNLABELED);
        parentH.assign(n, NIL);
        bridgeH.assign(n, NIL);
        dirH.assign(n, 0);
        even_timeH.assign(n, 0);
        tH = 0;
        max_delta_used = 0;
        contracted_into.resize(n);

        /* Build edge list with dedup */
        std::vector<std::pair<int,int>> sorted_edges;
        for (auto& e : edges) {
            int u = e.first, v = e.second;
            if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
                if (u > v) std::swap(u, v);
                sorted_edges.push_back({u, v});
            }
        }
        std::sort(sorted_edges.begin(), sorted_edges.end());
        sorted_edges.erase(std::unique(sorted_edges.begin(), sorted_edges.end()),
                           sorted_edges.end());

        m = (int)sorted_edges.size();
        esrc.resize(m);
        etgt.resize(m);
        for (int i = 0; i < m; i++) {
            esrc[i] = sorted_edges[i].first;
            etgt[i] = sorted_edges[i].second;
            adj[esrc[i]].push_back(i);
            adj[etgt[i]].push_back(i);
        }
        is_H.assign(m, false);
    }

    /* ---- d(v): dual variable ---- */
    int d(int v) {
        int bv = find_base(v);
        if (label[bv] == UNLABELED) return 1;
        if (label[bv] == EVEN) return bd[v] - (Delta - bDelta[v]);
        return bd[v] + (Delta - bDelta[v]);
    }

    /* ---- scan_edge: schedule edge into PQ ---- */
    void scan_edge(int eid, int z) {
        int u = opposite(z, eid);
        if (mate[u] == z || label[find_base(u)] == ODD) return;
        int p = d(z) + d(u);
        int tight_at;
        if (label[find_base(u)] == UNLABELED)
            tight_at = Delta + p;
        else
            tight_at = Delta + p / 2;
        if (tight_at >= 0 && tight_at < max_pq) {
            L[tight_at].push_back(eid);
            if (tight_at > max_delta_used) max_delta_used = tight_at;
        }
    }

    /* ---- shrink_path ---- */
    void shrink_path(int b, int x, int y,
                     std::vector<std::pair<int,int>>& dunions) {
        int v = find_base(x);
        while (v != b) {
            base_par[v] = b;
            dunions.push_back({v, b});
            v = mate[v];
            base_par[v] = b;
            dunions.push_back({v, b});
            base_par[b] = b;
            source_bridge[v] = x;
            target_bridge[v] = y;
            bd[v] = bd[v] + (Delta - bDelta[v]);
            bDelta[v] = Delta;
            for (int eid : adj[v])
                scan_edge(eid, v);
            v = find_base(parent[v]);
        }
        dunions.push_back({b, b});
    }

    int max_delta_used;  /* track highest Delta used for selective L clear */

    /* Previous tree nodes for selective reset */
    std::vector<int> prev_tree_nodes;

    /* Free vertex tracking */
    std::vector<int> free_vertices;
    bool free_list_built = false;

    void build_free_list() {
        free_vertices.clear();
        for (int v = 0; v < n; v++)
            if (mate[v] == NIL) free_vertices.push_back(v);
        free_list_built = true;
    }

    void update_free_list() {
        /* Remove vertices that got matched since last time */
        int j = 0;
        for (int i = 0; i < (int)free_vertices.size(); i++)
            if (mate[free_vertices[i]] == NIL)
                free_vertices[j++] = free_vertices[i];
        free_vertices.resize(j);
    }

    /* ================================================================ */
    /*                          PHASE 1                                 */
    /* ================================================================ */
    bool phase_1() {
        Delta = 0;
        tree_nodes.clear();
        /* Only clear L entries up to max_delta_used from previous phase */
        int clear_limit = std::min(max_delta_used + 1, max_pq);
        for (int i = 0; i < clear_limit; i++) L[i].clear();
        max_delta_used = 0;
        std::vector<std::pair<int,int>> dunions;

        /* Reset only previous tree nodes.
         * Constructor already set correct defaults for all arrays,
         * so first call (prev_tree_nodes empty) needs no full reset. */
        for (int v : prev_tree_nodes) {
            base_par[v] = v;
            dbase_par[v] = v;
            dbase_rank[v] = 0;
            label[v] = UNLABELED;
            parent[v] = NIL;
            source_bridge[v] = NIL;
            target_bridge[v] = NIL;
            bd[v] = 1;
            bDelta[v] = 0;
            for (int eid : adj[v])
                is_H[eid] = false;
        }

        /* Build or update free vertex list */
        if (!free_list_built)
            build_free_list();
        else
            update_free_list();

        /* Label free vertices EVEN, then scan */
        for (int v : free_vertices) {
            label[v] = EVEN;
            tree_nodes.push_back(v);
        }
        for (int v : free_vertices) {
            for (int eid : adj[v])
                scan_edge(eid, v);
        }

        bool found_sap = false;

        while (Delta <= max_delta_used) {
            /* Skip empty levels */
            while (Delta <= max_delta_used && L[Delta].empty()) Delta++;
            if (Delta > max_delta_used) break;

            int qi = 0;
            while (qi < (int)L[Delta].size()) {
                int eid = L[Delta][qi++];
                int x = esrc[eid], y = etgt[eid];

                /* Stale-entry guard: the priority queue schedules edges
                 * based on predicted d-trajectories. A later label change
                 * on either endpoint can invalidate the prediction. Every
                 * label change re-scans and enqueues a fresh correct entry,
                 * so discarding stale entries here loses no correct SAP. */
                if (d(x) + d(y) != w_edge(eid)) continue;

                if (label[find_base(x)] != EVEN) std::swap(x, y);
                if (y == mate[x] || find_base(x) == find_base(y) ||
                    label[find_base(y)] == ODD) continue;

                if (label[find_base(y)] == UNLABELED) {
                    int z = mate[y];
                    bd[y] = 1; bDelta[y] = Delta;
                    bd[z] = 1; bDelta[z] = Delta;
                    parent[z] = y;
                    parent[y] = x;
                    label[y] = ODD;
                    label[z] = EVEN;
                    tree_nodes.push_back(y);
                    tree_nodes.push_back(z);
                    for (int e2 : adj[z])
                        scan_edge(e2, z);

                } else if (label[find_base(y)] == EVEN) {
                    strue += 1.0;
                    int hx = find_base(x), hy = find_base(y);
                    path1[hx] = strue; path2[hy] = strue;
                    int lca = NIL;
                    while (true) {
                        if (path1[hy] == strue) { lca = hy; break; }
                        if (path2[hx] == strue) { lca = hx; break; }
                        bool hxr = (mate[hx] == NIL || parent[mate[hx]] == NIL);
                        bool hyr = (mate[hy] == NIL || parent[mate[hy]] == NIL);
                        if (hxr && hyr) break;
                        if (!hxr) { hx = find_base(parent[mate[hx]]); path1[hx] = strue; }
                        if (!hyr) { hy = find_base(parent[mate[hy]]); path2[hy] = strue; }
                    }
                    if (lca != NIL) {
                        shrink_path(lca, x, y, dunions);
                        shrink_path(lca, y, x, dunions);
                    } else {
                        found_sap = true;
                    }
                }
            }
            L[Delta].clear();

            if (found_sap) {
                /* Build H */
                for (int v : tree_nodes) {
                    contracted_into[find_dbase(v)].push_back(v);
                    mateH[v] = NIL;
                }
                /* Mark tight edges */
                for (int u : tree_nodes) {
                    int uh = find_dbase(u);
                    for (int eid : adj[u]) {
                        int v = opposite(u, eid);
                        int vh = find_dbase(v);
                        if (uh != vh && d(u) + d(v) == w_edge(eid)) {
                            is_H[eid] = true;
                            if (w_edge(eid) == 2) {
                                mateH[uh] = vh;
                                mateH[vh] = uh;
                            }
                        }
                    }
                }
                prev_tree_nodes = tree_nodes;
                return true;
            }

            for (auto& [a, b] : dunions) {
                if (a == b) make_rep_dbase(a);
                else union_dbase(a, b);
            }
            dunions.clear();
            Delta++;
        }
        prev_tree_nodes = tree_nodes;
        return false;
    }

    /* ================================================================ */
    /*                          PHASE 2                                 */
    /* ================================================================ */

    /* find_apHG: recursive DFS on H.
     * Mechanical translation of LEDA. Edge IDs give exact semantics. */
    int find_apHG(int vh) {
        for (int v : contracted_into[vh]) {
            for (int eid : adj[v]) {
                if (!is_H[eid]) continue;
                int w = opposite(v, eid);
                /* LEDA: uh = rep[G.opposite(v, eh)] */
                int uh = rep[w];
                if (mateH[vh] == uh) continue;

                if (labelH[uh] == UNLABELED) {
                    int muh = mateH[uh];
                    if (muh == NIL) {
                        labelH[uh] = ODD;
                        parentH[uh] = eid;
                        return uh;
                    }
                    labelH[uh] = ODD;
                    labelH[muh] = EVEN;
                    parentH[uh] = eid;
                    even_timeH[muh] = tH++;
                    int s = find_apHG(muh);
                    if (s != NIL) return s;

                } else {
                    int bh = find_dbase(vh);
                    int zh = find_dbase(uh);
                    if (even_timeH[bh] < even_timeH[zh]) {
                        std::vector<int> tmp;
                        std::vector<int> endpoints;
                        while (zh != bh) {
                            endpoints.push_back(zh);
                            zh = mateH[zh];
                            endpoints.push_back(zh);
                            tmp.insert(tmp.begin(), zh);
                            /* LEDA: zh = dbase(rep[rep[G.source(parentHG[zh])] == zh ?
                             *   G.target(parentHG[zh]) : G.source(parentHG[zh])]) */
                            int pe = parentH[zh];
                            zh = find_dbase(rep[rep[esrc[pe]] == zh ? etgt[pe] : esrc[pe]]);
                        }
                        for (int nd : endpoints) union_dbase(nd, bh);
                        make_rep_dbase(bh);
                        for (int odd_node : tmp) {
                            bridgeH[odd_node] = eid;
                            /* LEDA: dirHG[zh] = (G.target(eh) == v ? 1 : -1) */
                            dirH[odd_node] = (etgt[eid] == v) ? 1 : -1;
                        }
                        for (int odd_node : tmp) {
                            int s = find_apHG(odd_node);
                            if (s != NIL) return s;
                        }
                    }
                }
            }
        }
        return NIL;
    }

    /* find_path_in_HG: trace augmenting path in H from vh to uh.
     * Exact LEDA translation using edge IDs. */
    void find_path_in_HG(std::vector<int>& path, int vh, int uh) {
        if (vh == uh) return;
        if (labelH[vh] == EVEN) {
            int mvh = mateH[vh];
            int pe = parentH[mvh];
            path.push_back(pe);
            /* LEDA: rep[rep[G.source(e)] == mvh ? G.target(e) : G.source(e)] */
            int next = rep[rep[esrc[pe]] == mvh ? etgt[pe] : esrc[pe]];
            find_path_in_HG(path, next, uh);
        } else {
            /* ODD: use bridge */
            int be = bridgeH[vh];
            /* LEDA: dir==1 ? source : target goes to mate side
             *        dir==1 ? target : source goes to uh side */
            int mate_side, uh_side;
            if (dirH[vh] == 1) {
                mate_side = rep[esrc[be]];
                uh_side = rep[etgt[be]];
            } else {
                mate_side = rep[etgt[be]];
                uh_side = rep[esrc[be]];
            }
            int mt = (mateH[vh] != NIL) ? rep[mateH[vh]] : vh;
            find_path_in_HG(path, mate_side, mt);
            path.push_back(be);
            find_path_in_HG(path, uh_side, uh);
        }
    }

    /* find_path_in_G: unfold within Phase 1 blossom */
    void find_path_in_G(std::vector<std::pair<int,int>>& pairs, int v, int u) {
        if (v == u) return;
        if (label[v] == EVEN) {
            pairs.push_back({mate[v], parent[mate[v]]});
            find_path_in_G(pairs, parent[mate[v]], u);
        } else {
            find_path_in_G(pairs, source_bridge[v], mate[v]);
            pairs.push_back({source_bridge[v], target_bridge[v]});
            find_path_in_G(pairs, target_bridge[v], u);
        }
    }


    /* augmentG: unfold H-path edges to G, augment matching */
    void augmentG(const std::vector<int>& h_edge_ids) {
        std::vector<std::pair<int,int>> pairs;
        for (int eid : h_edge_ids) {
            int u = esrc[eid], v = etgt[eid];
            pairs.push_back({u, v});
            find_path_in_G(pairs, u, rep[u]);
            find_path_in_G(pairs, v, rep[v]);
        }
        for (auto& [a, b] : pairs) {
            mate[a] = b; mate[b] = a;
        }
        size_of_M++;
    }

    void phase_2() {
        tH = 0;
        /* DO NOT reset dbase here. LEDA keeps Phase 1 dbase state
         * and adds Phase 2 unions on top. The dbase.split(T) only
         * happens at the start of the NEXT Phase 1. */
        for (int v : tree_nodes) {
            rep[v] = find_dbase(v);
            labelH[v] = UNLABELED;
            parentH[v] = NIL;
            bridgeH[v] = NIL;
            dirH[v] = 0;
            even_timeH[v] = 0;
        }

        std::vector<std::vector<int>> all_paths;

        for (int vh : tree_nodes) {
            if (vh != rep[vh]) continue;
            if (labelH[vh] == UNLABELED && mateH[vh] == NIL) {
                labelH[vh] = EVEN;
                even_timeH[vh] = tH++;
                int found = find_apHG(vh);
                if (found != NIL) {
                    std::vector<int> path;
                    int pe = parentH[found];
                    path.push_back(pe);
                    int next = rep[rep[esrc[pe]] == found ? etgt[pe] : esrc[pe]];
                    find_path_in_HG(path, next, vh);
                    all_paths.push_back(std::move(path));
                }
            }
        }

        for (auto& p : all_paths)
            augmentG(p);


        for (int v : tree_nodes) {
            contracted_into[v].clear();
            mateH[v] = NIL;
        }
    }

    /* ================================================================ */
    int greedy_init() {
        int cnt = 0;
        for (int u = 0; u < n; u++) {
            if (mate[u] != NIL) continue;
            for (int eid : adj[u]) {
                int v = opposite(u, eid);
                if (mate[v] == NIL) { mate[u] = v; mate[v] = u; cnt++; break; }
            }
        }
        return cnt;
    }

    int greedy_init_md() {
        int cnt = 0;
        std::vector<int> deg(n, 0);
        for (int i = 0; i < m; i++) { deg[esrc[i]]++; deg[etgt[i]]++; }
        std::vector<int> order(n);
        for (int i = 0; i < n; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b){
            return deg[a] < deg[b] || (deg[a] == deg[b] && a < b);
        });
        for (int u : order) {
            if (mate[u] != NIL) continue;
            int best = -1, best_deg = INT_MAX;
            for (int eid : adj[u]) {
                int v = opposite(u, eid);
                if (mate[v] == NIL && deg[v] < best_deg) { best = v; best_deg = deg[v]; }
            }
            if (best >= 0) { mate[u] = best; mate[best] = u; cnt++; }
        }
        return cnt;
    }

    std::vector<std::pair<int,int>> solve(int greedy_mode = 0) {
        if (greedy_mode == 1) { greedy_size = greedy_init(); size_of_M = greedy_size; }
        else if (greedy_mode == 2) { greedy_size = greedy_init_md(); size_of_M = greedy_size; }

        int phase_count = 0;
        while (true) {
            bool has_sap = phase_1();
            if (!has_sap) break;
            phase_2();
            phase_count++;
        }
        printf("Phases: %d\n", phase_count);

        std::vector<std::pair<int,int>> result;
        for (int u = 0; u < n; u++)
            if (mate[u] != NIL && mate[u] > u)
                result.push_back({u, mate[u]});
        std::sort(result.begin(), result.end());
        return result;
    }
};

/* ================================================================ */
void validate_matching(int n, const std::vector<std::vector<int>>& adj,
                       const std::vector<int>& esrc, const std::vector<int>& etgt,
                       const std::vector<std::pair<int,int>>& matching) {
    std::vector<int> deg(n, 0);
    int errors = 0;
    for (auto& [u, v] : matching) {
        bool found = false;
        for (int eid : adj[u]) {
            if ((esrc[eid] == u && etgt[eid] == v) || (esrc[eid] == v && etgt[eid] == u))
                { found = true; break; }
        }
        if (!found) { fprintf(stderr, "ERROR: Edge (%d,%d) not in graph!\n", u, v); errors++; }
        deg[u]++; deg[v]++;
    }
    for (int i = 0; i < n; i++)
        if (deg[i] > 1) { fprintf(stderr, "ERROR: Vertex %d in %d edges!\n", i, deg[i]); errors++; }
    printf("\n=== Validation Report ===\n");
    printf("Matching size: %d\n", (int)matching.size());
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

int main(int argc, char* argv[]) {
    printf("Gabow MCM (duals + edge IDs) - C++\n");
    printf("===================================\n\n");
    if (argc < 2) { printf("Usage: %s <filename> [--greedy|--greedy-md]\n", argv[0]); return 1; }
    int greedy_mode = 0;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedy_mode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedy_mode = 2;
    }
    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    int n, m_in;
    if (fscanf(f, "%d %d", &n, &m_in) != 2) { fprintf(stderr, "Bad header\n"); fclose(f); return 1; }
    std::vector<std::pair<int,int>> edges;
    edges.reserve(m_in);
    for (int i = 0; i < m_in; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edges.push_back({u, v});
    }
    fclose(f);
    printf("Graph: %d vertices, %d edges (input)\n", n, (int)edges.size());
    auto t0 = std::chrono::high_resolution_clock::now();
    GabowMCM gabow(n, edges);
    printf("Graph: %d vertices, %d edges (deduped)\n", gabow.n, gabow.m);
    auto matching = gabow.solve(greedy_mode);
    auto t1 = std::chrono::high_resolution_clock::now();
    validate_matching(n, gabow.adj, gabow.esrc, gabow.etgt, matching);
    printf("Matching size: %d\n", (int)matching.size());
    if (greedy_mode > 0) {
        printf("Greedy init size: %d\n", gabow.greedy_size);
        if ((int)matching.size() > 0)
            printf("Greedy/Final: %.2f%%\n", 100.0 * gabow.greedy_size / matching.size());
    }
    printf("Time: %ld ms\n",
        (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
    return 0;
}

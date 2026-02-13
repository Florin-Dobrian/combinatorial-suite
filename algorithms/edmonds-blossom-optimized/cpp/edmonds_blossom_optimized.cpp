/*
 * Edmonds' Blossom Algorithm (Optimized) â€” Unweighted Maximum Cardinality Matching
 *
 * Forest BFS: each stage labels ALL free vertices as S-roots simultaneously
 * and grows a search forest. An augmenting path is found when two different
 * trees meet (S-S edge across trees). One augmentation per stage, then
 * expand all blossoms and repeat until no augmenting path exists.
 *
 * Same blossom machinery as edmonds-simple (NetworkX-derived), just with
 * forest search instead of single-source tree search.
 *
 * Blossom IDs reset to n each stage. All indices are 32-bit signed int.
 *
 * Complexity: O(V * E) worst case (one stage per augmentation, each stage O(E)).
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <climits>
#include <chrono>
#include <cassert>

struct Solver {
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<int> mate; // mate[v] = matched partner, or -1

    // Blossom storage. IDs 0..n-1 are trivial (one vertex each, no data).
    // Non-trivial blossoms have id in [n, nblos). Reset each BFS.
    struct Blos {
        std::vector<int> childs;               // sub-blossom IDs in cycle order
        std::vector<std::pair<int,int>> edges;  // edges[i] connects childs[i] to childs[(i+1)%k]
    };
    std::vector<Blos> blos;
    int nblos;                     // next blossom ID to allocate

    std::vector<int> inblossom;    // inblossom[v] = top-level blossom containing v
    std::vector<int> blossomparent; // blossomparent[b] = parent blossom, or -1
    std::vector<int> blossombase;  // blossombase[b] = base vertex of blossom b

    // Per-search state (sized to nblos, reset each BFS)
    std::vector<int> label;                    // 0=unlabeled, 1=S, 2=T (5=breadcrumb)
    std::vector<std::pair<int,int>> labeledge; // label edge for tree structure
    std::vector<int> queue;                    // BFS queue of S-vertices

    Solver(int nn, const std::vector<std::pair<int,int>>& edges_in) : n(nn) {
        adj.resize(n);
        for (auto& e : edges_in) {
            int u = e.first, v = e.second;
            if (u != v && u >= 0 && u < n && v >= 0 && v < n) {
                adj[u].push_back(v);
                adj[v].push_back(u);
            }
        }
        for (int i = 0; i < n; i++) {
            std::sort(adj[i].begin(), adj[i].end());
            adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());
        }
        mate.assign(n, -1);
        inblossom.resize(n);
        blossombase.resize(n);
        blossomparent.resize(n, -1);
        for (int i = 0; i < n; i++) { inblossom[i] = i; blossombase[i] = i; }
        nblos = n;
    }

    // Ensure arrays cover blossom id b
    void ensure(int b) {
        if (b < (int)label.size()) return;
        int old = (int)label.size();
        label.resize(b + 1, 0);
        labeledge.resize(b + 1, {-1, -1});
        blossomparent.resize(b + 1, -1);
        blossombase.resize(b + 1, -1);
    }

    bool isBlossom(int b) const { return b >= n; }

    void leaves(int b, std::vector<int>& out) {
        if (!isBlossom(b)) { out.push_back(b); return; }
        for (int c : blos[b].childs) leaves(c, out);
    }

    // Reset blossom state for a new BFS. All blossoms from the previous search
    // have already been expanded, so every vertex is its own top-level blossom.
    void resetBlossoms() {
        nblos = n;
        blos.resize(n); // shrink away any leftover non-trivial blossom data
        for (int i = 0; i < n; i++) {
            inblossom[i] = i;
            blossombase[i] = i;
            blossomparent[i] = -1;
        }
        label.assign(n, 0);
        labeledge.assign(n, {-1, -1});
        queue.clear();
    }

    // ---- Tree building ----

    void assignLabel(int w, int t, int v) {
        int b = inblossom[w];
        ensure(b);
        label[b] = t;
        label[w] = t;
        if (v != -1) {
            labeledge[w] = labeledge[b] = {v, w};
        } else {
            labeledge[w] = labeledge[b] = {-1, -1};
        }
        if (t == 1) {
            // S-blossom: add its leaves to the BFS queue
            std::vector<int> lv;
            leaves(b, lv);
            for (int u : lv) queue.push_back(u);
        } else if (t == 2) {
            // T-blossom: label the mate of its base as S
            int base = blossombase[b];
            assignLabel(mate[base], 1, base);
        }
    }

    // ---- Blossom detection ----

    // Trace from two S-vertices to find their LCA (blossom base).
    // Returns base vertex, or -2 if they belong to different trees
    // (augmenting path â€” should not occur in single-source mode).
    int scanBlossom(int v, int w) {
        std::vector<int> path;
        int base = -2;
        while (v != -2 || w != -2) {
            if (v != -2) {
                int b = inblossom[v];
                if (label[b] & 4) { base = blossombase[b]; break; }
                path.push_back(b);
                label[b] = 5; // breadcrumb
                auto& le = labeledge[b];
                if (le.first == -1) {
                    v = -2; // reached root
                } else {
                    v = le.first;
                    int bt = inblossom[v];
                    v = labeledge[bt].first;
                }
                if (w != -2) std::swap(v, w);
            } else {
                std::swap(v, w);
            }
        }
        for (int b : path) label[b] = 1; // restore breadcrumbs
        return base;
    }

    // ---- Blossom contraction ----

    void addBlossom(int base, int v, int w) {
        int bb = inblossom[base];
        int bv = inblossom[v];
        int bw = inblossom[w];

        int bid = nblos++;
        if (bid >= (int)blos.size()) blos.push_back(Blos());
        else { blos[bid].childs.clear(); blos[bid].edges.clear(); }
        ensure(bid);
        blossombase[bid] = base;
        blossomparent[bid] = -1;
        blossomparent[bb] = bid;

        auto& childs = blos[bid].childs;
        auto& edges = blos[bid].edges;
        edges.push_back({v, w}); // bridge edge

        // Trace from v back to base
        while (bv != bb) {
            blossomparent[bv] = bid;
            childs.push_back(bv);
            edges.push_back(labeledge[bv]);
            v = labeledge[bv].first;
            bv = inblossom[v];
        }
        childs.push_back(bb);
        std::reverse(childs.begin(), childs.end());
        std::reverse(edges.begin(), edges.end());

        // Trace from w back to base
        while (bw != bb) {
            blossomparent[bw] = bid;
            childs.push_back(bw);
            auto le = labeledge[bw];
            edges.push_back({le.second, le.first}); // reversed
            w = labeledge[bw].first;
            bw = inblossom[w];
        }

        label[bid] = 1;
        labeledge[bid] = labeledge[bb];

        // Relabel: T-vertices inside the blossom become S
        std::vector<int> lv;
        leaves(bid, lv);
        for (int u : lv) {
            if (label[inblossom[u]] == 2) queue.push_back(u);
            inblossom[u] = bid;
        }
    }

    // ---- Blossom expansion ----

    void expandBlossom(int b, bool endstage) {
        struct Frame { int b; bool endstage; int idx; };
        std::vector<Frame> stack;
        stack.push_back({b, endstage, 0});

        while (!stack.empty()) {
            auto& f = stack.back();
            auto& bl = blos[f.b];
            if (f.idx < (int)bl.childs.size()) {
                int s = bl.childs[f.idx];
                f.idx++;
                blossomparent[s] = -1;
                if (isBlossom(s)) {
                    if (f.endstage) {
                        // Recursively expand sub-blossoms at end of stage
                        stack.push_back({s, true, 0});
                        continue;
                    } else {
                        std::vector<int> lv;
                        leaves(s, lv);
                        for (int u : lv) inblossom[u] = s;
                    }
                } else {
                    inblossom[s] = s;
                }
            } else {
                // All children processed
                if (!f.endstage && label[f.b] == 2) {
                    // Mid-stage T-blossom expansion: relabel children
                    auto& bl2 = blos[f.b];
                    int entrychild = inblossom[labeledge[f.b].second];
                    int k = (int)bl2.childs.size();
                    int j = 0;
                    for (; j < k; j++) if (bl2.childs[j] == entrychild) break;
                    int jstep;
                    if (j & 1) { j -= k; jstep = 1; } else { jstep = -1; }
                    int lv_ = labeledge[f.b].first, lw_ = labeledge[f.b].second;
                    while (j != 0) {
                        int pp, qq;
                        if (jstep == 1) {
                            pp = bl2.edges[((j % k) + k) % k].first;
                            qq = bl2.edges[((j % k) + k) % k].second;
                        } else {
                            int ei = (((j - 1) % k) + k) % k;
                            qq = bl2.edges[ei].first;
                            pp = bl2.edges[ei].second;
                        }
                        label[lw_] = 0;
                        label[qq] = 0;
                        assignLabel(lw_, 2, lv_);
                        j += jstep;
                        if (jstep == 1) {
                            lv_ = bl2.edges[((j % k) + k) % k].first;
                            lw_ = bl2.edges[((j % k) + k) % k].second;
                        } else {
                            int ei = (((j - 1) % k) + k) % k;
                            lw_ = bl2.edges[ei].first;
                            lv_ = bl2.edges[ei].second;
                        }
                        j += jstep;
                    }
                    int bwi = bl2.childs[((j % k) + k) % k];
                    ensure(bwi);
                    label[lw_] = label[bwi] = 2;
                    labeledge[lw_] = labeledge[bwi] = {lv_, lw_};
                    j += jstep;
                    while (bl2.childs[((j % k) + k) % k] != entrychild) {
                        int bvi = bl2.childs[((j % k) + k) % k];
                        ensure(bvi);
                        if (label[bvi] == 1) { j += jstep; continue; }
                        int found_v = -1;
                        if (isBlossom(bvi)) {
                            std::vector<int> lvs;
                            leaves(bvi, lvs);
                            for (int u : lvs) if (label[u]) { found_v = u; break; }
                        } else {
                            found_v = bvi;
                        }
                        if (found_v != -1 && label[found_v]) {
                            label[found_v] = 0;
                            label[mate[blossombase[bvi]]] = 0;
                            assignLabel(found_v, 2, labeledge[found_v].first);
                        }
                        j += jstep;
                    }
                }
                label[f.b] = 0;
                bl.childs.clear();
                bl.edges.clear();
                stack.pop_back();
            }
        }
    }

    // ---- Augmentation through blossoms ----

    void augmentBlossom(int b, int v) {
        struct Frame { int b; int v; int phase; int i; int j; int jstep; };
        std::vector<Frame> stack;
        stack.push_back({b, v, 0, 0, 0, 0});

        while (!stack.empty()) {
            auto& f = stack.back();
            if (f.phase == 0) {
                // Find sub-blossom containing v
                int t = f.v;
                while (blossomparent[t] != f.b) t = blossomparent[t];
                auto& bl = blos[f.b];
                int k = (int)bl.childs.size();
                f.i = 0;
                for (; f.i < k; f.i++) if (bl.childs[f.i] == t) break;
                if (isBlossom(t)) {
                    f.phase = 1;
                    stack.push_back({t, f.v, 0, 0, 0, 0});
                    continue;
                }
                f.phase = 2;
                if (f.i & 1) { f.j = f.i - k; f.jstep = 1; }
                else          { f.j = f.i;     f.jstep = -1; }
                continue;
            }
            if (f.phase == 1) {
                // After recursion into sub-blossom
                f.phase = 2;
                int k = (int)blos[f.b].childs.size();
                if (f.i & 1) { f.j = f.i - k; f.jstep = 1; }
                else          { f.j = f.i;     f.jstep = -1; }
                continue;
            }
            if (f.phase == 2) {
                // Main loop: walk from position i toward position 0
                auto& bl = blos[f.b];
                int k = (int)bl.childs.size();
                if (f.j == 0) {
                    // Done: rotate childs/edges so new base is first
                    if (f.i > 0) {
                        std::vector<int> nc(bl.childs.begin() + f.i, bl.childs.end());
                        nc.insert(nc.end(), bl.childs.begin(), bl.childs.begin() + f.i);
                        std::vector<std::pair<int,int>> ne(bl.edges.begin() + f.i, bl.edges.end());
                        ne.insert(ne.end(), bl.edges.begin(), bl.edges.begin() + f.i);
                        bl.childs = nc;
                        bl.edges = ne;
                    }
                    blossombase[f.b] = f.v;
                    stack.pop_back();
                    continue;
                }
                // Step to next pair of sub-blossoms
                f.j += f.jstep;
                int idx1 = ((f.j % k) + k) % k;
                int c1 = bl.childs[idx1];
                int ww, xx;
                if (f.jstep == 1) {
                    ww = bl.edges[idx1].first;
                    xx = bl.edges[idx1].second;
                } else {
                    int ei = (((f.j - 1) % k) + k) % k;
                    xx = bl.edges[ei].first;
                    ww = bl.edges[ei].second;
                }
                if (isBlossom(c1)) {
                    f.phase = 3;
                    stack.push_back({c1, ww, 0, 0, 0, 0});
                    continue;
                }
                f.phase = 3;
            }
            if (f.phase == 3) {
                // After optional recursion for c1, step to c2
                auto& bl = blos[f.b];
                int k = (int)bl.childs.size();
                int idx1 = ((f.j % k) + k) % k;
                int ww, xx;
                if (f.jstep == 1) {
                    ww = bl.edges[idx1].first;
                    xx = bl.edges[idx1].second;
                } else {
                    int ei = (((f.j - 1) % k) + k) % k;
                    xx = bl.edges[ei].first;
                    ww = bl.edges[ei].second;
                }
                f.j += f.jstep;
                int idx2 = ((f.j % k) + k) % k;
                int c2 = bl.childs[idx2];
                if (isBlossom(c2)) {
                    f.phase = 4;
                    stack.push_back({c2, xx, 0, 0, 0, 0});
                    continue;
                }
                f.phase = 4;
            }
            if (f.phase == 4) {
                // After optional recursion for c2, set mate pair
                auto& bl = blos[f.b];
                int k = (int)bl.childs.size();
                int prev_j = f.j - f.jstep;
                int idx1 = ((prev_j % k) + k) % k;
                int ww, xx;
                if (f.jstep == 1) {
                    ww = bl.edges[idx1].first;
                    xx = bl.edges[idx1].second;
                } else {
                    int ei = (((prev_j - 1) % k) + k) % k;
                    xx = bl.edges[ei].first;
                    ww = bl.edges[ei].second;
                }
                mate[ww] = xx;
                mate[xx] = ww;
                f.phase = 2; // continue loop
            }
        }
    }

    // ---- Augmenting path: trace both sides back to their roots ----

    void augmentMatching(int v, int w) {
        // v and w are S-vertices in different trees. Edge (v,w) completes
        // an augmenting path. Trace from each side back to its root,
        // flipping matched/unmatched edges and augmenting through blossoms.
        for (int iter = 0; iter < 2; iter++) {
            int s = (iter == 0) ? v : w;
            int j = (iter == 0) ? w : v;
            while (true) {
                int bs = inblossom[s];
                if (isBlossom(bs)) augmentBlossom(bs, s);
                mate[s] = j;
                auto& le = labeledge[bs];
                if (le.first == -1) break; // root
                int t = le.first;          // T-vertex
                int bt = inblossom[t];
                auto& le2 = labeledge[bt];
                s = le2.first;
                j = le2.second;
                if (isBlossom(bt)) augmentBlossom(bt, j);
                mate[j] = s;
            }
        }
    }

    // ---- Greedy initialization ----

    int greedy_size = 0;

    int greedy_init() {
        int cnt = 0;
        for (int u = 0; u < n; u++) {
            if (mate[u] != -1) continue;
            for (int v : adj[u]) {
                if (mate[v] == -1) { mate[u] = v; mate[v] = u; cnt++; break; }
            }
        }
        return cnt;
    }

    int greedy_init_md() {
        int cnt = 0;
        std::vector<int> deg(n, 0);
        for (int u = 0; u < n; u++) for (int v : adj[u]) deg[v]++;
        std::vector<int> order(n);
        for (int i = 0; i < n; i++) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return deg[a] < deg[b] || (deg[a] == deg[b] && a < b); });
        for (int u : order) {
            if (mate[u] != -1) continue;
            int best = -1, bd = INT_MAX;
            for (int v : adj[u])
                if (mate[v] == -1 && deg[v] < bd) { best = v; bd = deg[v]; }
            if (best >= 0) { mate[u] = best; mate[best] = u; cnt++; }
        }
        return cnt;
    }

    // ---- Main solver ----

    std::vector<std::pair<int,int>> solve(int greedy_mode = 0) {
        if (greedy_mode == 1) greedy_size = greedy_init();
        else if (greedy_mode == 2) greedy_size = greedy_init_md();

        while (true) {
            // New stage: reset all blossom state
            resetBlossoms();

            // Label ALL free vertices as S-roots
            for (int v = 0; v < n; v++) {
                if (mate[v] == -1 && label[inblossom[v]] == 0) {
                    assignLabel(v, 1, -1);
                }
            }

            // BFS: grow forest until augmenting path or exhaustion
            bool augmented = false;
            while (!queue.empty() && !augmented) {
                int v = queue.back(); queue.pop_back();
                if (label[inblossom[v]] != 1) continue; // stale
                for (int w : adj[v]) {
                    int bv = inblossom[v];
                    int bw = inblossom[w];
                    if (bv == bw) continue;
                    ensure(bw);
                    if (label[bw] == 0) {
                        // w is unlabeled: grow the tree
                        assignLabel(w, 2, v);
                    } else if (label[bw] == 1) {
                        // S-S edge: blossom or augmenting path
                        int base = scanBlossom(v, w);
                        if (base >= 0) {
                            addBlossom(base, v, w);
                        } else {
                            // base == -2: two different trees met â†’ augmenting path
                            augmentMatching(v, w);
                            augmented = true;
                            break;
                        }
                    }
                    // label[bw]==2: T-blossom edge, ignore
                }
            }

            // Expand all remaining blossoms (end of stage)
            for (int b = n; b < nblos; b++) {
                if (!blos[b].childs.empty() && blossomparent[b] == -1) {
                    expandBlossom(b, true);
                }
            }

            if (!augmented) break; // no augmenting path found: matching is maximum
        }

        std::vector<std::pair<int,int>> result;
        for (int u = 0; u < n; u++)
            if (mate[u] > u) result.push_back({u, mate[u]});
        std::sort(result.begin(), result.end());
        return result;
    }
};

// ---- Validation and main ----

void validate_matching(int n, const std::vector<std::vector<int>>& graph,
                       const std::vector<std::pair<int,int>>& matching) {
    std::vector<int> deg(n, 0);
    int errors = 0;
    for (auto& e : matching) {
        if (!std::binary_search(graph[e.first].begin(), graph[e.first].end(), e.second)) {
            fprintf(stderr, "ERROR: Edge (%d,%d) not in graph!\n", e.first, e.second);
            errors++;
        }
        deg[e.first]++;
        deg[e.second]++;
    }
    for (int i = 0; i < n; i++)
        if (deg[i] > 1) {
            fprintf(stderr, "ERROR: Vertex %d in %d edges!\n", i, deg[i]);
            errors++;
        }
    int mt = 0;
    for (int i = 0; i < n; i++) if (deg[i] > 0) mt++;
    printf("\n=== Validation Report ===\n"
           "Matching size: %d\nMatched vertices: %d\n%s\n"
           "=========================\n\n",
           (int)matching.size(), mt,
           errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
}

int main(int argc, char* argv[]) {
    printf("Edmonds' Blossom Algorithm (Simple) - C++ Implementation\n"
           "=========================================================\n\n");
    if (argc < 2) {
        printf("Usage: %s <filename> [--greedy|--greedy-md]\n", argv[0]);
        return 1;
    }
    int gm = 0;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") gm = 1;
        else if (std::string(argv[i]) == "--greedy-md") gm = 2;
    }
    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "Cannot open: %s\n", argv[1]); return 1; }
    int nn, m;
    if (fscanf(f, "%d %d", &nn, &m) != 2) { fclose(f); return 1; }
    std::vector<std::pair<int,int>> edges;
    edges.reserve(m);
    for (int i = 0; i < m; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edges.push_back({u, v});
    }
    fclose(f);
    printf("Graph: %d vertices, %d edges\n", nn, (int)edges.size());
    auto t0 = std::chrono::high_resolution_clock::now();
    Solver sol(nn, edges);
    auto matching = sol.solve(gm);
    auto t1 = std::chrono::high_resolution_clock::now();
    validate_matching(nn, sol.adj, matching);
    printf("Matching size: %d\n", (int)matching.size());
    if (gm > 0) {
        printf("Greedy init size: %d\n", sol.greedy_size);
        if (!matching.empty())
            printf("Greedy/Final: %.2f%%\n", 100.0 * sol.greedy_size / matching.size());
    }
    printf("Time: %ld ms\n",
           (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
    return 0;
}

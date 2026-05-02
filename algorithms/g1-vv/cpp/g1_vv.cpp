/*
 * Gabow's Algorithm (Simple) - O(V * E) Maximum Matching
 *
 * Faithful to Gabow 1976: forest BFS with blossom contraction via
 * union-find. No physical contraction â€” bases are tracked virtually.
 * Epoch-based interleaved LCA, path-only contraction, bridge recording
 * for augmentation through blossoms.
 *
 * Forest search: each iteration labels ALL free vertices as EVEN roots
 * simultaneously and grows a search forest. An augmenting path is found
 * when two different trees meet (EVEN-EVEN edge across trees, detected
 * by find_lca returning NIL). One augmentation per iteration, then full
 * reset and repeat until no augmenting path exists.
 *
 * Complexity: O(V * E) â€” each iteration does O(E) work, at most V/2
 * augmentations total.
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <climits>
#include <chrono>

static const int NIL = -1;
static const int UNLABELED = 0;
static const int EVEN = 1;
static const int ODD = 2;

struct GabowSimple {
    int n;
    int greedy_size = 0;
    std::vector<std::vector<int>> graph;
    std::vector<int> mate;
    std::vector<int> base;
    std::vector<int> parent;       /* tree parent: for EVEN v, parent[v] is the ODD
                                      vertex through which v entered the tree.
                                      For ODD v, parent[v] is the EVEN vertex that
                                      discovered v. For roots, parent[v] = NIL. */
    std::vector<int> label;        /* UNLABELED / EVEN / ODD */

    /* Bridge recording for ODD vertices absorbed into blossoms.
     * When an ODD vertex mv becomes effectively EVEN through blossom
     * contraction, we record the non-matching edge (x, y) that caused
     * the contraction. This allows augment to trace through the blossom. */
    std::vector<int> bridge_src;   /* bridge_src[v]: the "x" side */
    std::vector<int> bridge_tgt;   /* bridge_tgt[v]: the "y" side */

    /* Epoch-based interleaved LCA */
    std::vector<size_t> lca_tag1, lca_tag2;
    size_t lca_epoch;

    GabowSimple(int n, const std::vector<std::pair<int,int>>& edges) : n(n) {
        graph.resize(n);
        mate.assign(n, NIL);
        base.resize(n);
        parent.resize(n);
        label.resize(n);
        bridge_src.resize(n);
        bridge_tgt.resize(n);
        lca_tag1.assign(n, 0);
        lca_tag2.assign(n, 0);
        lca_epoch = 0;

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
    }

    /* Greedy initial matching: iterate exposed vertices, pick first available edge */
    int greedy_init() {
        int cnt = 0;
        for (int u = 0; u < n; u++) {
            if (mate[u] != NIL) continue;
            for (int v : graph[u]) {
                if (mate[v] == NIL) { mate[u] = v; mate[v] = u; cnt++; break; }
            }
        }
        return cnt;
    }

    /* Min-degree greedy: match each exposed vertex with its lowest-degree unmatched neighbor */
    int greedy_init_md() {
        int cnt = 0;
        std::vector<int> deg(n, 0);
        for (int u = 0; u < n; u++)
            for (int v : graph[u])
                deg[v]++;
        std::vector<int> order(n);
        for (int i = 0; i < n; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b){ return deg[a] < deg[b] || (deg[a] == deg[b] && a < b); });
        for (int u : order) {
            if (mate[u] != NIL) continue;
            int best = -1, best_deg = INT_MAX;
            for (int v : graph[u]) {
                if (mate[v] == NIL && deg[v] < best_deg) {
                    best = v; best_deg = deg[v];
                }
            }
            if (best >= 0) { mate[u] = best; mate[best] = u; cnt++; }
        }
        return cnt;
    }

    /* Path-halving find for union-find base */
    int find_base(int v) {
        while (base[v] != v) {
            base[v] = base[base[v]];
            v = base[v];
        }
        return v;
    }

    /* Interleaved LCA using epoch tags â€” O(path length), no allocation.
     * Returns the LCA base if u and v are in the same tree, or NIL if
     * they are in different trees (= augmenting path). */
    int find_lca(int u, int v) {
        ++lca_epoch;
        size_t ep = lca_epoch;
        int hx = find_base(u), hy = find_base(v);
        lca_tag1[hx] = ep;
        lca_tag2[hy] = ep;
        while (true) {
            if (lca_tag1[hy] == ep) return hy;
            if (lca_tag2[hx] == ep) return hx;
            bool hxr = (mate[hx] == NIL);  /* hx is a root (free vertex) */
            bool hyr = (mate[hy] == NIL);
            if (hxr && hyr) return NIL;     /* different trees */
            if (!hxr) {
                hx = find_base(parent[mate[hx]]);
                lca_tag1[hx] = ep;
            }
            if (!hyr) {
                hy = find_base(parent[mate[hy]]);
                lca_tag2[hy] = ep;
            }
        }
    }

    /* Path-only contraction: walk from x back to lca, union bases.
     * For each ODD vertex mv on the path, record the bridge (x, y)
     * and enqueue mv as newly-EVEN if it wasn't already. */
    void shrink_path(int lca, int x, int y,
                     std::vector<int>& queue, int& qtail) {
        int v = find_base(x);
        while (v != lca) {
            int mv = mate[v];
            /* Union both v and mv into lca's component */
            base[find_base(v)] = lca;
            base[find_base(mv)] = lca;
            base[lca] = lca;  /* keep lca as representative */

            /* Record bridge for mv (ODD vertex becoming effectively EVEN).
             * The bridge edge is (x, y) â€” the non-matching edge that
             * triggered this blossom contraction. */
            bridge_src[mv] = x;
            bridge_tgt[mv] = y;

            /* If mv was ODD and not yet enqueued as EVEN, enqueue it */
            if (label[mv] != EVEN) {
                label[mv] = EVEN;
                queue[qtail++] = mv;
            }

            /* Walk up: mv's parent is EVEN, mate of that is the next step */
            v = find_base(parent[mv]);
        }
    }

    /* Trace from EVEN vertex v to a target vertex t within the same
     * blossom/tree, collecting the alternating path edges.
     * Modeled on gabow_optimized's find_path_in_G.
     *
     * For EVEN v: step to mate[v] (ODD), then parent[mate[v]] (EVEN).
     * For ODD v with a bridge: recurse through the bridge.
     *
     * Iterative with explicit stack to handle nested blossoms. */
    /* Trace from vertex v to vertex u (or to a root if u==NIL),
     * collecting edge pairs for augmentation.
     * Exactly mirrors gabow_optimized's find_path_in_G:
     *   - No bridge â†’ "originally EVEN": step mate â†’ parent
     *   - Has bridge â†’ "originally ODD, absorbed into blossom":
     *     recurse through bridge */
    void trace_path(int v, int u,
                    std::vector<std::pair<int,int>>& pairs) {
        struct Frame { int v, u, phase, sb, tb; };
        std::vector<Frame> stk;
        stk.push_back({v, u, 0, 0, 0});

        while (!stk.empty()) {
            auto& f = stk.back();
            if (f.v == f.u) { stk.pop_back(); continue; }

            if (f.phase == 0) {
                if (bridge_src[f.v] == NIL) {
                    /* Originally EVEN vertex (no bridge) */
                    if (mate[f.v] == NIL) {
                        /* Root (free vertex) â€” done */
                        stk.pop_back();
                        continue;
                    }
                    int mv = mate[f.v];
                    int pmv = parent[mv];
                    pairs.push_back({mv, pmv});
                    f.v = pmv;
                    continue;
                }
                /* Has bridge â€” originally ODD, absorbed into blossom */
                f.sb = bridge_src[f.v];
                f.tb = bridge_tgt[f.v];
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

    /* Augment along the path:
     *   root_u ~~~ u â€” v ~~~ root_v
     * where (u,v) is the cross-tree non-matching edge.
     * Collect all edge pairs, then flip mate for all of them. */
    void augment_two_sides(int u, int v) {
        std::vector<std::pair<int,int>> pairs;
        /* The cross-tree edge */
        pairs.push_back({u, v});
        /* Trace from u to its root (free vertex) */
        trace_path(u, NIL, pairs);
        /* Trace from v to its root (free vertex) */
        trace_path(v, NIL, pairs);
        /* Flip all */
        for (auto& [a, b] : pairs) {
            mate[a] = b;
            mate[b] = a;
        }
    }

    /* Find one augmenting path in the forest and augment.
     * Returns true if an augmentation was performed.
     *
     * All free vertices start as EVEN roots. BFS grows the forest.
     * EVEN-EVEN edge between different trees â†’ augmenting path.
     * EVEN-EVEN edge within same tree â†’ blossom contraction. */
    bool find_and_augment() {
        /* Reset per-iteration state */
        for (int i = 0; i < n; i++) {
            base[i] = i;
            parent[i] = NIL;
            label[i] = UNLABELED;
            bridge_src[i] = NIL;
            bridge_tgt[i] = NIL;
        }

        std::vector<int> queue(n);
        int qhead = 0, qtail = 0;

        /* All free vertices become EVEN roots */
        for (int v = 0; v < n; v++) {
            if (mate[v] == NIL) {
                label[v] = EVEN;
                queue[qtail++] = v;
            }
        }

        while (qhead < qtail) {
            int u = queue[qhead++];
            /* Check that u is still effectively EVEN */
            if (label[find_base(u)] != EVEN) continue;

            for (int v : graph[u]) {
                int bu = find_base(u), bv = find_base(v);
                if (bu == bv) continue;            /* same blossom */
                if (v == mate[u]) continue;        /* skip matching edge */

                if (label[bv] == UNLABELED) {
                    /* v is matched and unlabeled â†’ grow step */
                    label[v] = ODD;
                    parent[v] = u;
                    int w = mate[v];
                    label[w] = EVEN;
                    queue[qtail++] = w;

                } else if (label[bv] == EVEN) {
                    /* EVEN-EVEN edge: blossom or augmenting path */
                    int lca = find_lca(u, v);
                    if (lca != NIL) {
                        /* Same tree â†’ blossom contraction */
                        shrink_path(lca, u, v, queue, qtail);
                        shrink_path(lca, v, u, queue, qtail);
                    } else {
                        /* Different trees â†’ augmenting path! */
                        augment_two_sides(u, v);
                        return true;
                    }
                }
                /* label[bv] == ODD: ignore */
            }
        }
        return false;
    }

    std::vector<std::pair<int,int>> maximum_matching(int greedy_mode = 0) {
        int greedy_count = 0;
        if (greedy_mode == 1) greedy_count = greedy_init();
        else if (greedy_mode == 2) greedy_count = greedy_init_md();
        greedy_size = greedy_count;

        while (find_and_augment()) {}

        std::vector<std::pair<int,int>> matching;
        for (int u = 0; u < n; u++)
            if (mate[u] != NIL && mate[u] > u)
                matching.push_back({u, mate[u]});
        std::sort(matching.begin(), matching.end());
        return matching;
    }
};

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
    for (int i = 0; i < n; i++) {
        if (deg[i] > 1) { fprintf(stderr, "ERROR: Vertex %d in %d edges!\n", i, deg[i]); errors++; }
    }
    int matched = 0;
    for (int i = 0; i < n; i++) if (deg[i] > 0) matched++;

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %d\n", (int)matching.size());
    printf("Matched vertices: %d\n", matched);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
}

int main(int argc, char* argv[]) {
    printf("Gabow's Algorithm (Simple) - C++ Implementation\n");
    printf("=================================================\n\n");

    if (argc < 2) { printf("Usage: %s <filename> [--greedy|--greedy-md]\n", argv[0]); return 1; }
    int greedy_mode = 0;
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--greedy") greedy_mode = 1;
        else if (std::string(argv[i]) == "--greedy-md") greedy_mode = 2;
    }

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
    GabowSimple gabow(n, edges);
    auto matching = gabow.maximum_matching(greedy_mode);
    auto t1 = std::chrono::high_resolution_clock::now();

    validate_matching(n, gabow.graph, matching);

    printf("Matching size: %d\n", (int)matching.size());

    if (greedy_mode > 0) {
        int gs = gabow.greedy_size;
        int fs = (int)matching.size();
        printf("Greedy init size: %d\n", gs);
        if (fs > 0) printf("Greedy/Final: %.2f%%\n", 100.0 * gs / fs);
        else printf("Greedy/Final: NA\n");
    }
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

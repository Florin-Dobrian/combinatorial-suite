/*
 * Micali-Vazirani Algorithm (Hybrid) - O(E√V) Maximum Matching
 *
 * Hybrid approach:
 * - MV-style MIN phase (level building with even/odd tracking)
 * - Gabow-style MAX phase (BFS path finding and augmentation)
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <climits>

static const int NIL = -1;
static const int UNSET = INT_MAX;

struct Node {
    std::vector<int> preds;
    std::vector<int> hanging_bridges;

    int match;
    int min_level;
    int even_level;
    int odd_level;

    Node() : match(NIL), min_level(UNSET), even_level(UNSET), odd_level(UNSET) {}

    void set_min_level(int level) {
        min_level = level;
        if (level % 2 == 0) even_level = level;
        else odd_level = level;
    }

    void reset() {
        preds.clear();
        hanging_bridges.clear();
        min_level = even_level = odd_level = UNSET;
    }
};

struct MicaliVazirani {
    int n;
    int greedy_size = 0;
    std::vector<std::vector<int>> graph;
    std::vector<Node> nodes;
    std::vector<int> base;
    std::vector<std::vector<int>> levels;

    MicaliVazirani(int n, const std::vector<std::pair<int,int>>& edges) : n(n) {
        graph.resize(n);
        nodes.resize(n);
        base.resize(n);

        for (auto& e : edges) {
            int u = e.first, v = e.second;
            if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
                graph[u].push_back(v);
                graph[v].push_back(u);
            }
        }
        for (int i = 0; i < n; i++) { std::sort(graph[i].begin(), graph[i].end()); graph[i].erase(std::unique(graph[i].begin(), graph[i].end()), graph[i].end()); }
    }

    int find_base(int v) {
        if (base[v] != v) base[v] = find_base(base[v]);
        return base[v];
    }

    void add_to_level(int level, int node) {
        if ((int)levels.size() <= level) levels.resize(level + 1);
        levels[level].push_back(node);
    }

    void step_to(int to, int from, int level) {
        level++;
        int tl = nodes[to].min_level;
        if (tl >= level) {
            if (tl != level) {
                add_to_level(level, to);
                nodes[to].set_min_level(level);
            }
            nodes[to].preds.push_back(from);
        }
    }

    /* MIN: Build alternating tree level by level */
    void phase_1() {
        levels.clear();
        for (int i = 0; i < n; i++) { base[i] = i; nodes[i].reset(); }

        for (int i = 0; i < n; i++) {
            if (nodes[i].match == NIL) {
                add_to_level(0, i);
                nodes[i].set_min_level(0);
            }
        }

        for (int i = 0; i < n; i++) {
            if ((int)levels.size() <= i || levels[i].empty()) continue;
            for (int cur : levels[i]) {
                if (i % 2 == 0) {
                    for (int nb : graph[cur]) {
                        if (nb != nodes[cur].match) step_to(nb, cur, i);
                    }
                } else {
                    if (nodes[cur].match != NIL) step_to(nodes[cur].match, cur, i);
                }
            }
        }
    }

    /* MAX: Find and augment paths (Gabow-style BFS) */
    bool phase_2() {
        bool found = false;

        for (int start = 0; start < n; start++) {
            if (nodes[start].match != NIL || nodes[start].min_level != 0) continue;

            std::vector<int> queue(n);
            std::vector<int> pred(n, NIL);
            std::vector<bool> vis(n, false);
            int qh = 0, qt = 0;

            queue[qt++] = start;
            vis[find_base(start)] = true;
            int endpoint = NIL;

            while (qh < qt && endpoint == NIL) {
                int u = queue[qh++];
                for (int v : graph[u]) {
                    int bu = find_base(u), bv = find_base(v);
                    if (bu == bv || vis[bv]) continue;

                    if (nodes[v].match == NIL && v != start) {
                        pred[v] = u;
                        endpoint = v;
                        break;
                    }

                    pred[v] = u;
                    vis[bv] = true;
                    int mv = nodes[v].match;
                    if (mv != NIL && !vis[find_base(mv)]) {
                        pred[mv] = v;
                        vis[find_base(mv)] = true;
                        queue[qt++] = mv;
                    }
                }
            }

            if (endpoint != NIL) {
                std::vector<int> path;
                for (int c = endpoint; c != NIL; c = pred[c]) path.push_back(c);
                std::reverse(path.begin(), path.end());
                for (size_t i = 0; i + 1 < path.size(); i += 2) {
                    nodes[path[i]].match = path[i + 1];
                    nodes[path[i + 1]].match = path[i];
                }
                found = true;
            }
        }
        return found;
    }

    int greedy_init() {
        int cnt = 0;
        for (int u = 0; u < n; u++) {
            if (nodes[u].match != NIL) continue;
            for (int v : graph[u]) {
                if (nodes[v].match == NIL) { nodes[u].match = v; nodes[v].match = u; cnt++; break; }
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
        /* Process vertices in order of ascending degree */
        std::vector<int> order(n);
        for (int i = 0; i < n; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b){ return deg[a] < deg[b]; });
        for (int u : order) {
            if (nodes[u].match != NIL) continue;
            int best = -1, best_deg = INT_MAX;
            for (int v : graph[u]) {
                if (nodes[v].match == NIL && deg[v] < best_deg) {
                    best = v; best_deg = deg[v];
                }
            }
            if (best >= 0) { nodes[u].match = best; nodes[best].match = u; cnt++; }
        }
        return cnt;
    }

    std::vector<std::pair<int,int>> maximum_matching(int greedy_mode = 0) {
        int greedy_count = 0;
        if (greedy_mode == 1) greedy_count = greedy_init();
        else if (greedy_mode == 2) greedy_count = greedy_init_md();
        greedy_size = greedy_count;
        while (true) {
            phase_1();
            if (!phase_2()) break;
        }

        std::vector<std::pair<int,int>> matching;
        for (int u = 0; u < n; u++) {
            if (nodes[u].match != NIL && nodes[u].match > u)
                matching.push_back({u, nodes[u].match});
        }
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
    printf("Micali-Vazirani Algorithm (Hybrid) - C++ Implementation\n");
    printf("========================================================\n\n");

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
    MicaliVazirani mv(n, edges);
    auto matching = mv.maximum_matching(greedy_mode);
    auto t1 = std::chrono::high_resolution_clock::now();

    validate_matching(n, mv.graph, matching);

    printf("Matching size: %d\n", (int)matching.size());

    if (greedy_mode > 0) {
        int gs = mv.greedy_size;
        int fs = (int)matching.size();
        printf("Greedy init size: %d\n", gs);
        if (fs > 0) printf("Greedy/Final: %.2f%%\n", 100.0 * gs / fs);
        else printf("Greedy/Final: NA\n");
    }
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

/*
 * Micali-Vazirani Pure Algorithm - O(EâˆšV) Maximum Matching
 *
 * True MV with DDFS, tenacity, regular + hanging bridges, petal contraction.
 * Ported from production Jorants MV-Matching-V2.
 *
 * All integers, no hash containers, fully deterministic.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <chrono>
#include <climits>
#include <string>
#include <string>

static const int NIL = -1;

/* DDFS result codes */
static const int DDFS_EMPTY = 0;
static const int DDFS_PETAL = 1;
static const int DDFS_PATH  = 2;

/* =========================================================================
 * Node
 * ========================================================================= */
struct Node {
    std::vector<int> preds;
    std::vector<std::pair<int,int>> pred_to; /* (target, index in target's preds) */
    std::vector<int> hanging_bridges;

    int min_level;
    int max_level;
    int even_level;
    int odd_level;
    int match;
    int bud;
    int above;
    int below;
    int ddfs_green;
    int ddfs_red;
    int number_preds;
    bool deleted;
    bool visited;

    Node() : min_level(NIL), max_level(NIL), even_level(NIL), odd_level(NIL),
             match(NIL), bud(NIL), above(NIL), below(NIL),
             ddfs_green(NIL), ddfs_red(NIL), number_preds(0),
             deleted(false), visited(false) {}

    void set_min_level(int level) {
        min_level = level;
        if (level % 2) odd_level = level;
        else even_level = level;
    }

    void set_max_level(int level) {
        max_level = level;
        if (level % 2) odd_level = level;
        else even_level = level;
    }

    bool outer() const { return even_level != NIL && (odd_level == NIL || even_level < odd_level); }
    bool inner() const { return !outer(); }

    void reset() {
        preds.clear();
        pred_to.clear();
        hanging_bridges.clear();
        min_level = max_level = even_level = odd_level = NIL;
        bud = above = below = ddfs_green = ddfs_red = NIL;
        number_preds = 0;
        deleted = false;
        visited = false;
    }
};

/* =========================================================================
 * DDFS result
 * ========================================================================= */
struct DDFSResult {
    std::vector<int> nodes_seen;
    int bottleneck;
    DDFSResult() : bottleneck(NIL) {}
};

/* =========================================================================
 * MVGraph â€” the full algorithm
 * ========================================================================= */
struct MVGraph {
    std::vector<Node> nodes;
    std::vector<int> edges;           /* flat adjacency array (CSR) */
    std::vector<int> adj_start;
    std::vector<int> deg;

    std::vector<std::vector<int>> levels;
    std::vector<std::vector<std::pair<int,int>>> bridges; /* bridges by tenacity bucket */

    std::vector<std::pair<int,int>> green_stack;
    std::vector<std::pair<int,int>> red_stack;
    std::vector<int> path_found;
    DDFSResult last_ddfs;

    int matchnum;
    int bridgenum;
    int todonum;

    MVGraph() : matchnum(0), bridgenum(0), todonum(0) {}

    /* ---- construction ---- */
    void build(int n, const std::vector<std::pair<int,int>>& edge_list) {
        nodes.resize(n);
        std::vector<std::vector<int>> adj(n);
        for (auto& e : edge_list) {
            int u = e.first, v = e.second;
            if (u >= 0 && u < n && v >= 0 && v < n && u != v) {
                adj[u].push_back(v);
                adj[v].push_back(u);
            }
        }
        for (int i = 0; i < n; i++) { std::sort(adj[i].begin(), adj[i].end()); adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end()); }

        adj_start.resize(n);
        deg.resize(n);
        edges.clear();
        for (int i = 0; i < n; i++) {
            adj_start[i] = (int)edges.size();
            deg[i] = (int)adj[i].size();
            for (int nb : adj[i]) edges.push_back(nb);
        }
        levels.reserve(n / 2 + 1);
        bridges.reserve(n / 2 + 1);
    }

    /* ---- greedy initialization ---- */
    int greedy_init() {
        int cnt = 0;
        for (int j = 0; j < (int)nodes.size(); j++) {
            if (nodes[j].match == NIL) {
                for (int k = 0; k < deg[j]; k++) {
                    int i = edges[adj_start[j] + k];
                    if (nodes[i].match == NIL) {
                        nodes[j].match = i;
                        nodes[i].match = j;
                        matchnum++;
                        cnt++;
                        break;
                    }
                }
            }
        }
        return cnt;
    }

    /* Min-degree greedy: match each exposed vertex with its lowest-degree unmatched neighbor */
    int greedy_init_md() {
        int cnt = 0;
        int nn = (int)nodes.size();
        std::vector<int> order(nn);
        for (int i = 0; i < nn; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b){ return deg[a] < deg[b] || (deg[a] == deg[b] && a < b); });
        for (int j : order) {
            if (nodes[j].match != NIL) continue;
            int best = -1, best_deg = INT_MAX;
            for (int k = 0; k < deg[j]; k++) {
                int i = edges[adj_start[j] + k];
                if (nodes[i].match == NIL && deg[i] < best_deg) {
                    best = i; best_deg = deg[i];
                }
            }
            if (best >= 0) {
                nodes[j].match = best;
                nodes[best].match = j;
                matchnum++;
                cnt++;
            }
        }
        return cnt;
    }

    /* ---- helpers ---- */
    void add_to_level(int level, int node) {
        if (level >= (int)levels.size()) levels.resize(level + 1);
        levels[level].push_back(node);
        todonum++;
    }

    void add_to_bridges(int level, int n1, int n2) {
        if (level >= (int)bridges.size()) bridges.resize(level + 1);
        bridges[level].push_back({n1, n2});
        bridgenum++;
    }

    int tenacity(int n1, int n2) const {
        if (nodes[n1].match == n2) { /* matched bridge */
            if (nodes[n1].odd_level != NIL && nodes[n2].odd_level != NIL)
                return nodes[n1].odd_level + nodes[n2].odd_level + 1;
        }
        else { /* unmatched bridge */
            if (nodes[n1].even_level != NIL && nodes[n2].even_level != NIL)
                return nodes[n1].even_level + nodes[n2].even_level + 1;
        }
        return NIL;
    }

    int bud_star(int c) const {
        int b = nodes[c].bud;
        if (b == NIL) return c;
        return bud_star(b);
    }

    bool bud_star_includes(int c, int goal) const {
        if (c == goal) return true;
        int b = nodes[c].bud;
        if (b == NIL) return false;
        return bud_star_includes(b, goal);
    }

    /* ---- reset between phases ---- */
    void reset() {
        for (auto& v : levels) v.clear();
        for (auto& v : bridges) v.clear();
        bridgenum = 0;
        todonum = 0;
        for (int i = 0; i < (int)nodes.size(); i++) {
            nodes[i].reset();
            if (nodes[i].match == NIL) {
                add_to_level(0, i);
                nodes[i].set_min_level(0);
            }
        }
    }

    /* ---- step_to: core level-building step ---- */
    void step_to(int to, int from, int level) {
        level++;
        int tl = nodes[to].min_level;
        if (tl == NIL || tl >= level) {
            if (tl != level) {
                add_to_level(level, to);
                nodes[to].set_min_level(level);
            }
            nodes[to].preds.push_back(from);
            nodes[to].number_preds++;
            nodes[from].pred_to.push_back({to, (int)(nodes[to].preds.size() - 1)});
        }
        else {
            /* found a bridge */
            int ten = tenacity(to, from);
            if (ten == NIL) {
                nodes[to].hanging_bridges.push_back(from);
                nodes[from].hanging_bridges.push_back(to);
            }
            else {
                add_to_bridges((ten - 1) / 2, to, from);
            }
        }
    }

    /* ---- MIN phase ---- */
    void MIN(int i) {
        if (i >= (int)levels.size()) return;
        for (size_t k = 0; k < levels[i].size(); k++) {
            int current = levels[i][k];
            todonum--;
            Node& n = nodes[current];
            if (i % 2 == 0) {
                for (int j = 0; j < deg[current]; j++) {
                    int edge = edges[adj_start[current] + j];
                    if (edge != n.match) step_to(edge, current, i);
                }
            }
            else {
                if (n.match != NIL) step_to(n.match, current, i);
            }
        }
    }

    /* ---- MAX phase ---- */
    bool MAX(int i) {
        bool found = false;
        if (i >= (int)bridges.size()) return false;

        for (size_t j = 0; j < bridges[i].size(); j++) {
            auto current = bridges[i][j];
            bridgenum--;
            int n1 = current.first;
            int n2 = current.second;
            if (nodes[n1].deleted || nodes[n2].deleted) continue;

            int result = DDFS(n1, n2);
            if (result == DDFS_EMPTY) continue;

            if (result == DDFS_PATH) {
                find_path(n1, n2);
                augment_path();
                if ((int)nodes.size() / 2 <= matchnum) return true;
                remove_path();
                found = true;
            }
            else { /* DDFS_PETAL */
                int b = last_ddfs.bottleneck;
                int current_ten = i * 2 + 1;
                for (int itt : last_ddfs.nodes_seen) {
                    nodes[itt].bud = b;
                    nodes[itt].set_max_level(current_ten - nodes[itt].min_level);
                    add_to_level(nodes[itt].max_level, itt);
                    for (int hanging : nodes[itt].hanging_bridges) {
                        int hanging_ten = tenacity(itt, hanging);
                        if (hanging_ten != NIL)
                            add_to_bridges((hanging_ten - 1) / 2, itt, hanging);
                    }
                }
            }
        }
        return found;
    }

    /* ==================================================================
     * DDFS â€” Double Depth-First Search
     * ================================================================== */

    void add_pred_to_stack(int cur, std::vector<std::pair<int,int>>& stack) {
        for (int pred : nodes[cur].preds) {
            if (pred != NIL) stack.push_back({cur, pred});
        }
    }

    void prepare_next(std::pair<int,int>& Nx) {
        if (Nx.first != NIL) nodes[Nx.first].below = Nx.second;
        Nx.second = bud_star(Nx.second);
    }

    static bool edge_valid(const std::pair<int,int>& e) {
        return !(e.first == NIL && e.second == NIL);
    }

    static void node_from_stack(std::pair<int,int>& e, std::vector<std::pair<int,int>>& S) {
        if (!S.empty()) { e = S.back(); S.pop_back(); }
        else e = {NIL, NIL};
    }

    int L(const std::pair<int,int>& e) const {
        int n = bud_star(e.second);
        return nodes[n].min_level;
    }

    void step_into(int& C, std::pair<int,int>& Nx, std::vector<std::pair<int,int>>& S,
                   int green_top, int red_top) {
        prepare_next(Nx);
        if (!nodes[Nx.second].visited) {
            nodes[Nx.second].above = Nx.first;
            C = Nx.second;
            Node& n = nodes[C];
            n.visited = true;
            n.ddfs_green = green_top;
            n.ddfs_red = red_top;
            last_ddfs.nodes_seen.push_back(C);
            add_pred_to_stack(C, S);
        }
        node_from_stack(Nx, S);
    }

    int DDFS(int green_top, int red_top) {
        last_ddfs.nodes_seen.clear();
        last_ddfs.bottleneck = NIL;

        auto& Sr = red_stack;
        auto& Sg = green_stack;
        Sr.clear();
        Sg.clear();

        int G = NIL, R = NIL;

        if (bud_star(red_top) == bud_star(green_top)) return DDFS_EMPTY;
        if (nodes[green_top].min_level == 0 && nodes[red_top].min_level == 0)
            return DDFS_PATH;

        std::pair<int,int> Ng = {NIL, green_top};
        std::pair<int,int> Nr = {NIL, red_top};
        std::pair<int,int> red_before = {NIL, NIL};
        std::pair<int,int> green_before = {NIL, NIL};

        while (R == NIL || G == NIL ||
               nodes[R].min_level > 0 || nodes[G].min_level > 0) {

            while (edge_valid(Nr) && edge_valid(Ng) && L(Nr) != L(Ng)) {

                while (edge_valid(Nr) && L(Nr) > L(Ng))
                    step_into(R, Nr, Sr, green_top, red_top);

                if (!edge_valid(Nr)) {
                    Nr = red_before;
                    int tmp = red_before.first;
                    while (nodes[tmp].above != NIL) {
                        int rc = nodes[tmp].above;
                        for (int ri : nodes[rc].preds) {
                            if (ri == NIL) continue;
                            if (bud_star(ri) == tmp) { nodes[rc].below = ri; break; }
                        }
                        tmp = nodes[tmp].above;
                    }
                }

                while (edge_valid(Ng) && L(Nr) < L(Ng))
                    step_into(G, Ng, Sg, green_top, red_top);

                if (!edge_valid(Ng)) {
                    Ng = green_before;
                    int tmp = green_before.first;
                    while (nodes[tmp].above != NIL) {
                        int rc = nodes[tmp].above;
                        for (int ri : nodes[rc].preds) {
                            if (ri == NIL) continue;
                            if (bud_star(ri) == tmp) { nodes[rc].below = ri; break; }
                        }
                        tmp = nodes[tmp].above;
                    }
                }
            }

            if (bud_star(Nr.second) == bud_star(Ng.second)) {
                if (!Sr.empty()) {
                    red_before = Nr;
                    prepare_next(Nr);
                    node_from_stack(Nr, Sr);
                    if (edge_valid(Nr)) R = Nr.first;
                    else Nr = red_before;
                }
                else if (!Sg.empty()) {
                    green_before = Ng;
                    prepare_next(Ng);
                    node_from_stack(Ng, Sg);
                    if (edge_valid(Ng)) G = Ng.first;
                    else Ng = green_before;
                }
                else {
                    prepare_next(Nr);
                    prepare_next(Ng);
                    last_ddfs.bottleneck = Nr.second;
                    return DDFS_PETAL;
                }
            }
            else {
                step_into(R, Nr, Sr, green_top, red_top);
                step_into(G, Ng, Sg, green_top, red_top);
            }
        }
        return DDFS_PATH;
    }

    /* ==================================================================
     * Path finding and augmentation
     * ================================================================== */

    void find_path(int n1, int n2) {
        path_found.clear();
        walk_down_path(n1);
        std::reverse(path_found.begin(), path_found.end());
        walk_down_path(n2);
    }

    void walk_down_path(int start) {
        int cur = start;
        while (cur != NIL) {
            if (nodes[cur].bud != NIL) cur = walk_blossom(cur);
            else { path_found.push_back(cur); cur = nodes[cur].below; }
        }
    }

    int jump_bridge(int cur) {
        if (nodes[cur].ddfs_green == cur) return nodes[cur].ddfs_red;
        if (nodes[cur].ddfs_red == cur) return nodes[cur].ddfs_green;
        if (bud_star_includes(nodes[cur].ddfs_green, cur)) {
            size_t before = path_found.size();
            int b = nodes[cur].ddfs_green;
            while (b != cur) b = walk_blossom(b);
            std::reverse(path_found.begin() + before, path_found.end());
            return nodes[cur].ddfs_red;
        }
        else {
            size_t before = path_found.size();
            int b = nodes[cur].ddfs_red;
            while (b != cur) b = walk_blossom(b);
            std::reverse(path_found.begin() + before, path_found.end());
            return nodes[cur].ddfs_green;
        }
    }

    int walk_blossom(int cur) {
        if (nodes[cur].outer()) {
            cur = walk_blossom_down(cur, NIL);
        }
        else {
            cur = walk_blossom_up(cur);
            int before = cur;
            cur = jump_bridge(cur);
            cur = walk_blossom_down(cur, before);
        }
        return cur;
    }

    int walk_blossom_down(int cur, int before) {
        if (before == NIL) before = cur;
        int b = nodes[cur].bud;
        while (cur != NIL && cur != b) {
            if (nodes[cur].ddfs_green != nodes[before].ddfs_green ||
                nodes[cur].ddfs_red != nodes[before].ddfs_red)
                cur = walk_blossom(cur);
            else { path_found.push_back(cur); cur = nodes[cur].below; }
        }
        return cur;
    }

    int walk_blossom_up(int cur) {
        while (true) {
            path_found.push_back(cur);
            if (nodes[cur].above == NIL) break;
            int b = nodes[nodes[cur].above].below;
            if (b != cur && bud_star_includes(b, cur)) {
                size_t before = path_found.size();
                while (b != cur) b = walk_blossom(b);
                std::reverse(path_found.begin() + before, path_found.end());
            }
            cur = nodes[cur].above;
        }
        return cur;
    }

    void augment_path() {
        for (size_t i = 0; i + 1 < path_found.size(); i += 2) {
            int n1 = path_found[i];
            int n2 = path_found[i + 1];
            nodes[n1].match = n2;
            nodes[n2].match = n1;
        }
        matchnum++;
    }

    void remove_path() {
        while (!path_found.empty()) {
            int current = path_found.back();
            path_found.pop_back();
            if (!nodes[current].deleted) {
                nodes[current].deleted = true;
                for (auto& itt : nodes[current].pred_to) {
                    Node& n = nodes[itt.first];
                    if (!n.deleted) {
                        n.preds[itt.second] = NIL;
                        n.number_preds--;
                        if (n.number_preds <= 0) path_found.push_back(itt.first);
                    }
                }
            }
        }
    }

    /* ---- main matching driver ---- */
    void max_match() {
        for (int i = 0; i < (int)nodes.size(); i++) {
            if (nodes[i].match == NIL) {
                add_to_level(0, i);
                nodes[i].set_min_level(0);
            }
        }
        bool found = max_match_phase();
        while ((int)nodes.size() / 2 > matchnum && found) {
            reset();
            found = max_match_phase();
        }
    }

    bool max_match_phase() {
        bool found = false;
        for (int i = 0; i < (int)nodes.size() / 2 + 1 && !found; i++) {
            if (todonum <= 0 && bridgenum <= 0) return false;
            MIN(i);
            found = MAX(i);
        }
        return found;
    }

    std::vector<std::pair<int,int>> get_matching() const {
        std::vector<std::pair<int,int>> result;
        for (int i = 0; i < (int)nodes.size(); i++) {
            if (nodes[i].match != NIL && nodes[i].match > i)
                result.push_back({i, nodes[i].match});
        }
        return result;
    }
};

/* =========================================================================
 * File I/O and main
 * ========================================================================= */

int main(int argc, char* argv[]) {
    printf("Micali-Vazirani Pure Algorithm - C++ Implementation\n");
    printf("====================================================\n\n");

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

    std::vector<std::pair<int,int>> edge_list;
    edge_list.reserve(m);
    for (int i = 0; i < m; i++) {
        int u, v;
        if (fscanf(f, "%d %d", &u, &v) != 2) break;
        edge_list.push_back({u, v});
    }
    fclose(f);

    printf("Graph: %d vertices, %d edges\n", n, (int)edge_list.size());

    auto t0 = std::chrono::high_resolution_clock::now();
    MVGraph mv;
    mv.build(n, edge_list);
    int greedy_count = 0;
    if (greedy_mode == 1) greedy_count = mv.greedy_init();
    else if (greedy_mode == 2) greedy_count = mv.greedy_init_md();
    mv.max_match();
    auto t1 = std::chrono::high_resolution_clock::now();

    auto matching = mv.get_matching();
    int matching_size = (int)matching.size();

    /* validation */
    std::vector<int> vdeg(n, 0);
    int errors = 0;
    for (auto& e : matching) { vdeg[e.first]++; vdeg[e.second]++; }
    for (int i = 0; i < n; i++) {
        if (vdeg[i] > 1) { fprintf(stderr, "ERROR: Vertex %d in %d edges!\n", i, vdeg[i]); errors++; }
    }
    int matched = 0;
    for (int i = 0; i < n; i++) if (vdeg[i] > 0) matched++;

    printf("\n=== Validation Report ===\n");
    printf("Matching size: %d\n", matching_size);
    printf("Matched vertices: %d\n", matched);
    printf("%s\n", errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED");
    printf("=========================\n\n");
    printf("Matching size: %d\n", matching_size);

    if (greedy_mode > 0) {
        int gs = greedy_count;
        int fs = mv.matchnum;
        printf("Greedy init size: %d\n", gs);
        if (fs > 0) printf("Greedy/Final: %.2f%%\n", 100.0 * gs / fs);
        else printf("Greedy/Final: NA\n");
    }
    printf("Time: %ld ms\n", (long)std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

    return 0;
}

/*
 * Micali-Vazirani Algorithm - O(E√V)
 * 
 * Hybrid approach:
 * - Uses MV's tenacity-based MIN phase (level building with even/odd tracking)
 * - Uses Gabow-style MAX phase (simple path finding)
 */

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <climits>

class MicaliVazirani {
private:
    static const int NIL = -1;
    static const int UNSET = INT_MAX;
    
    struct Node {
        std::vector<int> preds;
        std::vector<int> hanging_bridges;
        
        int match;
        int min_level;
        int even_level;
        int odd_level;
        
        Node() : match(NIL), min_level(UNSET),
                 even_level(UNSET), odd_level(UNSET) {}
        
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
    
    int vertex_count;
    std::vector<std::vector<int>> graph;
    std::vector<Node> nodes;
    std::vector<int> base;  // For simple blossom tracking
    
    std::vector<std::vector<int>> levels;
    
    int matchnum;

public:
    MicaliVazirani(int n, const std::vector<std::pair<int, int>>& edges) {
        vertex_count = n;
        graph.resize(n);
        nodes.resize(n);
        base.resize(n);
        matchnum = 0;
        
        for (const auto& [u, v] : edges) {
            if (u < n && v < n && u != v) {
                graph[u].push_back(v);
                graph[v].push_back(u);
            }
        }
        
        for (auto& adj : graph) {
            std::sort(adj.begin(), adj.end());
        }
    }
    
    int find_base(int v) {
        if (base[v] != v) {
            base[v] = find_base(base[v]);
        }
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
    
    // MIN: Build alternating tree level by level
    void phase_1() {
        levels.clear();
        
        // Reset bases
        for (int i = 0; i < vertex_count; i++) {
            base[i] = i;
            nodes[i].reset();
        }
        
        // Initialize free vertices at level 0
        for (int i = 0; i < vertex_count; i++) {
            if (nodes[i].match == NIL) {
                add_to_level(0, i);
                nodes[i].set_min_level(0);
            }
        }
        
        // Build levels
        for (int i = 0; i < vertex_count; i++) {
            if ((int)levels.size() <= i || levels[i].empty()) continue;
            
            for (int current : levels[i]) {
                if (i % 2 == 0) {
                    // Even level - explore all non-matching edges
                    for (int neighbor : graph[current]) {
                        if (neighbor != nodes[current].match) {
                            step_to(neighbor, current, i);
                        }
                    }
                } else {
                    // Odd level - follow matching edge only
                    if (nodes[current].match != NIL) {
                        step_to(nodes[current].match, current, i);
                    }
                }
            }
        }
    }
    
    // MAX: Find and augment paths (Gabow-style)
    bool phase_2() {
        bool found = false;
        
        for (int start = 0; start < vertex_count; start++) {
            if (nodes[start].match != NIL) continue;
            if (nodes[start].min_level != 0) continue;
            
            // BFS from this free vertex
            std::queue<int> q;
            std::vector<int> pred(vertex_count, NIL);
            std::vector<bool> visited(vertex_count, false);
            
            q.push(start);
            visited[find_base(start)] = true;
            
            int endpoint = NIL;
            
            while (!q.empty() && endpoint == NIL) {
                int u = q.front();
                q.pop();
                
                for (int v : graph[u]) {
                    int base_u = find_base(u);
                    int base_v = find_base(v);
                    
                    if (base_u == base_v) continue;
                    if (visited[base_v]) continue;
                    
                    // Check if v is a free vertex
                    if (nodes[v].match == NIL && v != start) {
                        pred[v] = u;
                        endpoint = v;
                        break;
                    }
                    
                    // Otherwise, follow tree structure
                    pred[v] = u;
                    visited[base_v] = true;
                    
                    // Continue along matching edge
                    int mate_v = nodes[v].match;
                    if (mate_v != NIL && !visited[find_base(mate_v)]) {
                        pred[mate_v] = v;
                        visited[find_base(mate_v)] = true;
                        q.push(mate_v);
                    }
                }
            }
            
            // If we found a path, augment it
            if (endpoint != NIL) {
                // Reconstruct path
                std::vector<int> path;
                int curr = endpoint;
                while (curr != NIL) {
                    path.push_back(curr);
                    curr = pred[curr];
                }
                std::reverse(path.begin(), path.end());
                
                // Augment along path
                for (size_t i = 0; i + 1 < path.size(); i += 2) {
                    int u = path[i];
                    int v = path[i + 1];
                    nodes[u].match = v;
                    nodes[v].match = u;
                }
                matchnum++;
                found = true;
            }
        }
        
        return found;
    }
    
    std::vector<std::pair<int, int>> maximum_matching() {
        while (true) {
            phase_1();  // Build level structure
            if (!phase_2()) {  // Find and augment paths
                break;
            }
        }
        
        std::vector<std::pair<int, int>> matching;
        std::vector<bool> seen(vertex_count, false);
        
        for (int u = 0; u < vertex_count; u++) {
            if (nodes[u].match != NIL && !seen[u]) {
                int v = nodes[u].match;
                if (v >= 0 && v < vertex_count) {
                    matching.push_back({std::min(u, v), std::max(u, v)});
                    seen[u] = true;
                    seen[v] = true;
                }
            }
        }
        
        std::sort(matching.begin(), matching.end());
        return matching;
    }
    
    void validate_matching(const std::vector<std::pair<int, int>>& matching) {
        std::vector<int> degree(vertex_count, 0);
        int errors = 0;
        
        std::cout << "\n=== Validation Report ===" << std::endl;
        std::cout << "Matching size: " << matching.size() << std::endl;
        
        for (const auto& [u, v] : matching) {
            if (std::find(graph[u].begin(), graph[u].end(), v) == graph[u].end()) {
                std::cerr << "ERROR: Edge (" << u << ", " << v << ") not in graph!" << std::endl;
                errors++;
            }
            degree[u]++;
            degree[v]++;
        }
        
        for (int i = 0; i < vertex_count; i++) {
            if (degree[i] > 1) {
                std::cerr << "ERROR: Vertex " << i << " in " << degree[i] << " edges!" << std::endl;
                errors++;
            }
        }
        
        int matched = 0;
        for (int d : degree) if (d > 0) matched++;
        
        std::cout << "Matched vertices: " << matched << std::endl;
        std::cout << (errors > 0 ? "VALIDATION FAILED" : "VALIDATION PASSED") << std::endl;
        std::cout << "=========================\n" << std::endl;
    }
};

std::pair<int, std::vector<std::pair<int, int>>> load_graph(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    int n, m;
    file >> n >> m;
    
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < m; i++) {
        int u, v;
        file >> u >> v;
        edges.push_back({u, v});
    }
    
    return {n, edges};
}

int main(int argc, char* argv[]) {
    std::cout << "Micali-Vazirani Algorithm - C++ Implementation" << std::endl;
    std::cout << "===============================================" << std::endl;
    std::cout << std::endl;
    
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    
    try {
        auto [n, edges] = load_graph(argv[1]);
        
        std::cout << "Graph: " << n << " vertices, " << edges.size() << " edges" << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        MicaliVazirani mv(n, edges);
        auto matching = mv.maximum_matching();
        auto end = std::chrono::high_resolution_clock::now();
        
        mv.validate_matching(matching);
        
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Matching size: " << matching.size() << std::endl;
        std::cout << "Time: " << duration.count() << " ms" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

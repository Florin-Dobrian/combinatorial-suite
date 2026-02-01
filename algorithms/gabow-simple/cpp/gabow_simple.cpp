/*
 * Gabow's Algorithm for Maximum Cardinality Matching (Simple Version)
 * Time complexity: O(V * E)
 * 
 * This is the "more implementable" version of Gabow's matching algorithm.
 * Uses efficient blossom handling with path compression.
 * Simpler than Micali-Vazirani but still correct for all general graphs.
 * 
 * Key features:
 * - Integer vertices (0 to n-1)
 * - Deterministic (sorted adjacency lists)
 * - Efficient blossom detection and contraction
 * - Clear separation between search and augmentation phases
 */

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>
#include <chrono>

class GabowSimple {
private:
    int vertex_count;
    std::vector<std::vector<int>> graph;
    std::vector<int> mate;      // mate[v] = matched vertex, or -1
    std::vector<int> base;      // base[v] = current base of blossom containing v
    std::vector<int> parent;    // parent[v] = predecessor in alternating tree
    std::vector<bool> blossom;  // blossom[v] = true if v is a blossom base
    std::vector<bool> visited;  // visited during current search
    std::queue<int> q;          // BFS queue
    
    const int NIL = -1;

public:
    GabowSimple(int n, const std::vector<std::pair<int, int>>& edges) {
        vertex_count = n;
        graph.resize(n);
        mate.assign(n, NIL);
        base.resize(n);
        parent.resize(n);
        blossom.resize(n);
        visited.resize(n);
        
        // Build adjacency list
        for (const auto& [u, v] : edges) {
            if (u < n && v < n && u != v) {
                graph[u].push_back(v);
                graph[v].push_back(u);
            }
        }
        
        // Sort for determinism
        for (auto& adj : graph) {
            std::sort(adj.begin(), adj.end());
        }
    }
    
    // Find base of blossom containing v (with path compression)
    int find_base(int v) {
        if (base[v] != v) {
            base[v] = find_base(base[v]);
        }
        return base[v];
    }
    
    // Find lowest common ancestor in alternating tree
    int find_lca(int u, int v) {
        std::vector<bool> path(vertex_count, false);
        
        // Mark path from u to root
        int safety = 0;
        while (safety < vertex_count) {
            u = find_base(u);
            path[u] = true;
            if (mate[u] == NIL) break;
            if (parent[mate[u]] == NIL) break;
            u = parent[mate[u]];
            safety++;
        }
        
        // Find first common ancestor from v
        safety = 0;
        while (safety < vertex_count) {
            v = find_base(v);
            if (path[v]) return v;
            if (mate[v] == NIL) break;
            if (parent[mate[v]] == NIL) break;
            v = parent[mate[v]];
            safety++;
        }
        
        return NIL;
    }
    
    // Mark vertices in blossom
    void mark_blossom(int u, int lca) {
        int safety = 0;
        while (find_base(u) != lca && safety < vertex_count) {
            int bv = find_base(u);
            int bw = find_base(mate[u]);
            
            blossom[bv] = blossom[bw] = true;
            
            if (!visited[bw]) {
                visited[bw] = true;
                q.push(bw);
            }
            
            if (mate[u] == NIL || parent[mate[u]] == NIL) break;
            u = parent[mate[u]];
            safety++;
        }
    }
    
    // Contract blossom
    void contract_blossom(int u, int v) {
        int lca = find_lca(u, v);
        
        blossom.assign(vertex_count, false);
        mark_blossom(u, lca);
        mark_blossom(v, lca);
        
        // Update bases
        for (int i = 0; i < vertex_count; i++) {
            if (blossom[find_base(i)]) {
                base[i] = lca;
                if (!visited[i]) {
                    visited[i] = true;
                    q.push(i);
                }
            }
        }
    }
    
    // Find augmenting path from start using BFS
    bool find_augmenting_path(int start) {
        // Initialize for this search
        for (int i = 0; i < vertex_count; i++) {
            base[i] = i;
            parent[i] = NIL;
        }
        visited.assign(vertex_count, false);
        
        // Clear queue and start BFS from unmatched vertex
        while (!q.empty()) q.pop();
        q.push(start);
        visited[start] = true;
        
        int iterations = 0;
        while (!q.empty() && iterations < vertex_count * vertex_count) {
            iterations++;
            int u = q.front();
            q.pop();
            
            for (int v : graph[u]) {
                int base_u = find_base(u);
                int base_v = find_base(v);
                
                if (base_u == base_v) continue;  // Same blossom
                
                if (mate[v] == NIL) {
                    // Found augmenting path!
                    parent[v] = u;
                    return true;
                }
                
                if (!visited[base_v]) {
                    // v is matched, extend alternating tree
                    parent[v] = u;
                    visited[base_v] = true;
                    
                    int w = mate[v];
                    visited[find_base(w)] = true;
                    q.push(w);
                } else {
                    // Both in tree - potential blossom
                    int root_u = base_u;
                    int safety = 0;
                    while (mate[root_u] != NIL && safety < vertex_count) {
                        if (parent[mate[root_u]] == NIL) break;
                        root_u = find_base(parent[mate[root_u]]);
                        safety++;
                    }
                    
                    int root_v = base_v;
                    safety = 0;
                    while (mate[root_v] != NIL && safety < vertex_count) {
                        if (parent[mate[root_v]] == NIL) break;
                        root_v = find_base(parent[mate[root_v]]);
                        safety++;
                    }
                    
                    if (root_u == root_v) {
                        // Same tree - this is a blossom!
                        contract_blossom(u, v);
                    }
                }
            }
        }
        
        if (iterations >= vertex_count * vertex_count) {
            std::cerr << "Warning: BFS loop timeout in find_augmenting_path" << std::endl;
        }
        
        return false;
    }
    
    // Augment along path
    void augment_path(int v) {
        while (v != NIL) {
            int pv = parent[v];
            int ppv = mate[pv];
            mate[v] = pv;
            mate[pv] = v;
            v = ppv;
        }
    }
    
    // Find maximum matching
    std::vector<std::pair<int, int>> maximum_matching() {
        // Try to find augmenting paths from each unmatched vertex
        bool found = true;
        int iterations = 0;
        
        while (found) {
            found = false;
            iterations++;
            
            for (int v = 0; v < vertex_count; v++) {
                if (mate[v] == NIL) {
                    if (find_augmenting_path(v)) {
                        // Find the endpoint of the augmenting path
                        for (int u = 0; u < vertex_count; u++) {
                            if (mate[u] == NIL && parent[u] != NIL) {
                                augment_path(u);
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            
            // Safety check
            if (iterations > vertex_count) {
                std::cerr << "Warning: Too many iterations, stopping" << std::endl;
                break;
            }
        }
        
        // Build result
        std::vector<std::pair<int, int>> matching;
        std::vector<bool> seen(vertex_count, false);
        
        for (int u = 0; u < vertex_count; u++) {
            if (mate[u] != NIL && !seen[u]) {
                int v = mate[u];
                matching.push_back({std::min(u, v), std::max(u, v)});
                seen[u] = true;
                seen[v] = true;
            }
        }
        
        std::sort(matching.begin(), matching.end());
        return matching;
    }
    
    // Validate matching
    void validate_matching(const std::vector<std::pair<int, int>>& matching) {
        std::vector<int> degree(vertex_count, 0);
        int errors = 0;
        
        std::cout << "\n=== Validation Report ===" << std::endl;
        std::cout << "Matching size: " << matching.size() << std::endl;
        
        for (const auto& [u, v] : matching) {
            // Check edge exists
            if (std::find(graph[u].begin(), graph[u].end(), v) == graph[u].end()) {
                std::cerr << "ERROR: Edge (" << u << ", " << v << ") not in graph!" << std::endl;
                errors++;
            }
            degree[u]++;
            degree[v]++;
        }
        
        // Check each vertex in at most one edge
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

// Load graph from file
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
    std::cout << "Gabow's Algorithm (Simple Version) - C++ Implementation" << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << std::endl;
    
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    
    try {
        auto [n, edges] = load_graph(argv[1]);
        
        std::cout << "Graph: " << n << " vertices, " << edges.size() << " edges" << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        GabowSimple gabow(n, edges);
        auto matching = gabow.maximum_matching();
        auto end = std::chrono::high_resolution_clock::now();
        
        gabow.validate_matching(matching);
        
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Matching size: " << matching.size() << std::endl;
        std::cout << "Time: " << duration.count() << " ms" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

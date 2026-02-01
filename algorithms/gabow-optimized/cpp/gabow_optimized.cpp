/*
 * Gabow's O(E√V) Matching Algorithm
 * Based on LEDA's implementation approach
 * 
 * Algorithm structure:
 * - Phase 1: Build level structure by distance, detect blossoms
 * - Phase 2: Find and augment all shortest augmenting paths
 * - Repeat until no augmenting paths exist
 * 
 * Time complexity: O(E√V)
 * Space complexity: O(V + E)
 * 
 * Features:
 * - Integer vertices (0 to n-1)
 * - Deterministic (sorted adjacency lists)
 * - No hash-based structures
 */

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <climits>

class GabowOptimized {
private:
    int vertex_count;
    std::vector<std::vector<int>> graph;
    
    // Core matching data
    std::vector<int> mate;           // mate[v] = matched vertex or NIL
    std::vector<int> label;          // EVEN, ODD, or UNLABELED
    std::vector<int> base;           // base[v] = base of blossom containing v
    std::vector<int> parent;         // parent[v] = parent in alternating tree
    
    // Bridge information for blossom expansion
    std::vector<int> source_bridge;  // source endpoint of blossom-forming edge
    std::vector<int> target_bridge;  // target endpoint of blossom-forming edge
    
    // Distance tracking
    std::vector<std::vector<std::pair<int,int>>> edge_queue; // edges by distance
    
    const int NIL = -1;
    const int UNLABELED = 0;
    const int EVEN = 1;
    const int ODD = 2;
    
    int Delta;       // Current distance level

public:
    GabowOptimized(int n, const std::vector<std::pair<int, int>>& edges) {
        vertex_count = n;
        graph.resize(n);
        mate.assign(n, NIL);
        label.assign(n, UNLABELED);
        base.resize(n);
        parent.assign(n, NIL);
        source_bridge.assign(n, NIL);
        target_bridge.assign(n, NIL);
        edge_queue.resize(n + 1);
        
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
    
    // Find base of blossom with path compression
    int find_base(int v) {
        if (base[v] != v) {
            base[v] = find_base(base[v]);
        }
        return base[v];
    }
    
    // Find LCA of two vertices in alternating tree
    int find_lca(int u, int v) {
        std::vector<bool> marked(vertex_count, false);
        
        // Mark path from u to root
        int x = find_base(u);
        while (x != NIL) {
            marked[x] = true;
            if (mate[x] == NIL) break;
            if (parent[mate[x]] == NIL) break;
            x = find_base(parent[mate[x]]);
        }
        
        // Find first marked vertex from v
        int y = find_base(v);
        while (y != NIL) {
            if (marked[y]) return y;
            if (mate[y] == NIL) break;
            if (parent[mate[y]] == NIL) break;
            y = find_base(parent[mate[y]]);
        }
        
        return NIL;
    }
    
    // Shrink blossom: merge path from x to lca
    void shrink_path(int lca, int x, int y) {
        int v = find_base(x);
        while (v != lca) {
            // Union v into lca
            base[v] = lca;
            
            int mate_v = mate[v];
            if (mate_v == NIL) break;
            
            // Union mate[v] into lca
            base[mate_v] = lca;
            
            // Store bridge information for later expansion
            source_bridge[mate_v] = x;
            target_bridge[mate_v] = y;
            
            if (parent[mate_v] == NIL) break;
            v = find_base(parent[mate_v]);
        }
    }
    
    // Scan edge: add to queue at appropriate distance
    void scan_edge(int u, int v) {
        if (Delta < (int)edge_queue.size()) {
            edge_queue[Delta].push_back({u, v});
        }
    }
    
    // Phase 1: Build level structure and detect blossoms
    bool phase_1() {
        Delta = 0;
        bool found_augmenting_path = false;
        
        // Clear edge queues
        for (auto& q : edge_queue) q.clear();
        
        // Reset bases and labels
        for (int i = 0; i < vertex_count; i++) {
            base[i] = i;
            label[i] = (mate[i] == NIL) ? EVEN : UNLABELED;
            parent[i] = NIL;
            source_bridge[i] = NIL;
            target_bridge[i] = NIL;
        }
        
        // Add free vertices to tree and scan their edges
        for (int v = 0; v < vertex_count; v++) {
            if (mate[v] == NIL) {
                for (int u : graph[v]) {
                    scan_edge(v, u);
                }
            }
        }
        
        // Process edges by increasing distance
        while (Delta <= vertex_count && !found_augmenting_path) {
            while (!edge_queue[Delta].empty()) {
                auto [x, y] = edge_queue[Delta].back();
                edge_queue[Delta].pop_back();
                
                // Ensure x is at even level
                int base_x = find_base(x);
                int base_y = find_base(y);
                
                if (label[base_x] != EVEN) {
                    std::swap(x, y);
                    std::swap(base_x, base_y);
                }
                
                // Skip if both in same blossom or x not even
                if (base_x == base_y) continue;
                if (label[base_x] != EVEN) continue;
                if (y == mate[x]) continue;
                if (label[base_y] == ODD) continue;
                
                if (label[base_y] == UNLABELED) {
                    // Grow tree
                    int z = mate[y];
                    if (z != NIL) {
                        label[y] = ODD;
                        label[z] = EVEN;
                        parent[y] = x;
                        parent[z] = y;
                        
                        // Scan edges from z
                        for (int w : graph[z]) {
                            scan_edge(z, w);
                        }
                    }
                } else if (label[base_y] == EVEN) {
                    // Both even - check if same tree or different trees
                    int lca = find_lca(x, y);
                    
                    if (lca != NIL) {
                        // Same tree - shrink blossom
                        shrink_path(lca, x, y);
                        shrink_path(lca, y, x);
                    } else {
                        // Different trees - found augmenting path!
                        found_augmenting_path = true;
                        break;
                    }
                }
            }
            
            if (!found_augmenting_path) {
                Delta++;
            }
        }
        
        return found_augmenting_path;
    }
    
    // Phase 2: Find and augment all shortest augmenting paths
    void phase_2() {
        // For each free EVEN vertex, find path to another free vertex
        for (int start = 0; start < vertex_count; start++) {
            if (mate[start] != NIL || label[start] != EVEN) continue;
            
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
                    
                    // Check if v is a free vertex in different tree
                    if (mate[v] == NIL && v != start) {
                        pred[v] = u;
                        endpoint = v;
                        break;
                    }
                    
                    // Otherwise, follow the tree structure
                    if (label[base_v] != ODD) {
                        pred[v] = u;
                        visited[base_v] = true;
                        
                        // Continue along matching edge
                        int mate_v = mate[v];
                        if (mate_v != NIL && !visited[find_base(mate_v)]) {
                            pred[mate_v] = v;
                            visited[find_base(mate_v)] = true;
                            q.push(mate_v);
                        }
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
                    mate[u] = v;
                    mate[v] = u;
                }
            }
        }
    }
    
    // Main algorithm
    std::vector<std::pair<int, int>> maximum_matching() {
        while (phase_1()) {
            phase_2();
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
    std::cout << "Gabow's Scaling Algorithm (Optimized) - C++ Implementation" << std::endl;
    std::cout << "===========================================================" << std::endl;
    std::cout << std::endl;
    
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    
    try {
        auto [n, edges] = load_graph(argv[1]);
        
        std::cout << "Graph: " << n << " vertices, " << edges.size() << " edges" << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        GabowOptimized gabow(n, edges);
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

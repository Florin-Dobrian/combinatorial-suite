/*
 * Hopcroft-Karp Algorithm for Maximum Bipartite Matching (C++ Implementation)
 * Time complexity: O(E * sqrt(V))
 * 
 * Uses integers for vertices and deterministic data structures.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <chrono>
#include <climits>

class HopcroftKarp {
private:
    int left_count;
    int right_count;
    std::vector<std::vector<int>> graph;  // graph[left_node] = list of right nodes
    std::vector<int> pair_left;   // pair_left[u] = matched right node (or -1)
    std::vector<int> pair_right;  // pair_right[v] = matched left node (or -1)
    std::vector<int> dist;        // distance in BFS
    
    const int NIL = -1;
    
    bool bfs() {
        std::queue<int> q;
        
        // Initialize distances and queue with unmatched left nodes
        for (int u = 0; u < left_count; u++) {
            if (pair_left[u] == NIL) {
                dist[u] = 0;
                q.push(u);
            } else {
                dist[u] = INT_MAX;
            }
        }
        
        dist[left_count] = INT_MAX;  // NIL node is at index left_count
        
        // BFS
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            if (dist[u] < dist[left_count]) {
                for (int v : graph[u]) {
                    int paired_node = pair_right[v];
                    
                    if (paired_node == NIL) {
                        paired_node = left_count;  // Use NIL index
                    }
                    
                    if (dist[paired_node] == INT_MAX) {
                        dist[paired_node] = dist[u] + 1;
                        if (pair_right[v] != NIL) {
                            q.push(pair_right[v]);
                        }
                    }
                }
            }
        }
        
        return dist[left_count] != INT_MAX;
    }
    
    bool dfs(int u) {
        if (u == NIL) {
            return true;
        }
        
        for (int v : graph[u]) {
            int paired_node = pair_right[v];
            if (paired_node == NIL) {
                paired_node = left_count;  // Use NIL index
            }
            
            if (dist[paired_node] == dist[u] + 1) {
                if (dfs(pair_right[v])) {
                    pair_right[v] = u;
                    pair_left[u] = v;
                    return true;
                }
            }
        }
        
        dist[u] = INT_MAX;
        return false;
    }

public:
    HopcroftKarp(int left_count, int right_count, const std::vector<std::pair<int, int>>& edges) {
        this->left_count = left_count;
        this->right_count = right_count;
        
        graph.resize(left_count);
        
        // Build adjacency list
        for (const auto& [u, v] : edges) {
            if (u >= 0 && u < left_count && v >= 0 && v < right_count) {
                graph[u].push_back(v);
            }
        }
        
        // Sort adjacency lists for deterministic iteration
        for (auto& adj : graph) {
            std::sort(adj.begin(), adj.end());
        }
        
        // Initialize pairs to NIL
        pair_left.assign(left_count, NIL);
        pair_right.assign(right_count, NIL);
        
        // Distance array (size left_count + 1 to include NIL at index left_count)
        dist.resize(left_count + 1);
    }
    
    std::vector<std::pair<int, int>> maximum_matching() {
        // While there exist augmenting paths
        while (bfs()) {
            for (int u = 0; u < left_count; u++) {
                if (pair_left[u] == NIL) {
                    dfs(u);
                }
            }
        }
        
        // Build the matching vector
        std::vector<std::pair<int, int>> matching;
        for (int u = 0; u < left_count; u++) {
            if (pair_left[u] != NIL) {
                matching.push_back({u, pair_left[u]});
            }
        }
        
        std::sort(matching.begin(), matching.end());
        return matching;
    }
    
    void validate_matching(const std::vector<std::pair<int, int>>& matching) {
        std::vector<int> left_degree(left_count, 0);
        std::vector<int> right_degree(right_count, 0);
        int errors = 0;
        
        std::cerr << "\n=== Validation Report ===" << std::endl;
        std::cerr << "Matching size (claimed): " << matching.size() << std::endl;
        
        for (const auto& [u, v] : matching) {
            // Check if edge exists in original graph
            if (std::find(graph[u].begin(), graph[u].end(), v) == graph[u].end()) {
                std::cerr << "ERROR: Edge (" << u << ", " << v 
                         << ") in matching but NOT in original graph!" << std::endl;
                errors++;
            }
            
            left_degree[u]++;
            right_degree[v]++;
        }
        
        for (int i = 0; i < left_count; i++) {
            if (left_degree[i] > 1) {
                std::cerr << "ERROR: Left node " << i << " appears in " 
                         << left_degree[i] << " edges (should be at most 1)!" << std::endl;
                errors++;
            }
        }
        
        for (int i = 0; i < right_count; i++) {
            if (right_degree[i] > 1) {
                std::cerr << "ERROR: Right node " << i << " appears in " 
                         << right_degree[i] << " edges (should be at most 1)!" << std::endl;
                errors++;
            }
        }
        
        int unique_left = 0, unique_right = 0;
        for (int deg : left_degree) if (deg > 0) unique_left++;
        for (int deg : right_degree) if (deg > 0) unique_right++;
        
        std::cerr << "Number of edges in matching: " << matching.size() << std::endl;
        std::cerr << "Left nodes matched: " << unique_left << std::endl;
        std::cerr << "Right nodes matched: " << unique_right << std::endl;
        
        if (errors > 0) {
            std::cerr << "VALIDATION FAILED: " << errors << " errors found" << std::endl;
        } else {
            std::cerr << "VALIDATION PASSED: Matching is valid" << std::endl;
        }
        std::cerr << "=========================\n" << std::endl;
    }
};

bool load_graph_from_file(const std::string& filename,
                          int& left_count,
                          int& right_count,
                          std::vector<std::pair<int, int>>& edges) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    int edge_count;
    file >> left_count >> right_count >> edge_count;
    
    edges.clear();
    for (int i = 0; i < edge_count; i++) {
        int u, v;
        file >> u >> v;
        edges.push_back({u, v});
    }
    
    file.close();
    return true;
}

int main(int argc, char* argv[]) {
    std::cout << "Hopcroft-Karp Algorithm - C++ Implementation" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << std::endl;
    
    if (argc > 1) {
        std::string filename = argv[1];
        std::cout << "Loading graph from: " << filename << std::endl;
        
        int left_count, right_count;
        std::vector<std::pair<int, int>> edges;
        
        if (!load_graph_from_file(filename, left_count, right_count, edges)) {
            std::cerr << "Error: Could not open file '" << filename << "'\n";
            return 1;
        }
        
        std::cout << "File: " << filename << std::endl;
        std::cout << "Graph: " << left_count << " left nodes, " 
                  << right_count << " right nodes, " 
                  << edges.size() << " edges" << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        HopcroftKarp hk(left_count, right_count, edges);
        auto matching = hk.maximum_matching();
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        hk.validate_matching(matching);
        
        std::cout << "Matching size: " << matching.size() << std::endl;
        std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    } else {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    
    return 0;
}

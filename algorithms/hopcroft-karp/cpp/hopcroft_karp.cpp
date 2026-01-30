/*
 * Hopcroft-Karp Algorithm for Maximum Bipartite Matching (C++ Implementation)
 * Time complexity: O(E * sqrt(V))
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <chrono>
#include <algorithm>

class HopcroftKarp {
private:
    std::unordered_set<std::string> left;
    std::unordered_set<std::string> right;
    std::unordered_map<std::string, std::vector<std::string>> graph;
    std::unordered_map<std::string, std::string> pair_left;
    std::unordered_map<std::string, std::string> pair_right;
    std::unordered_map<std::string, int> dist;
    
    const std::string NIL = "";
    
    bool bfs() {
        std::queue<std::string> q;
        
        // Initialize distances and queue with unmatched left nodes
        for (const auto& u : left) {
            if (pair_left[u].empty()) {
                dist[u] = 0;
                q.push(u);
            } else {
                dist[u] = INT32_MAX;
            }
        }
        
        dist[NIL] = INT32_MAX;
        
        // BFS
        while (!q.empty()) {
            std::string u = q.front();
            q.pop();
            
            if (dist[u] < dist[NIL]) {
                for (const auto& v : graph[u]) {
                    std::string paired_node = pair_right[v];
                    
                    if (dist[paired_node] == INT32_MAX) {
                        dist[paired_node] = dist[u] + 1;
                        if (!paired_node.empty()) {
                            q.push(paired_node);
                        }
                    }
                }
            }
        }
        
        return dist[NIL] != INT32_MAX;
    }
    
    bool dfs(const std::string& u) {
        if (u.empty()) {
            return true;
        }
        
        for (const auto& v : graph[u]) {
            std::string paired_node = pair_right[v];
            
            if (dist[paired_node] == dist[u] + 1) {
                if (dfs(paired_node)) {
                    pair_right[v] = u;
                    pair_left[u] = v;
                    return true;
                }
            }
        }
        
        dist[u] = INT32_MAX;
        return false;
    }

public:
    HopcroftKarp(const std::vector<std::string>& left_nodes,
                 const std::vector<std::string>& right_nodes,
                 const std::vector<std::pair<std::string, std::string>>& edges) {
        
        left.insert(left_nodes.begin(), left_nodes.end());
        right.insert(right_nodes.begin(), right_nodes.end());
        
        // Build adjacency list
        for (const auto& edge : edges) {
            if (left.count(edge.first) && right.count(edge.second)) {
                graph[edge.first].push_back(edge.second);
            }
        }
        
        // Initialize pairs
        for (const auto& u : left) {
            pair_left[u] = "";
        }
        for (const auto& v : right) {
            pair_right[v] = "";
        }
    }
    
    std::vector<std::pair<std::string, std::string>> maximum_matching() {
        int matching_size = 0;
        
        // While there exist augmenting paths
        while (bfs()) {
            for (const auto& u : left) {
                if (pair_left[u].empty()) {
                    if (dfs(u)) {
                        matching_size++;
                    }
                }
            }
        }
        
        // Build the matching vector
        std::vector<std::pair<std::string, std::string>> matching;
        for (const auto& u : left) {
            if (!pair_left[u].empty()) {
                matching.push_back({u, pair_left[u]});
            }
        }
        
        return matching;
    }
    
    int get_matching_size() const {
        int size = 0;
        for (const auto& u : left) {
            if (!pair_left.at(u).empty()) {
                size++;
            }
        }
        return size;
    }
};

// Load graph from file
bool load_graph_from_file(const std::string& filename,
                          std::vector<std::string>& left,
                          std::vector<std::string>& right,
                          std::vector<std::pair<std::string, std::string>>& edges) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    int left_count, right_count, edge_count;
    file >> left_count >> right_count >> edge_count;
    
    left.clear();
    right.clear();
    edges.clear();
    
    for (int i = 0; i < left_count; i++) {
        left.push_back("L" + std::to_string(i));
    }
    for (int i = 0; i < right_count; i++) {
        right.push_back("R" + std::to_string(i));
    }
    
    for (int i = 0; i < edge_count; i++) {
        int u, v;
        file >> u >> v;
        edges.push_back({"L" + std::to_string(u), "R" + std::to_string(v)});
    }
    
    file.close();
    return true;
}

// Generate a large random bipartite graph for benchmarking
std::vector<std::pair<std::string, std::string>> generate_large_graph(
    int left_size, int right_size, int edges_per_left_node) {
    
    std::vector<std::pair<std::string, std::string>> edges;
    
    for (int i = 0; i < left_size; i++) {
        for (int j = 0; j < edges_per_left_node; j++) {
            int right_idx = (i * edges_per_left_node + j) % right_size;
            edges.push_back({"L" + std::to_string(i), "R" + std::to_string(right_idx)});
        }
    }
    
    return edges;
}

void run_example(const std::vector<std::string>& left,
                 const std::vector<std::string>& right,
                 const std::vector<std::pair<std::string, std::string>>& edges,
                 const std::string& description) {
    std::cout << description << "\n";
    std::cout << "Graph: " << left.size() << " left nodes, " 
              << right.size() << " right nodes, " 
              << edges.size() << " edges\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    HopcroftKarp hk(left, right, edges);
    auto matching = hk.maximum_matching();
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Matching size: " << matching.size() << "\n";
    if (matching.size() <= 10) {
        std::cout << "Matching: ";
        for (const auto& edge : matching) {
            std::cout << "(" << edge.first << "," << edge.second << ") ";
        }
        std::cout << "\n";
    }
    std::cout << "Execution time: " << duration.count() << " ms\n";
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    std::cout << "C++ Hopcroft-Karp Implementation\n";
    std::cout << "=================================\n\n";
    
    // Check if a file was provided
    if (argc > 1) {
        std::string filename = argv[1];
        std::cout << "Loading graph from: " << filename << "\n";
        
        std::vector<std::string> left, right;
        std::vector<std::pair<std::string, std::string>> edges;
        
        if (!load_graph_from_file(filename, left, right, edges)) {
            std::cerr << "Error: Could not open file '" << filename << "'\n";
            return 1;
        }
        
        run_example(left, right, edges, "File: " + filename);
    } else {
        // Run built-in examples
        std::cout << "Running built-in examples (use: ./hopcroft_karp_cpp <filename> to load from file)\n\n";
        
        // Small example
        std::vector<std::string> left = {"A", "B", "C", "D"};
        std::vector<std::string> right = {"1", "2", "3", "4"};
        std::vector<std::pair<std::string, std::string>> edges = {
            {"A", "1"}, {"A", "2"},
            {"B", "2"}, {"B", "3"},
            {"C", "3"}, {"C", "4"},
            {"D", "4"}
        };
        
        run_example(left, right, edges, "Small example:");
        
        // Benchmark with larger graph
        std::cout << "Benchmarking with larger graph...\n";
        int left_size = 1000;
        int right_size = 1000;
        int edges_per_node = 10;
        
        std::vector<std::string> large_left, large_right;
        for (int i = 0; i < left_size; i++) {
            large_left.push_back("L" + std::to_string(i));
        }
        for (int i = 0; i < right_size; i++) {
            large_right.push_back("R" + std::to_string(i));
        }
        
        auto large_edges = generate_large_graph(left_size, right_size, edges_per_node);
        
        run_example(large_left, large_right, large_edges, "Large benchmark:");
    }
    
    return 0;
}

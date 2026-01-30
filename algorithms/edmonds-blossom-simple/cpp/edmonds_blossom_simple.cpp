/*
 * Edmonds' Blossom Algorithm (Simple Version) for Maximum Cardinality Matching
 * Time complexity: O(V⁴)
 * 
 * This is the straightforward implementation that finds one augmenting path per iteration.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <string>
#include <chrono>
#include <cstdlib>
#include <ctime>

class EdmondsBlossomSimple {
private:
    std::unordered_set<std::string> vertices;
    std::unordered_map<std::string, std::unordered_set<std::string>> graph;
    std::unordered_map<std::string, std::string> mate;
    
    // For finding augmenting paths
    std::unordered_map<std::string, std::string> parent;
    std::unordered_map<std::string, std::string> base;
    std::unordered_map<std::string, std::string> in_blossom;
    
    bool is_outer(const std::string& v) {
        if (parent.find(v) == parent.end()) {
            return false;
        }
        
        std::string current = v;
        int distance = 0;
        std::unordered_set<std::string> visited;
        
        while (current != parent[current] && visited.find(current) == visited.end()) {
            visited.insert(current);
            current = parent[current];
            distance++;
        }
        
        return distance % 2 == 0;
    }
    
    std::string find_blossom(const std::string& v, const std::string& w) {
        // Mark all ancestors of v
        std::unordered_set<std::string> path_v;
        std::string current = v;
        
        while (parent.find(current) != parent.end()) {
            path_v.insert(base[current]);
            if (current == parent[current]) break;
            current = parent[current];
        }
        
        // Find first common ancestor from w
        current = w;
        while (parent.find(current) != parent.end()) {
            if (path_v.find(base[current]) != path_v.end()) {
                return base[current];
            }
            if (current == parent[current]) break;
            current = parent[current];
        }
        
        return v; // Fallback
    }
    
    void contract_blossom(const std::string& blossom_base, const std::string& v, const std::string& w) {
        // Path from v to base
        std::string current = v;
        while (base[current] != blossom_base) {
            base[current] = blossom_base;
            if (!mate[current].empty()) {
                base[mate[current]] = blossom_base;
                current = parent[mate[current]];
            } else {
                break;
            }
        }
        
        // Path from w to base
        current = w;
        while (base[current] != blossom_base) {
            base[current] = blossom_base;
            if (!mate[current].empty()) {
                base[mate[current]] = blossom_base;
                current = parent[mate[current]];
            } else {
                break;
            }
        }
        
        base[v] = blossom_base;
        base[w] = blossom_base;
    }
    
    std::vector<std::string> build_path(const std::string& v, const std::string& w) {
        std::vector<std::string> path = {w, v};
        
        std::string current = v;
        while (parent.find(current) != parent.end() && parent[current] != current) {
            std::string p = parent[current];
            path.push_back(p);
            current = p;
        }
        
        return path;
    }
    
    void augment(const std::vector<std::string>& path) {
        for (size_t i = 0; i + 1 < path.size(); i += 2) {
            mate[path[i]] = path[i + 1];
            mate[path[i + 1]] = path[i];
        }
    }

public:
    EdmondsBlossomSimple(const std::vector<std::string>& vertex_list,
                         const std::vector<std::pair<std::string, std::string>>& edges) {
        
        vertices.insert(vertex_list.begin(), vertex_list.end());
        
        // Build adjacency list (undirected graph)
        for (const auto& edge : edges) {
            if (vertices.count(edge.first) && vertices.count(edge.second) && edge.first != edge.second) {
                graph[edge.first].insert(edge.second);
                graph[edge.second].insert(edge.first);
            }
        }
        
        // Initialize matching
        for (const auto& v : vertices) {
            mate[v] = "";
        }
    }
    
    std::vector<std::string> find_augmenting_path(const std::string& start) {
        // Initialize structures
        parent.clear();
        base.clear();
        in_blossom.clear();
        
        parent[start] = start;
        for (const auto& v : vertices) {
            base[v] = v;
            in_blossom[v] = v;
        }
        
        std::queue<std::string> q;
        q.push(start);
        
        while (!q.empty()) {
            std::string v = q.front();
            q.pop();
            
            for (const auto& w : graph[v]) {
                // Skip if already in tree at same level
                if (base[v] == base[w]) {
                    continue;
                }
                
                // If w is unmatched, we found an augmenting path!
                if (mate[w].empty()) {
                    return build_path(v, w);
                }
                
                // If w is not in tree yet
                if (parent.find(w) == parent.end()) {
                    parent[w] = v;
                    std::string mate_w = mate[w];
                    if (!mate_w.empty()) {
                        parent[mate_w] = w;
                        q.push(mate_w);
                    }
                }
                // If w is in tree at even distance (blossom found)
                else if (is_outer(w)) {
                    std::string blossom_base = find_blossom(v, w);
                    contract_blossom(blossom_base, v, w);
                    q.push(blossom_base);
                }
            }
        }
        
        return std::vector<std::string>();
    }
    
    std::vector<std::pair<std::string, std::string>> maximum_matching() {
        bool improved = true;
        
        while (improved) {
            improved = false;
            
            for (const auto& v : vertices) {
                if (mate[v].empty()) {
                    std::vector<std::string> path = find_augmenting_path(v);
                    if (!path.empty()) {
                        augment(path);
                        improved = true;
                        break;
                    }
                }
            }
        }
        
        // Build the matching set
        std::vector<std::pair<std::string, std::string>> matching;
        std::unordered_set<std::string> seen;
        
        for (const auto& u : vertices) {
            std::string v = mate[u];
            if (!v.empty() && seen.find(v) == seen.end()) {
                matching.push_back({u, v});
                seen.insert(u);
            }
        }
        
        return matching;
    }
};

// Load graph from file
bool load_graph_from_file(const std::string& filename,
                          std::vector<std::string>& vertices,
                          std::vector<std::pair<std::string, std::string>>& edges) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    int vertex_count, edge_count;
    file >> vertex_count >> edge_count;
    
    vertices.clear();
    edges.clear();
    
    for (int i = 0; i < vertex_count; i++) {
        vertices.push_back("V" + std::to_string(i));
    }
    
    for (int i = 0; i < edge_count; i++) {
        int u, v;
        file >> u >> v;
        edges.push_back({"V" + std::to_string(u), "V" + std::to_string(v)});
    }
    
    file.close();
    return true;
}

// Generate random test graph
std::vector<std::pair<std::string, std::string>> generate_test_graph(
    int n, double edge_probability, std::vector<std::string>& vertices) {
    
    vertices.clear();
    for (int i = 0; i < n; i++) {
        vertices.push_back("V" + std::to_string(i));
    }
    
    std::vector<std::pair<std::string, std::string>> edges;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if ((double)rand() / RAND_MAX < edge_probability) {
                edges.push_back({vertices[i], vertices[j]});
            }
        }
    }
    
    return edges;
}

void run_example(const std::vector<std::string>& vertices,
                 const std::vector<std::pair<std::string, std::string>>& edges,
                 const std::string& description) {
    std::cout << description << "\n";
    std::cout << "Graph: " << vertices.size() << " vertices, " << edges.size() << " edges\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    EdmondsBlossomSimple eb(vertices, edges);
    auto matching = eb.maximum_matching();
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
    std::cout << "Edmonds' Blossom Algorithm (Simple) - C++ Implementation\n";
    std::cout << "========================================================\n\n";
    
    srand(time(0));
    
    // Check if a file was provided
    if (argc > 1) {
        std::string filename = argv[1];
        std::cout << "Loading graph from: " << filename << "\n";
        
        std::vector<std::string> vertices;
        std::vector<std::pair<std::string, std::string>> edges;
        
        if (!load_graph_from_file(filename, vertices, edges)) {
            std::cerr << "Error: Could not open file '" << filename << "'\n";
            return 1;
        }
        
        run_example(vertices, edges, "File: " + filename);
    } else {
        // Run built-in examples
        std::cout << "Running built-in examples (use: ./edmonds_blossom_simple_cpp <filename> to load from file)\n\n";
        
        // Example 1: Simple graph with triangle
        std::vector<std::string> vertices = {"A", "B", "C", "D", "E"};
        std::vector<std::pair<std::string, std::string>> edges = {
            {"A", "B"}, {"B", "C"}, {"C", "A"},  // Triangle
            {"C", "D"}, {"D", "E"}
        };
        run_example(vertices, edges, "Example 1: Graph with triangle (blossom)");
        
        // Example 2: Larger random graph
        vertices.clear();
        edges = generate_test_graph(20, 0.2, vertices);
        run_example(vertices, edges, "Example 2: Random graph (20 vertices)");
        
        // Benchmark
        std::cout << "Benchmarking with larger graph...\n";
        vertices.clear();
        edges = generate_test_graph(100, 0.1, vertices);
        run_example(vertices, edges, "Benchmark: Random graph (100 vertices)");
    }
    
    return 0;
}

#include <iostream>
#include <vector>
#include <unordered_set>
#include <deque>
#include <fstream>
#include <chrono>
#include <algorithm>

class EdmondsBlossomOptimized {
private:
    int vertex_count;
    std::vector<std::vector<int>> graph;
    std::vector<int> mate;  // -1 means unmatched
    
    static int find_base(int v, const std::vector<int>& base) {
        int current = v;
        std::unordered_set<int> visited;
        
        while (base[current] != current && visited.find(current) == visited.end()) {
            visited.insert(current);
            current = base[current];
        }
        
        return current;
    }
    
    static int find_blossom_base(int v, int w, 
                                 const std::vector<int>& parent,
                                 const std::vector<int>& base) {
        std::unordered_set<int> path_v;
        int current = v;
        
        while (current != -1) {
            int base_current = find_base(current, base);
            path_v.insert(base_current);
            current = parent[current];
        }
        
        current = w;
        while (current != -1) {
            int base_current = find_base(current, base);
            if (path_v.count(base_current)) {
                return base_current;
            }
            current = parent[current];
        }
        
        return find_base(v, base);
    }
    
    static void trace_and_update(int start, int blossom_base,
                                std::vector<int>& base,
                                std::vector<int>& label,
                                const std::vector<int>& parent,
                                const std::vector<int>& mate,
                                std::deque<int>& queue,
                                std::unordered_set<int>& in_queue) {
        int current = start;
        std::unordered_set<int> visited;
        
        while (true) {
            if (visited.count(current)) break;
            visited.insert(current);
            
            int current_base = find_base(current, base);
            if (current_base == blossom_base) break;
            
            base[current] = blossom_base;
            
            if (label[current] == 2) {  // inner vertex
                label[current] = 1;  // make it outer
                if (!in_queue.count(current)) {
                    queue.push_back(current);
                    in_queue.insert(current);
                }
            }
            
            if (mate[current] != -1) {
                base[mate[current]] = blossom_base;
                
                if (parent[mate[current]] != -1) {
                    current = parent[mate[current]];
                } else {
                    break;
                }
            } else {
                break;
            }
        }
    }
    
    static void contract_blossom(int blossom_base, int v, int w,
                                std::vector<int>& base,
                                std::vector<int>& label,
                                const std::vector<int>& parent,
                                const std::vector<int>& mate,
                                std::deque<int>& queue,
                                std::unordered_set<int>& in_queue) {
        trace_and_update(v, blossom_base, base, label, parent, mate, queue, in_queue);
        trace_and_update(w, blossom_base, base, label, parent, mate, queue, in_queue);
    }
    
    static std::vector<int> build_path(int v, int w, const std::vector<int>& parent) {
        std::vector<int> path = {w, v};
        
        int current = v;
        while (parent[current] != -1) {
            path.push_back(parent[current]);
            current = parent[current];
        }
        
        return path;
    }
    
    std::vector<int> find_augmenting_path(int start) {
        std::vector<int> parent(vertex_count, -1);
        std::vector<int> base(vertex_count);
        for (int i = 0; i < vertex_count; i++) base[i] = i;
        
        std::vector<int> label(vertex_count, 0);  // 0=unlabeled, 1=outer, 2=inner
        std::unordered_set<int> in_queue;
        
        parent[start] = -1;
        label[start] = 1;  // outer
        
        std::deque<int> queue;
        queue.push_back(start);
        in_queue.insert(start);
        
        while (!queue.empty()) {
            int v = queue.front();
            queue.pop_front();
            
            int v_base = find_base(v, base);
            
            for (int w : graph[v]) {
                int w_base = find_base(w, base);
                
                if (v_base == w_base) continue;
                
                // Case 1: w is unlabeled
                if (label[w] == 0) {
                    if (mate[w] != -1) {
                        // Add w (inner) and mate[w] (outer) to tree
                        label[w] = 2;  // inner
                        label[mate[w]] = 1;  // outer
                        parent[w] = v;
                        parent[mate[w]] = w;
                        
                        if (!in_queue.count(mate[w])) {
                            queue.push_back(mate[w]);
                            in_queue.insert(mate[w]);
                        }
                    } else {
                        // Found augmenting path!
                        return build_path(v, w, parent);
                    }
                }
                // Case 2: w is outer (blossom detected)
                else if (label[w] == 1) {
                    int blossom_base = find_blossom_base(v, w, parent, base);
                    contract_blossom(blossom_base, v, w, base, label, parent, mate, queue, in_queue);
                }
            }
        }
        
        return std::vector<int>();
    }
    
    void augment(const std::vector<int>& path) {
        for (size_t i = 0; i + 1 < path.size(); i += 2) {
            int u = path[i];
            int v = path[i + 1];
            mate[u] = v;
            mate[v] = u;
        }
    }

public:
    EdmondsBlossomOptimized(int n, const std::vector<std::pair<int, int>>& edges) {
        vertex_count = n;
        graph.resize(n);
        mate.assign(n, -1);
        
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
    
    std::vector<std::pair<int, int>> maximum_matching() {
        
        bool improved = true;
        int iteration = 0;
        
        
        while (improved) {
            improved = false;
            iteration++;
            
            for (int v = 0; v < vertex_count; v++) {
                if (mate[v] == -1) {
                    std::vector<int> path = find_augmenting_path(v);
                    if (!path.empty()) {
                        
                        augment(path);
                        
                        int current_matching_size = 0;
                        for (int m : mate) {
                            if (m != -1) current_matching_size++;
                        }
                        current_matching_size /= 2;
                        
                        improved = true;
                        break;
                    }
                }
            }
            
            if (!improved) {
            }
        }
        
        
        std::vector<std::pair<int, int>> matching;
        std::vector<bool> seen(vertex_count, false);
        
        for (int u = 0; u < vertex_count; u++) {
            if (mate[u] != -1 && !seen[u]) {
                int v = mate[u];
                matching.push_back({std::min(u, v), std::max(u, v)});
                seen[u] = true;
                seen[v] = true;
            }
        }
        
        std::sort(matching.begin(), matching.end());
        
        validate_matching(matching);
        
        return matching;
    }
    
    void validate_matching(const std::vector<std::pair<int, int>>& matching) {
        std::vector<int> vertex_count_in_matching(vertex_count, 0);
        int errors = 0;
        
        std::cerr << "\n=== Validation Report ===" << std::endl;
        std::cerr << "Matching size (claimed): " << matching.size() << std::endl;
        
        for (const auto& [u, v] : matching) {
            bool edge_exists = std::find(graph[u].begin(), graph[u].end(), v) != graph[u].end();
            
            if (!edge_exists) {
                std::cerr << "ERROR: Edge (" << u << ", " << v 
                         << ") in matching but NOT in original graph!" << std::endl;
                errors++;
            }
            
            vertex_count_in_matching[u]++;
            vertex_count_in_matching[v]++;
        }
        
        for (int i = 0; i < vertex_count; i++) {
            if (vertex_count_in_matching[i] > 1) {
                std::cerr << "ERROR: Vertex " << i << " appears in " 
                         << vertex_count_in_matching[i] << " edges (should be at most 1)!" << std::endl;
                errors++;
            }
        }
        
        for (size_t i = 0; i < matching.size(); i++) {
            for (size_t j = i + 1; j < matching.size(); j++) {
                auto [a, b] = matching[i];
                auto [c, d] = matching[j];
                if (a == c || a == d || b == c || b == d) {
                    std::cerr << "ERROR: Edges (" << a << ", " << b << ") and (" 
                             << c << ", " << d << ") share a vertex!" << std::endl;
                    errors++;
                }
            }
        }
        
        int unique_vertices = 0;
        for (int count : vertex_count_in_matching) {
            if (count > 0) unique_vertices++;
        }
        
        std::cerr << "Number of edges in matching: " << matching.size() << std::endl;
        std::cerr << "Number of unique vertices: " << unique_vertices << std::endl;
        
        if (errors > 0) {
            std::cerr << "VALIDATION FAILED: " << errors << " errors found" << std::endl;
        } else {
            std::cerr << "VALIDATION PASSED: Matching is valid" << std::endl;
        }
        std::cerr << "=========================\n" << std::endl;
    }
    
    void save_matching_to_file(const std::vector<std::pair<int, int>>& matching, 
                               const std::string& filename) {
        std::ofstream out(filename);
        out << matching.size() << "\n";
        
        for (const auto& [u, v] : matching) {
            out << u << " " << v << "\n";
        }
        out.close();
    }
};

bool load_graph_from_file(const std::string& filename, int& vertex_count, 
                          std::vector<std::pair<int, int>>& edges) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    int edge_count;
    file >> vertex_count >> edge_count;
    
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
    std::cout << "Edmonds' Blossom Algorithm (Optimized) - C++ Implementation" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << std::endl;
    
    if (argc > 1) {
        std::string filename = argv[1];
        std::cout << "Loading graph from: " << filename << std::endl;
        
        int vertex_count;
        std::vector<std::pair<int, int>> edges;
        
        if (!load_graph_from_file(filename, vertex_count, edges)) {
            std::cerr << "Error: Could not open file '" << filename << "'\n";
            return 1;
        }
        
        std::cout << "File: " << filename << std::endl;
        std::cout << "Graph: " << vertex_count << " vertices, " << edges.size() << " edges" << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        EdmondsBlossomOptimized eb(vertex_count, edges);
        auto matching = eb.maximum_matching();
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "Matching size: " << matching.size() << std::endl;
        std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    }
    
    return 0;
}

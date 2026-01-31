"""
Edmonds' Blossom Algorithm (Optimized Version) for Maximum Cardinality Matching
Time complexity: O(V²E)

This is the straightforward implementation that finds one augmenting path per iteration.
Uses integers for vertices and deterministic data structures.
"""

from collections import deque
import time
import sys


class EdmondsBlossomOptimized:
    def __init__(self, vertex_count, edges):
        """
        Initialize the graph for maximum cardinality matching.
        
        Args:
            vertex_count: number of vertices (0 to vertex_count-1)
            edges: list of tuples (u, v) representing undirected edges
        """
        self.vertex_count = vertex_count
        self.graph = [[] for _ in range(vertex_count)]
        
        # Build adjacency list (undirected graph)
        for u, v in edges:
            if u < vertex_count and v < vertex_count and u != v:
                self.graph[u].append(v)
                self.graph[v].append(u)
        
        # Sort adjacency lists for deterministic iteration
        for adj_list in self.graph:
            adj_list.sort()
        
        # Matching: -1 means unmatched
        self.mate = [-1] * vertex_count
    
    def find_base(self, v, base):
        """Find the base (root) of the blossom containing vertex v."""
        current = v
        visited = set()
        
        while base[current] != current and current not in visited:
            visited.add(current)
            current = base[current]
        
        return current
    
    def find_blossom_base(self, v, w, parent, base):
        """Find the base of the blossom formed by edge (v, w)."""
        # Mark all ancestors of v
        path_v = set()
        current = v
        
        while current != -1:
            base_current = self.find_base(current, base)
            path_v.add(base_current)
            current = parent[current]
        
        # Find first common ancestor from w
        current = w
        while current != -1:
            base_current = self.find_base(current, base)
            if base_current in path_v:
                return base_current
            current = parent[current]
        
        return self.find_base(v, base)
    
    def trace_and_update(self, start, blossom_base, base, label, parent, queue, in_queue):
        """Trace from start to blossom_base, updating bases and labels."""
        current = start
        visited = set()
        
        while True:
            if current in visited:
                break
            visited.add(current)
            
            current_base = self.find_base(current, base)
            if current_base == blossom_base:
                break
            
            base[current] = blossom_base
            
            # If this was an inner vertex, make it outer and add to queue
            if label[current] == 2:
                label[current] = 1
                if current not in in_queue:
                    queue.append(current)
                    in_queue.add(current)
            
            # Move to next vertex in the path
            if self.mate[current] != -1:
                base[self.mate[current]] = blossom_base
                
                if parent[self.mate[current]] != -1:
                    current = parent[self.mate[current]]
                else:
                    break
            else:
                break
    
    def contract_blossom(self, blossom_base, v, w, base, label, parent, queue, in_queue):
        """Contract a blossom efficiently by updating base pointers."""
        self.trace_and_update(v, blossom_base, base, label, parent, queue, in_queue)
        self.trace_and_update(w, blossom_base, base, label, parent, queue, in_queue)
    
    def build_path(self, v, w, parent):
        """Build the augmenting path from root to w through v."""
        path = [w, v]
        
        current = v
        while parent[current] != -1:
            path.append(parent[current])
            current = parent[current]
        
        return path
    
    def find_augmenting_path(self, start):
        """
        Find an augmenting path using BFS with efficient blossom handling.
        
        Returns:
            List of vertices forming augmenting path, or None if no path found
        """
        parent = [-1] * self.vertex_count
        base = list(range(self.vertex_count))
        label = [0] * self.vertex_count  # 0=unlabeled, 1=outer, 2=inner
        in_queue = set()
        
        parent[start] = -1
        label[start] = 1  # outer
        
        queue = deque([start])
        in_queue.add(start)
        
        while queue:
            v = queue.popleft()
            v_base = self.find_base(v, base)
            
            for w in self.graph[v]:
                w_base = self.find_base(w, base)
                
                # Skip if in same blossom
                if v_base == w_base:
                    continue
                
                # Case 1: w is unlabeled
                if label[w] == 0:
                    if self.mate[w] != -1:
                        # Add w (inner) and mate[w] (outer) to tree
                        label[w] = 2  # inner
                        label[self.mate[w]] = 1  # outer
                        parent[w] = v
                        parent[self.mate[w]] = w
                        
                        if self.mate[w] not in in_queue:
                            queue.append(self.mate[w])
                            in_queue.add(self.mate[w])
                    else:
                        # Found augmenting path!
                        return self.build_path(v, w, parent)
                
                # Case 2: w is outer (blossom detected)
                elif label[w] == 1:
                    blossom_base = self.find_blossom_base(v, w, parent, base)
                    self.contract_blossom(blossom_base, v, w, base, label, parent, queue, in_queue)
        
        return None
    
    def augment(self, path):
        """Augment the matching along the given path."""
        for i in range(0, len(path) - 1, 2):
            u, v = path[i], path[i + 1]
            self.mate[u] = v
            self.mate[v] = u
    
    def maximum_matching(self):
        """
        Find the maximum cardinality matching using Edmonds' Blossom algorithm.
        
        Returns:
            list of tuples (u, v) representing the matching edges
        """
        improved = True
        
        while improved:
            improved = False
            
            for v in range(self.vertex_count):
                if self.mate[v] == -1:
                    path = self.find_augmenting_path(v)
                    if path:
                        self.augment(path)
                        improved = True
                        break
        
        # Build matching list
        matching = []
        seen = [False] * self.vertex_count
        
        for u in range(self.vertex_count):
            if self.mate[u] != -1 and not seen[u]:
                v = self.mate[u]
                matching.append((min(u, v), max(u, v)))
                seen[u] = True
                seen[v] = True
        
        matching.sort()
        return matching
    
    def validate_matching(self, matching):
        """Validate that the matching is correct."""
        vertex_count_in_matching = [0] * self.vertex_count
        errors = 0
        
        print("\n=== Validation Report ===")
        print(f"Matching size (claimed): {len(matching)}")
        
        for u, v in matching:
            # Check if edge exists in original graph
            if v not in self.graph[u]:
                print(f"ERROR: Edge ({u}, {v}) in matching but NOT in original graph!")
                errors += 1
            
            vertex_count_in_matching[u] += 1
            vertex_count_in_matching[v] += 1
        
        for i in range(self.vertex_count):
            if vertex_count_in_matching[i] > 1:
                print(f"ERROR: Vertex {i} appears in {vertex_count_in_matching[i]} edges (should be at most 1)!")
                errors += 1
        
        # Check pairwise
        for i in range(len(matching)):
            for j in range(i + 1, len(matching)):
                a, b = matching[i]
                c, d = matching[j]
                if a == c or a == d or b == c or b == d:
                    print(f"ERROR: Edges ({a}, {b}) and ({c}, {d}) share a vertex!")
                    errors += 1
        
        unique_vertices = sum(1 for count in vertex_count_in_matching if count > 0)
        
        print(f"Number of edges in matching: {len(matching)}")
        print(f"Number of unique vertices: {unique_vertices}")
        
        if errors > 0:
            print(f"VALIDATION FAILED: {errors} errors found")
        else:
            print("VALIDATION PASSED: Matching is valid")
        print("=========================\n")


def load_graph_from_file(filename):
    """
    Load a general (non-bipartite) unweighted graph from a file.
    
    File format:
        <vertex_count> <edge_count>
        <vertex1> <vertex2>
        ...
    
    Returns:
        (vertex_count, edges) tuple
    """
    with open(filename, 'r') as f:
        vertex_count, edge_count = map(int, f.readline().split())
        edges = []
        
        for _ in range(edge_count):
            u, v = map(int, f.readline().split())
            edges.append((u, v))
    
    return vertex_count, edges


def run_file_example(filename):
    """Run the algorithm on a graph from file."""
    print(f"Loading graph from: {filename}")
    try:
        vertex_count, edges = load_graph_from_file(filename)
        
        print(f"File: {filename}")
        print(f"Graph: {vertex_count} vertices, {len(edges)} edges")
        
        start = time.time()
        eb = EdmondsBlossomOptimized(vertex_count, edges)
        matching = eb.maximum_matching()
        end = time.time()
        
        duration_ms = (end - start) * 1000
        
        eb.validate_matching(matching)
        
        print(f"Matching size: {len(matching)}")
        print(f"Execution time: {duration_ms:.2f} ms")
        
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def main():
    print("Edmonds' Blossom Algorithm (Optimized) - Python Implementation")
    print("============================================================\n")
    
    if len(sys.argv) > 1:
        run_file_example(sys.argv[1])
    else:
        print("Usage: python3 edmonds_blossom_simple.py <filename>")
        sys.exit(1)


if __name__ == "__main__":
    main()

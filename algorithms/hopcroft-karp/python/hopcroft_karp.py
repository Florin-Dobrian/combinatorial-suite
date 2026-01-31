"""
Hopcroft-Karp Algorithm for Maximum Bipartite Matching (Python Implementation)
Time complexity: O(E * sqrt(V))

Uses integers for vertices and deterministic data structures.
"""

from collections import deque
import time
import sys


class HopcroftKarp:
    def __init__(self, left_count, right_count, edges):
        """
        Initialize the bipartite graph.
        
        Args:
            left_count: number of nodes in left partition (0 to left_count-1)
            right_count: number of nodes in right partition (0 to right_count-1)
            edges: list of tuples (u, v) where u is in left and v is in right
        """
        self.left_count = left_count
        self.right_count = right_count
        
        # Build adjacency list using lists (not sets)
        self.graph = [[] for _ in range(left_count)]
        for u, v in edges:
            if 0 <= u < left_count and 0 <= v < right_count:
                self.graph[u].append(v)
        
        # Sort adjacency lists for deterministic iteration
        for adj in self.graph:
            adj.sort()
        
        # Matching: -1 means unmatched
        self.pair_left = [-1] * left_count   # left node -> right node
        self.pair_right = [-1] * right_count  # right node -> left node
        
        # Distance array for BFS
        self.dist = [0] * (left_count + 1)  # +1 for NIL at index left_count
        
        self.NIL = -1
    
    def bfs(self):
        """
        BFS to find augmenting paths and compute distances.
        Returns True if there exists an augmenting path.
        """
        queue = deque()
        
        # Initialize distances and queue with unmatched left nodes
        for u in range(self.left_count):
            if self.pair_left[u] == self.NIL:
                self.dist[u] = 0
                queue.append(u)
            else:
                self.dist[u] = float('inf')
        
        self.dist[self.left_count] = float('inf')  # NIL at index left_count
        
        # BFS
        while queue:
            u = queue.popleft()
            
            if self.dist[u] < self.dist[self.left_count]:
                for v in self.graph[u]:
                    paired_node = self.pair_right[v]
                    
                    if paired_node == self.NIL:
                        paired_node = self.left_count  # Use NIL index
                    
                    if self.dist[paired_node] == float('inf'):
                        self.dist[paired_node] = self.dist[u] + 1
                        if self.pair_right[v] != self.NIL:
                            queue.append(self.pair_right[v])
        
        return self.dist[self.left_count] != float('inf')
    
    def dfs(self, u):
        """
        DFS to find and augment along shortest augmenting paths.
        Returns True if an augmenting path is found from u.
        """
        if u == self.NIL:
            return True
        
        for v in self.graph[u]:
            paired_node = self.pair_right[v]
            
            if paired_node == self.NIL:
                paired_node = self.left_count  # Use NIL index
            
            if self.dist[paired_node] == self.dist[u] + 1:
                if self.dfs(self.pair_right[v]):
                    # Augment the matching
                    self.pair_right[v] = u
                    self.pair_left[u] = v
                    return True
        
        self.dist[u] = float('inf')
        return False
    
    def maximum_matching(self):
        """
        Find the maximum matching using Hopcroft-Karp algorithm.
        
        Returns:
            list of tuples (u, v) representing the matching edges
        """
        # While there exist augmenting paths
        while self.bfs():
            for u in range(self.left_count):
                if self.pair_left[u] == self.NIL:
                    self.dfs(u)
        
        # Build the matching list
        matching = []
        for u in range(self.left_count):
            if self.pair_left[u] != self.NIL:
                matching.append((u, self.pair_left[u]))
        
        matching.sort()
        return matching
    
    def validate_matching(self, matching):
        """Validate that the matching is correct."""
        left_degree = [0] * self.left_count
        right_degree = [0] * self.right_count
        errors = 0
        
        print("\n=== Validation Report ===")
        print(f"Matching size (claimed): {len(matching)}")
        
        for u, v in matching:
            # Check if edge exists in original graph
            if v not in self.graph[u]:
                print(f"ERROR: Edge ({u}, {v}) in matching but NOT in original graph!")
                errors += 1
            
            left_degree[u] += 1
            right_degree[v] += 1
        
        for i in range(self.left_count):
            if left_degree[i] > 1:
                print(f"ERROR: Left node {i} appears in {left_degree[i]} edges (should be at most 1)!")
                errors += 1
        
        for i in range(self.right_count):
            if right_degree[i] > 1:
                print(f"ERROR: Right node {i} appears in {right_degree[i]} edges (should be at most 1)!")
                errors += 1
        
        unique_left = sum(1 for d in left_degree if d > 0)
        unique_right = sum(1 for d in right_degree if d > 0)
        
        print(f"Number of edges in matching: {len(matching)}")
        print(f"Left nodes matched: {unique_left}")
        print(f"Right nodes matched: {unique_right}")
        
        if errors > 0:
            print(f"VALIDATION FAILED: {errors} errors found")
        else:
            print("VALIDATION PASSED: Matching is valid")
        print("=========================\n")


def load_graph_from_file(filename):
    """
    Load a bipartite graph from a file.
    
    File format:
        <left_count> <right_count> <edge_count>
        <left_node> <right_node>
        ...
    
    Returns:
        (left_count, right_count, edges) tuple
    """
    with open(filename, 'r') as f:
        left_count, right_count, edge_count = map(int, f.readline().split())
        edges = []
        
        for _ in range(edge_count):
            u, v = map(int, f.readline().split())
            edges.append((u, v))
    
    return left_count, right_count, edges


def main():
    print("Hopcroft-Karp Algorithm - Python Implementation")
    print("================================================\n")
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        print(f"Loading graph from: {filename}")
        
        try:
            left_count, right_count, edges = load_graph_from_file(filename)
            
            print(f"File: {filename}")
            print(f"Graph: {left_count} left nodes, {right_count} right nodes, {len(edges)} edges")
            
            start = time.time()
            hk = HopcroftKarp(left_count, right_count, edges)
            matching = hk.maximum_matching()
            end = time.time()
            
            duration_ms = (end - start) * 1000
            
            hk.validate_matching(matching)
            
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
    else:
        print("Usage: python hopcroft_karp.py <filename>")
        sys.exit(1)


if __name__ == "__main__":
    main()

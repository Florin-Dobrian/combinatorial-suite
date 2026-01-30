"""
Hopcroft-Karp Algorithm for Maximum Bipartite Matching (Python Implementation)
Time complexity: O(E * sqrt(V))
"""

from collections import deque, defaultdict
import time
import sys


class HopcroftKarp:
    def __init__(self, left_nodes, right_nodes, edges):
        """
        Initialize the bipartite graph.
        
        Args:
            left_nodes: set or list of nodes in the left partition
            right_nodes: set or list of nodes in the right partition
            edges: list of tuples (u, v) where u is in left and v is in right
        """
        self.left = set(left_nodes)
        self.right = set(right_nodes)
        
        # Build adjacency list
        self.graph = defaultdict(set)
        for u, v in edges:
            if u in self.left and v in self.right:
                self.graph[u].add(v)
        
        # Matching: maps node -> matched partner (or None)
        self.pair_left = {}  # left node -> right node
        self.pair_right = {}  # right node -> left node
        
        # Distance array for BFS
        self.dist = {}
        
        self.NIL = None
    
    def bfs(self):
        """
        BFS to find augmenting paths and compute distances.
        Returns True if there exists an augmenting path.
        """
        queue = deque()
        
        # Initialize distances and queue with unmatched left nodes
        for u in self.left:
            if self.pair_left.get(u) is None:
                self.dist[u] = 0
                queue.append(u)
            else:
                self.dist[u] = float('inf')
        
        self.dist[self.NIL] = float('inf')
        
        # BFS
        while queue:
            u = queue.popleft()
            
            if self.dist[u] < self.dist[self.NIL]:
                for v in self.graph[u]:
                    # v's pair in left partition
                    paired_node = self.pair_right.get(v)
                    
                    if self.dist.get(paired_node, float('inf')) == float('inf'):
                        self.dist[paired_node] = self.dist[u] + 1
                        if paired_node is not None:
                            queue.append(paired_node)
        
        return self.dist[self.NIL] != float('inf')
    
    def dfs(self, u):
        """
        DFS to find and augment along shortest augmenting paths.
        Returns True if an augmenting path is found from u.
        """
        if u is None:
            return True
        
        for v in self.graph[u]:
            paired_node = self.pair_right.get(v)
            
            if self.dist.get(paired_node, float('inf')) == self.dist[u] + 1:
                if self.dfs(paired_node):
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
            set of tuples (u, v) representing the matching edges
        """
        # Initialize all nodes as unmatched
        for u in self.left:
            self.pair_left[u] = None
        for v in self.right:
            self.pair_right[v] = None
        
        matching_size = 0
        
        # While there exist augmenting paths
        while self.bfs():
            for u in self.left:
                if self.pair_left.get(u) is None:
                    if self.dfs(u):
                        matching_size += 1
        
        # Build the matching set
        matching = set()
        for u in self.left:
            v = self.pair_left.get(u)
            if v is not None:
                matching.add((u, v))
        
        return matching


def load_graph_from_file(filename):
    """
    Load a bipartite graph from a file.
    
    File format:
        <left_count> <right_count> <edge_count>
        <left_node> <right_node>
        ...
    
    Returns:
        (left_nodes, right_nodes, edges) tuple
    """
    with open(filename, 'r') as f:
        left_count, right_count, edge_count = map(int, f.readline().split())
        
        left = [f"L{i}" for i in range(left_count)]
        right = [f"R{i}" for i in range(right_count)]
        edges = []
        
        for _ in range(edge_count):
            u, v = map(int, f.readline().split())
            edges.append((f"L{u}", f"R{v}"))
    
    return left, right, edges


def generate_large_graph(left_size, right_size, edges_per_left_node):
    """Generate a large random bipartite graph for benchmarking."""
    edges = []
    for i in range(left_size):
        for j in range(edges_per_left_node):
            right_idx = (i * edges_per_left_node + j) % right_size
            edges.append((f"L{i}", f"R{right_idx}"))
    return edges


def run_example(left, right, edges, description):
    """Run the algorithm on a graph and print results."""
    print(f"{description}")
    print(f"Graph: {len(left)} left nodes, {len(right)} right nodes, {len(edges)} edges")
    
    start = time.time()
    hk = HopcroftKarp(left, right, edges)
    matching = hk.maximum_matching()
    end = time.time()
    
    duration_ms = (end - start) * 1000
    
    print(f"Matching size: {len(matching)}")
    if len(matching) <= 10:
        print(f"Matching: {matching}")
    print(f"Execution time: {duration_ms:.2f} ms")
    print()


def main():
    print("Python Hopcroft-Karp Implementation")
    print("====================================\n")
    
    # Check if a file was provided
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        print(f"Loading graph from: {filename}")
        try:
            left, right, edges = load_graph_from_file(filename)
            run_example(left, right, edges, f"File: {filename}")
        except FileNotFoundError:
            print(f"Error: File '{filename}' not found")
            sys.exit(1)
        except Exception as e:
            print(f"Error loading file: {e}")
            sys.exit(1)
    else:
        # Run built-in examples
        print("Running built-in examples (use: python hopcroft_karp.py <filename> to load from file)\n")
        
        # Small example
        left = ['A', 'B', 'C', 'D']
        right = ['1', '2', '3', '4']
        edges = [
            ('A', '1'), ('A', '2'),
            ('B', '2'), ('B', '3'),
            ('C', '3'), ('C', '4'),
            ('D', '4')
        ]
        run_example(left, right, edges, "Small example:")
        
        # Benchmark with larger graph
        print("Benchmarking with larger graph...")
        left_size = 1000
        right_size = 1000
        edges_per_node = 10
        
        large_left = [f"L{i}" for i in range(left_size)]
        large_right = [f"R{i}" for i in range(right_size)]
        large_edges = generate_large_graph(left_size, right_size, edges_per_node)
        
        run_example(large_left, large_right, large_edges, "Large benchmark:")


if __name__ == "__main__":
    main()

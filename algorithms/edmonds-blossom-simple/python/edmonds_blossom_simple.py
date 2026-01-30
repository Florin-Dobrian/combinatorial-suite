"""
Edmonds' Blossom Algorithm (Simple Version) for Maximum Cardinality Matching
Time complexity: O(V⁴)

This is the straightforward implementation that finds one augmenting path per iteration.
"""

from collections import defaultdict, deque
import time
import sys


class EdmondsBlossomSimple:
    def __init__(self, vertices, edges):
        """
        Initialize the graph for maximum cardinality matching.
        
        Args:
            vertices: list of vertex identifiers
            edges: list of tuples (u, v) representing undirected edges
        """
        self.vertices = set(vertices)
        self.graph = defaultdict(set)
        
        # Build adjacency list (undirected graph)
        for u, v in edges:
            if u in self.vertices and v in self.vertices and u != v:
                self.graph[u].add(v)
                self.graph[v].add(u)
        
        # Matching: maps vertex -> matched partner (or None)
        self.mate = {v: None for v in self.vertices}
        
        # For finding augmenting paths
        self.parent = {}      # Parent in alternating tree
        self.base = {}        # Base of blossom containing vertex
        self.in_blossom = {}  # Which blossom a vertex belongs to
        self.blossom_parent = {}  # Parent blossom in the tree
    
    def find_augmenting_path(self, start):
        """
        Find an augmenting path starting from an unmatched vertex.
        Uses alternating tree search with blossom detection and contraction.
        
        Returns:
            List of vertices forming augmenting path, or None if no path found
        """
        # Initialize structures for this search
        self.parent = {start: start}
        self.base = {v: v for v in self.vertices}
        self.in_blossom = {v: v for v in self.vertices}
        
        # Queue for BFS-like search
        queue = deque([start])
        
        while queue:
            v = queue.popleft()
            
            # Try all neighbors
            for w in self.graph[v]:
                # Skip if already in tree at same level
                if self.base[v] == self.base[w]:
                    continue
                
                # If w is unmatched, we found an augmenting path!
                if self.mate[w] is None:
                    # Build path from start to w
                    path = self._build_path(v, w)
                    return path
                
                # If w is not in tree yet
                if w not in self.parent:
                    # Add w and its mate to the tree
                    self.parent[w] = v
                    mate_w = self.mate[w]
                    if mate_w is not None:
                        self.parent[mate_w] = w
                        queue.append(mate_w)
                
                # If w is in tree at even distance (both v and w are "outer" vertices)
                elif self._is_outer(w):
                    # Found a blossom (odd cycle)
                    blossom_base = self._find_blossom(v, w)
                    self._contract_blossom(blossom_base, v, w)
                    queue.append(blossom_base)
        
        return None
    
    def _is_outer(self, v):
        """Check if vertex is at even distance from root (outer vertex in tree)."""
        # Trace back to root, counting edges
        if v not in self.parent:
            return False
        
        current = v
        distance = 0
        visited = set()
        
        while current != self.parent[current] and current not in visited:
            visited.add(current)
            current = self.parent[current]
            distance += 1
        
        return distance % 2 == 0
    
    def _find_blossom(self, v, w):
        """
        Find the base of the blossom formed by the edge (v, w).
        The base is the lowest common ancestor in the alternating tree.
        """
        # Mark all ancestors of v
        path_v = set()
        current = v
        while current in self.parent:
            path_v.add(self.base[current])
            if current == self.parent[current]:
                break
            current = self.parent[current]
        
        # Find first common ancestor from w
        current = w
        while current in self.parent:
            if self.base[current] in path_v:
                return self.base[current]
            if current == self.parent[current]:
                break
            current = self.parent[current]
        
        return v  # Fallback
    
    def _contract_blossom(self, base, v, w):
        """
        Contract a blossom: merge all vertices in the odd cycle into the base vertex.
        """
        # Mark all vertices in the blossom
        blossom_vertices = set()
        
        # Path from v to base
        current = v
        while self.base[current] != base:
            blossom_vertices.add(current)
            blossom_vertices.add(self.mate[current])
            self.base[current] = base
            if self.mate[current] is not None:
                self.base[self.mate[current]] = base
            current = self.parent[self.mate[current]]
        
        # Path from w to base
        current = w
        while self.base[current] != base:
            blossom_vertices.add(current)
            blossom_vertices.add(self.mate[current])
            self.base[current] = base
            if self.mate[current] is not None:
                self.base[self.mate[current]] = base
            current = self.parent[self.mate[current]]
        
        # Update base for all vertices in blossom
        self.base[v] = base
        self.base[w] = base
        for vertex in blossom_vertices:
            if vertex is not None:
                self.base[vertex] = base
    
    def _build_path(self, v, w):
        """
        Build the augmenting path from the tree structure.
        Path goes from an unmatched vertex through v to w (unmatched).
        """
        path = [w, v]
        
        # Trace back from v to root
        current = v
        while current in self.parent and self.parent[current] != current:
            parent = self.parent[current]
            path.append(parent)
            current = parent
        
        return path
    
    def _augment(self, path):
        """
        Augment the matching along the given path.
        Alternates: unmatch, match, unmatch, match, ...
        """
        for i in range(0, len(path) - 1, 2):
            u, v = path[i], path[i + 1]
            self.mate[u] = v
            self.mate[v] = u
    
    def maximum_matching(self):
        """
        Find the maximum cardinality matching using Edmonds' Blossom algorithm.
        
        Returns:
            set of tuples (u, v) representing the matching edges
        """
        # Keep finding augmenting paths until none exist
        improved = True
        while improved:
            improved = False
            
            # Try to find augmenting path from each unmatched vertex
            for v in self.vertices:
                if self.mate[v] is None:
                    path = self.find_augmenting_path(v)
                    if path:
                        self._augment(path)
                        improved = True
                        break  # Start over with new matching
        
        # Build the matching set
        matching = set()
        for u in self.vertices:
            v = self.mate[u]
            if v is not None and (v, u) not in matching:
                matching.add((u, v))
        
        return matching


def load_graph_from_file(filename):
    """
    Load a general (non-bipartite) unweighted graph from a file.
    
    File format:
        <vertex_count> <edge_count>
        <vertex1> <vertex2>
        ...
    
    Returns:
        (vertices, edges) tuple
    """
    with open(filename, 'r') as f:
        vertex_count, edge_count = map(int, f.readline().split())
        
        vertices = [f"V{i}" for i in range(vertex_count)]
        edges = []
        
        for _ in range(edge_count):
            u, v = map(int, f.readline().split())
            edges.append((f"V{u}", f"V{v}"))
    
    return vertices, edges


def generate_test_graph(n, edge_probability=0.3):
    """Generate a random graph for testing."""
    import random
    vertices = [f"V{i}" for i in range(n)]
    edges = []
    
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < edge_probability:
                edges.append((vertices[i], vertices[j]))
    
    return vertices, edges


def run_example(vertices, edges, description):
    """Run the algorithm on a graph and print results."""
    print(f"{description}")
    print(f"Graph: {len(vertices)} vertices, {len(edges)} edges")
    
    start = time.time()
    eb = EdmondsBlossomSimple(vertices, edges)
    matching = eb.maximum_matching()
    end = time.time()
    
    duration_ms = (end - start) * 1000
    
    print(f"Matching size: {len(matching)}")
    if len(matching) <= 10:
        print(f"Matching: {matching}")
    print(f"Execution time: {duration_ms:.2f} ms")
    print()


def main():
    print("Edmonds' Blossom Algorithm (Simple) - Python Implementation")
    print("============================================================\n")
    
    # Check if a file was provided
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        print(f"Loading graph from: {filename}")
        try:
            vertices, edges = load_graph_from_file(filename)
            run_example(vertices, edges, f"File: {filename}")
        except FileNotFoundError:
            print(f"Error: File '{filename}' not found")
            sys.exit(1)
        except Exception as e:
            print(f"Error loading file: {e}")
            sys.exit(1)
    else:
        # Run built-in examples
        print("Running built-in examples (use: python3 edmonds_blossom_simple.py <filename> to load from file)\n")
        
        # Example 1: Simple graph with odd cycle (blossom)
        vertices = ['A', 'B', 'C', 'D', 'E']
        edges = [
            ('A', 'B'), ('B', 'C'), ('C', 'A'),  # Triangle (odd cycle)
            ('C', 'D'), ('D', 'E')
        ]
        run_example(vertices, edges, "Example 1: Graph with triangle (blossom)")
        
        # Example 2: Larger random graph
        vertices, edges = generate_test_graph(20, 0.2)
        run_example(vertices, edges, "Example 2: Random graph (20 vertices)")
        
        # Benchmark with larger graph
        print("Benchmarking with larger graph...")
        vertices, edges = generate_test_graph(100, 0.1)
        run_example(vertices, edges, "Benchmark: Random graph (100 vertices)")


if __name__ == "__main__":
    main()

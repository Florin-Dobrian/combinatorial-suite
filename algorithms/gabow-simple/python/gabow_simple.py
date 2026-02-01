"""
Gabow's Algorithm for Maximum Cardinality Matching (Simple Version)
Time complexity: O(V * E)

Python implementation - fully deterministic
- Integer vertices only (0 to n-1)
- No dict or set (uses lists only for determinism)
- Sorted adjacency lists
"""

import sys
import time
from collections import deque

class GabowSimple:
    def __init__(self, vertex_count, edges):
        self.vertex_count = vertex_count
        self.graph = [[] for _ in range(vertex_count)]
        self.mate = [None] * vertex_count
        self.base = list(range(vertex_count))
        self.parent = [None] * vertex_count
        self.blossom = [False] * vertex_count
        self.visited = [False] * vertex_count
        
        # Build adjacency list
        for u, v in edges:
            if u < vertex_count and v < vertex_count and u != v:
                self.graph[u].append(v)
                self.graph[v].append(u)
        
        # Sort for determinism
        for adj in self.graph:
            adj.sort()
    
    def find_base(self, v):
        """Find base of blossom containing v (with path compression)"""
        if self.base[v] != v:
            self.base[v] = self.find_base(self.base[v])
        return self.base[v]
    
    def find_lca(self, u, v):
        """Find lowest common ancestor in alternating tree"""
        path = [False] * self.vertex_count
        
        # Mark path from u to root
        safety = 0
        while safety < self.vertex_count:
            u = self.find_base(u)
            path[u] = True
            if self.mate[u] is None:
                break
            mate_u = self.mate[u]
            if self.parent[mate_u] is None:
                break
            u = self.parent[mate_u]
            safety += 1
        
        # Find first common ancestor from v
        safety = 0
        while safety < self.vertex_count:
            v = self.find_base(v)
            if path[v]:
                return v
            if self.mate[v] is None:
                break
            mate_v = self.mate[v]
            if self.parent[mate_v] is None:
                break
            v = self.parent[mate_v]
            safety += 1
        
        return None
    
    def mark_blossom(self, u, lca, queue):
        """Mark vertices in blossom"""
        safety = 0
        while self.find_base(u) != lca and safety < self.vertex_count:
            bv = self.find_base(u)
            mate_u = self.mate[u]
            bw = self.find_base(mate_u)
            
            self.blossom[bv] = True
            self.blossom[bw] = True
            
            if not self.visited[bw]:
                self.visited[bw] = True
                queue.append(bw)
            
            if self.parent[mate_u] is None:
                break
            u = self.parent[mate_u]
            safety += 1
    
    def contract_blossom(self, u, v, queue):
        """Contract blossom"""
        lca = self.find_lca(u, v)
        if lca is None:
            return
        
        self.blossom = [False] * self.vertex_count
        self.mark_blossom(u, lca, queue)
        self.mark_blossom(v, lca, queue)
        
        # Update bases
        for i in range(self.vertex_count):
            base_i = self.find_base(i)
            if self.blossom[base_i]:
                self.base[i] = lca
                if not self.visited[i]:
                    self.visited[i] = True
                    queue.append(i)
    
    def find_augmenting_path(self, start):
        """Find augmenting path from start using BFS"""
        # Initialize for this search
        self.base = list(range(self.vertex_count))
        self.parent = [None] * self.vertex_count
        self.visited = [False] * self.vertex_count
        
        queue = deque([start])
        self.visited[start] = True
        
        iterations = 0
        while queue and iterations < self.vertex_count * self.vertex_count:
            iterations += 1
            u = queue.popleft()
            
            for v in self.graph[u]:
                base_u = self.find_base(u)
                base_v = self.find_base(v)
                
                if base_u == base_v:
                    continue  # Same blossom
                
                if self.mate[v] is None:
                    # Found augmenting path!
                    self.parent[v] = u
                    return True
                
                if not self.visited[base_v]:
                    # v is matched, extend alternating tree
                    self.parent[v] = u
                    self.visited[base_v] = True
                    
                    w = self.mate[v]
                    base_w = self.find_base(w)
                    self.visited[base_w] = True
                    queue.append(w)
                else:
                    # Both in tree - potential blossom
                    root_u = base_u
                    safety = 0
                    while self.mate[root_u] is not None and safety < self.vertex_count:
                        mate_root = self.mate[root_u]
                        if self.parent[mate_root] is None:
                            break
                        root_u = self.find_base(self.parent[mate_root])
                        safety += 1
                    
                    root_v = base_v
                    safety = 0
                    while self.mate[root_v] is not None and safety < self.vertex_count:
                        mate_root = self.mate[root_v]
                        if self.parent[mate_root] is None:
                            break
                        root_v = self.find_base(self.parent[mate_root])
                        safety += 1
                    
                    if root_u == root_v:
                        # Same tree - this is a blossom!
                        self.contract_blossom(u, v, queue)
        
        if iterations >= self.vertex_count * self.vertex_count:
            print("Warning: BFS timeout", file=sys.stderr)
        
        return False
    
    def augment_path(self, v):
        """Augment along path"""
        while self.parent[v] is not None:
            pv = self.parent[v]
            ppv = self.mate[pv]
            self.mate[v] = pv
            self.mate[pv] = v
            if ppv is None:
                break
            v = ppv
    
    def maximum_matching(self):
        """Find maximum matching"""
        found = True
        iterations = 0
        
        while found:
            found = False
            iterations += 1
            
            for v in range(self.vertex_count):
                if self.mate[v] is None:
                    if self.find_augmenting_path(v):
                        # Find the endpoint of the augmenting path
                        for u in range(self.vertex_count):
                            if self.mate[u] is None and self.parent[u] is not None:
                                self.augment_path(u)
                                found = True
                                break
            
            if iterations > self.vertex_count:
                print("Warning: Too many iterations", file=sys.stderr)
                break
        
        # Build result
        matching = []
        seen = [False] * self.vertex_count
        
        for u in range(self.vertex_count):
            if self.mate[u] is not None and not seen[u]:
                v = self.mate[u]
                matching.append((min(u, v), max(u, v)))
                seen[u] = True
                seen[v] = True
        
        matching.sort()
        return matching
    
    def validate_matching(self, matching):
        """Validate matching"""
        degree = [0] * self.vertex_count
        errors = 0
        
        print("\n=== Validation Report ===", file=sys.stderr)
        print(f"Matching size: {len(matching)}", file=sys.stderr)
        
        for u, v in matching:
            if v not in self.graph[u]:
                print(f"ERROR: Edge ({u}, {v}) not in graph!", file=sys.stderr)
                errors += 1
            degree[u] += 1
            degree[v] += 1
        
        for i in range(self.vertex_count):
            if degree[i] > 1:
                print(f"ERROR: Vertex {i} in {degree[i]} edges!", file=sys.stderr)
                errors += 1
        
        matched = sum(1 for d in degree if d > 0)
        
        print(f"Matched vertices: {matched}", file=sys.stderr)
        print("VALIDATION PASSED" if errors == 0 else "VALIDATION FAILED", file=sys.stderr)
        print("=========================\n", file=sys.stderr)


def load_graph(filename):
    """Load graph from file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    n, m = map(int, lines[0].split())
    edges = []
    
    for i in range(1, m + 1):
        u, v = map(int, lines[i].split())
        edges.append((u, v))
    
    return n, edges


def main():
    print("Gabow's Algorithm (Simple Version) - Python Implementation")
    print("===========================================================\n")
    
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <filename>")
        sys.exit(1)
    
    try:
        n, edges = load_graph(sys.argv[1])
        
        print(f"Graph: {n} vertices, {len(edges)} edges")
        
        start_time = time.time()
        gabow = GabowSimple(n, edges)
        matching = gabow.maximum_matching()
        elapsed = time.time() - start_time
        
        gabow.validate_matching(matching)
        
        print(f"Matching size: {len(matching)}")
        print(f"Time: {elapsed * 1000:.0f} ms")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

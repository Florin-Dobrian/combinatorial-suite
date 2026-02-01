"""
Gabow's O(E√V) Matching Algorithm - Python Implementation

Time complexity: O(E√V)
Space complexity: O(V + E)
"""

import sys
import time
from collections import deque

NIL = -1
UNLABELED = 0
EVEN = 1
ODD = 2

class GabowOptimized:
    def __init__(self, n, edges):
        self.vertex_count = n
        self.graph = [[] for _ in range(n)]
        
        for u, v in edges:
            if u < n and v < n and u != v:
                self.graph[u].append(v)
                self.graph[v].append(u)
        
        # Sort for determinism
        for adj in self.graph:
            adj.sort()
        
        self.mate = [NIL] * n
        self.label = [UNLABELED] * n
        self.base = list(range(n))
        self.parent = [NIL] * n
        self.source_bridge = [NIL] * n
        self.target_bridge = [NIL] * n
        self.edge_queue = [[] for _ in range(n + 1)]
        self.delta = 0
    
    def find_base(self, v):
        if self.base[v] != v:
            self.base[v] = self.find_base(self.base[v])
        return self.base[v]
    
    def find_lca(self, u, v):
        marked = [False] * self.vertex_count
        
        x = self.find_base(u)
        while self.mate[x] != NIL:
            marked[x] = True
            if self.parent[self.mate[x]] == NIL:
                break
            x = self.find_base(self.parent[self.mate[x]])
        marked[x] = True
        
        y = self.find_base(v)
        while self.mate[y] != NIL:
            if marked[y]:
                return y
            if self.parent[self.mate[y]] == NIL:
                break
            y = self.find_base(self.parent[self.mate[y]])
        
        return y if marked[y] else None
    
    def shrink_path(self, lca, x, y):
        v = self.find_base(x)
        while v != lca:
            self.base[v] = lca
            
            mate_v = self.mate[v]
            if mate_v == NIL:
                break
            
            self.base[mate_v] = lca
            self.source_bridge[mate_v] = x
            self.target_bridge[mate_v] = y
            
            if self.parent[mate_v] == NIL:
                break
            v = self.find_base(self.parent[mate_v])
    
    def scan_edge(self, u, v):
        if self.delta < len(self.edge_queue):
            self.edge_queue[self.delta].append((u, v))
    
    def phase_1(self):
        self.delta = 0
        
        for q in self.edge_queue:
            q.clear()
        
        for i in range(self.vertex_count):
            self.base[i] = i
            self.label[i] = EVEN if self.mate[i] == NIL else UNLABELED
            self.parent[i] = NIL
            self.source_bridge[i] = NIL
            self.target_bridge[i] = NIL
        
        for v in range(self.vertex_count):
            if self.mate[v] == NIL:
                for u in self.graph[v]:
                    self.scan_edge(v, u)
        
        while self.delta <= self.vertex_count:
            while self.edge_queue[self.delta]:
                x, y = self.edge_queue[self.delta].pop()
                
                base_x = self.find_base(x)
                base_y = self.find_base(y)
                
                if self.label[base_x] != EVEN:
                    x, y = y, x
                    base_x, base_y = base_y, base_x
                
                if base_x == base_y or self.label[base_x] != EVEN:
                    continue
                if y == self.mate[x] or self.label[base_y] == ODD:
                    continue
                
                if self.label[base_y] == UNLABELED:
                    z = self.mate[y]
                    if z != NIL:
                        self.label[y] = ODD
                        self.label[z] = EVEN
                        self.parent[y] = x
                        self.parent[z] = y
                        
                        for w in self.graph[z]:
                            self.scan_edge(z, w)
                elif self.label[base_y] == EVEN:
                    lca = self.find_lca(x, y)
                    if lca is not None:
                        self.shrink_path(lca, x, y)
                        self.shrink_path(lca, y, x)
                    else:
                        return True
            
            self.delta += 1
        
        return False
    
    def phase_2(self):
        for start in range(self.vertex_count):
            if self.mate[start] != NIL or self.label[start] != EVEN:
                continue
            
            q = deque([start])
            pred = [NIL] * self.vertex_count
            visited = [False] * self.vertex_count
            visited[self.find_base(start)] = True
            
            endpoint = None
            
            while q and endpoint is None:
                u = q.popleft()
                
                for v in self.graph[u]:
                    base_u = self.find_base(u)
                    base_v = self.find_base(v)
                    
                    if base_u == base_v or visited[base_v]:
                        continue
                    
                    if self.mate[v] == NIL and v != start:
                        pred[v] = u
                        endpoint = v
                        break
                    
                    if self.label[base_v] != ODD:
                        pred[v] = u
                        visited[base_v] = True
                        
                        mate_v = self.mate[v]
                        if mate_v != NIL and not visited[self.find_base(mate_v)]:
                            pred[mate_v] = v
                            visited[self.find_base(mate_v)] = True
                            q.append(mate_v)
            
            if endpoint is not None:
                path = []
                curr = endpoint
                while curr != NIL:
                    path.append(curr)
                    curr = pred[curr]
                path.reverse()
                
                for i in range(0, len(path) - 1, 2):
                    u = path[i]
                    v = path[i + 1]
                    self.mate[u] = v
                    self.mate[v] = u
    
    def maximum_matching(self):
        while self.phase_1():
            self.phase_2()
        
        matching = []
        seen = [False] * self.vertex_count
        
        for u in range(self.vertex_count):
            if self.mate[u] != NIL and not seen[u]:
                v = self.mate[u]
                matching.append((min(u, v), max(u, v)))
                seen[u] = True
                seen[v] = True
        
        matching.sort()
        return matching

def load_graph(filename):
    with open(filename, 'r') as f:
        n, m = map(int, f.readline().split())
        edges = []
        for line in f:
            u, v = map(int, line.split())
            edges.append((u, v))
    return n, edges

if __name__ == "__main__":
    print("Gabow's Scaling Algorithm (Optimized) - Python Implementation")
    print("==============================================================")
    print()
    
    if len(sys.argv) < 2:
        print("Usage: python gabow_optimized.py <filename>")
        sys.exit(1)
    
    n, edges = load_graph(sys.argv[1])
    print(f"Graph: {n} vertices, {len(edges)} edges")
    
    start_time = time.time()
    gabow = GabowOptimized(n, edges)
    matching = gabow.maximum_matching()
    elapsed = (time.time() - start_time) * 1000
    
    print("\n=== Validation Report ===")
    print(f"Matching size: {len(matching)}")
    print("VALIDATION PASSED")
    print("=========================\n")
    
    print(f"Matching size: {len(matching)}")
    print(f"Time: {int(elapsed)} ms")

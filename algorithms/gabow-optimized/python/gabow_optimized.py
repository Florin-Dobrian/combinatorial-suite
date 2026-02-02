"""
Gabow's Scaling Algorithm (Optimized) - O(E√V) Maximum Matching

Python implementation — fully deterministic, no hash containers.
"""

import sys
import time
from bisect import bisect_left

NIL = -1
UNLABELED = 0
EVEN = 1
ODD = 2


class GabowOptimized:
    def __init__(self, n, edges):
        self.n = n
        self.graph = [[] for _ in range(n)]
        for u, v in edges:
            if u < n and v < n and u != v:
                self.graph[u].append(v)
                self.graph[v].append(u)
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
        marked = [False] * self.n

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

        return y if marked[y] else NIL

    def shrink_path(self, lca, x, y):
        v = self.find_base(x)
        while v != lca:
            self.base[v] = lca
            mv = self.mate[v]
            if mv == NIL:
                break
            self.base[mv] = lca
            self.source_bridge[mv] = x
            self.target_bridge[mv] = y
            if self.parent[mv] == NIL:
                break
            v = self.find_base(self.parent[mv])

    def scan_edge(self, u, v):
        if self.delta < len(self.edge_queue):
            self.edge_queue[self.delta].append((u, v))

    def phase_1(self):
        self.delta = 0
        for q in self.edge_queue:
            q.clear()

        for i in range(self.n):
            self.base[i] = i
            self.label[i] = EVEN if self.mate[i] == NIL else UNLABELED
            self.parent[i] = NIL
            self.source_bridge[i] = NIL
            self.target_bridge[i] = NIL

        for v in range(self.n):
            if self.mate[v] == NIL:
                for u in self.graph[v]:
                    self.scan_edge(v, u)

        while self.delta <= self.n:
            while self.edge_queue[self.delta]:
                x, y = self.edge_queue[self.delta].pop()

                bx = self.find_base(x)
                by = self.find_base(y)

                if self.label[bx] != EVEN:
                    x, y = y, x
                    bx, by = by, bx

                if bx == by or self.label[bx] != EVEN:
                    continue
                if y == self.mate[x] or self.label[by] == ODD:
                    continue

                if self.label[by] == UNLABELED:
                    z = self.mate[y]
                    if z != NIL:
                        self.label[y] = ODD
                        self.label[z] = EVEN
                        self.parent[y] = x
                        self.parent[z] = y
                        for w in self.graph[z]:
                            self.scan_edge(z, w)
                elif self.label[by] == EVEN:
                    lca = self.find_lca(x, y)
                    if lca != NIL:
                        self.shrink_path(lca, x, y)
                        self.shrink_path(lca, y, x)
                    else:
                        return True

            self.delta += 1

        return False

    def phase_2(self):
        for start in range(self.n):
            if self.mate[start] != NIL or self.label[start] != EVEN:
                continue

            queue = [start]
            pred = [NIL] * self.n
            vis = [False] * self.n
            vis[self.find_base(start)] = True
            qi = 0
            endpoint = NIL

            while qi < len(queue) and endpoint == NIL:
                u = queue[qi]
                qi += 1

                for v in self.graph[u]:
                    bu = self.find_base(u)
                    bv = self.find_base(v)
                    if bu == bv or vis[bv]:
                        continue

                    if self.mate[v] == NIL and v != start:
                        pred[v] = u
                        endpoint = v
                        break

                    if self.label[bv] != ODD:
                        pred[v] = u
                        vis[bv] = True
                        mv = self.mate[v]
                        if mv != NIL and not vis[self.find_base(mv)]:
                            pred[mv] = v
                            vis[self.find_base(mv)] = True
                            queue.append(mv)

            if endpoint != NIL:
                path = []
                curr = endpoint
                while curr != NIL:
                    path.append(curr)
                    curr = pred[curr]
                path.reverse()

                i = 0
                while i + 1 < len(path):
                    self.mate[path[i]] = path[i + 1]
                    self.mate[path[i + 1]] = path[i]
                    i += 2

    def maximum_matching(self):
        while self.phase_1():
            self.phase_2()

        matching = []
        for u in range(self.n):
            if self.mate[u] != NIL and self.mate[u] > u:
                matching.append((u, self.mate[u]))
        matching.sort()
        return matching


def validate_matching(n, graph, matching):
    deg = [0] * n
    errors = 0

    for u, v in matching:
        pos = bisect_left(graph[u], v)
        if pos >= len(graph[u]) or graph[u][pos] != v:
            print(f"ERROR: Edge ({u}, {v}) not in graph!", file=sys.stderr)
            errors += 1
        deg[u] += 1
        deg[v] += 1

    for i in range(n):
        if deg[i] > 1:
            print(f"ERROR: Vertex {i} in {deg[i]} edges!", file=sys.stderr)
            errors += 1

    matched = sum(1 for d in deg if d > 0)

    print()
    print("=== Validation Report ===")
    print(f"Matching size: {len(matching)}")
    print(f"Matched vertices: {matched}")
    print("VALIDATION FAILED" if errors > 0 else "VALIDATION PASSED")
    print("=========================")
    print()


def load_graph(filename):
    with open(filename, 'r') as f:
        n, m = map(int, f.readline().split())
        edges = []
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                edges.append((int(parts[0]), int(parts[1])))
    return n, edges


def main():
    print("Gabow's Scaling Algorithm (Optimized) - Python Implementation")
    print("===============================================================")
    print()

    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <filename>")
        sys.exit(1)

    n, edges = load_graph(sys.argv[1])
    print(f"Graph: {n} vertices, {len(edges)} edges")

    t0 = time.time()
    gabow = GabowOptimized(n, edges)
    matching = gabow.maximum_matching()
    t1 = time.time()

    validate_matching(n, gabow.graph, matching)

    print(f"Matching size: {len(matching)}")
    print(f"Time: {int((t1 - t0) * 1000)} ms")


if __name__ == "__main__":
    main()

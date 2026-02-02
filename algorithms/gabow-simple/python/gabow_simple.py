"""
Gabow's Algorithm (Simple) - O(VE) Maximum Matching

Python implementation — fully deterministic, no hash containers.
"""

import sys
import time
from bisect import bisect_left

NIL = -1


class GabowSimple:
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
        self.base = list(range(n))
        self.parent = [NIL] * n
        self.blossom = [False] * n
        self.visited = [False] * n

    def find_base(self, v):
        if self.base[v] != v:
            self.base[v] = self.find_base(self.base[v])
        return self.base[v]

    def find_lca(self, u, v):
        path = [False] * self.n

        for _ in range(self.n):
            u = self.find_base(u)
            path[u] = True
            if self.mate[u] == NIL:
                break
            mu = self.mate[u]
            if self.parent[mu] == NIL:
                break
            u = self.parent[mu]

        for _ in range(self.n):
            v = self.find_base(v)
            if path[v]:
                return v
            if self.mate[v] == NIL:
                break
            mv = self.mate[v]
            if self.parent[mv] == NIL:
                break
            v = self.parent[mv]

        return NIL

    def mark_blossom(self, u, lca, queue):
        for _ in range(self.n):
            if self.find_base(u) == lca:
                break
            bv = self.find_base(u)
            mu = self.mate[u]
            bw = self.find_base(mu)

            self.blossom[bv] = True
            self.blossom[bw] = True

            if not self.visited[bw]:
                self.visited[bw] = True
                queue.append(bw)

            if self.parent[mu] == NIL:
                break
            u = self.parent[mu]

    def contract_blossom(self, u, v, queue):
        lca = self.find_lca(u, v)
        if lca == NIL:
            return

        self.blossom = [False] * self.n
        self.mark_blossom(u, lca, queue)
        self.mark_blossom(v, lca, queue)

        for i in range(self.n):
            bi = self.find_base(i)
            if self.blossom[bi]:
                self.base[i] = lca
                if not self.visited[i]:
                    self.visited[i] = True
                    queue.append(i)

    def find_augmenting_path(self, start):
        for i in range(self.n):
            self.base[i] = i
            self.parent[i] = NIL
        self.visited = [False] * self.n

        queue = [start]
        self.visited[start] = True
        qi = 0

        while qi < len(queue):
            u = queue[qi]
            qi += 1

            for v in self.graph[u]:
                bu = self.find_base(u)
                bv = self.find_base(v)
                if bu == bv:
                    continue

                if self.mate[v] == NIL:
                    self.parent[v] = u
                    return True

                if not self.visited[bv]:
                    self.parent[v] = u
                    self.visited[bv] = True
                    w = self.mate[v]
                    bw = self.find_base(w)
                    self.visited[bw] = True
                    queue.append(w)
                else:
                    # check if same tree → blossom
                    ru = bu
                    for _ in range(self.n):
                        if self.mate[ru] == NIL:
                            break
                        mru = self.mate[ru]
                        if self.parent[mru] == NIL:
                            break
                        ru = self.find_base(self.parent[mru])

                    rv = bv
                    for _ in range(self.n):
                        if self.mate[rv] == NIL:
                            break
                        mrv = self.mate[rv]
                        if self.parent[mrv] == NIL:
                            break
                        rv = self.find_base(self.parent[mrv])

                    if ru == rv:
                        self.contract_blossom(u, v, queue)
        return False

    def augment_path(self, v):
        while self.parent[v] != NIL:
            pv = self.parent[v]
            ppv = self.mate[pv]
            self.mate[v] = pv
            self.mate[pv] = v
            if ppv == NIL:
                break
            v = ppv

    def maximum_matching(self):
        found = True
        while found:
            found = False
            for v in range(self.n):
                if self.mate[v] == NIL:
                    if self.find_augmenting_path(v):
                        for u in range(self.n):
                            if self.mate[u] == NIL and self.parent[u] != NIL:
                                self.augment_path(u)
                                found = True
                                break

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
    print("Gabow's Algorithm (Simple) - Python Implementation")
    print("====================================================")
    print()

    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <filename>")
        sys.exit(1)

    n, edges = load_graph(sys.argv[1])
    print(f"Graph: {n} vertices, {len(edges)} edges")

    t0 = time.time()
    gabow = GabowSimple(n, edges)
    matching = gabow.maximum_matching()
    t1 = time.time()

    validate_matching(n, gabow.graph, matching)

    print(f"Matching size: {len(matching)}")
    print(f"Time: {int((t1 - t0) * 1000)} ms")


if __name__ == "__main__":
    main()

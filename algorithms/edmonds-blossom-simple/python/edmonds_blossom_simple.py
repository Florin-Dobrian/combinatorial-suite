"""
Edmonds' Blossom Algorithm (Simple) - O(V^4) Maximum Matching

Python implementation — fully deterministic, no hash containers.
"""

import sys
import time
from bisect import bisect_left

NIL = -1


class EdmondsBlossomSimple:
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

    def find_base(self, v, base):
        c = v
        while base[c] != c:
            c = base[c]
        return c

    def find_blossom_base(self, v, w, parent, base):
        on_path = [False] * self.n
        c = v
        while c != NIL:
            on_path[self.find_base(c, base)] = True
            c = parent[c]
        c = w
        while c != NIL:
            bc = self.find_base(c, base)
            if on_path[bc]:
                return bc
            c = parent[c]
        return self.find_base(v, base)

    def trace_and_update(self, start, blossom_base, base, label, parent, queue, in_queue):
        c = start
        while True:
            cb = self.find_base(c, base)
            if cb == blossom_base:
                break
            base[c] = blossom_base
            if label[c] == 2:
                label[c] = 1
                if not in_queue[c]:
                    queue.append(c)
                    in_queue[c] = True
            if self.mate[c] == NIL:
                break
            mc = self.mate[c]
            base[mc] = blossom_base
            if parent[mc] == NIL:
                break
            c = parent[mc]

    def find_augmenting_path(self, start):
        parent = [NIL] * self.n
        base = list(range(self.n))
        label = [0] * self.n  # 0=unlabeled, 1=outer, 2=inner
        in_queue = [False] * self.n

        label[start] = 1
        queue = [start]
        in_queue[start] = True
        qi = 0

        while qi < len(queue):
            v = queue[qi]
            qi += 1
            vb = self.find_base(v, base)

            for w in self.graph[v]:
                wb = self.find_base(w, base)
                if vb == wb:
                    continue

                if label[w] == 0:
                    if self.mate[w] != NIL:
                        mw = self.mate[w]
                        label[w] = 2
                        label[mw] = 1
                        parent[w] = v
                        parent[mw] = w
                        if not in_queue[mw]:
                            queue.append(mw)
                            in_queue[mw] = True
                    else:
                        parent[w] = v
                        path = [w]
                        c = v
                        while c != NIL:
                            path.append(c)
                            c = parent[c]
                        return path
                elif label[w] == 1:
                    bb = self.find_blossom_base(v, w, parent, base)
                    self.trace_and_update(v, bb, base, label, parent, queue, in_queue)
                    self.trace_and_update(w, bb, base, label, parent, queue, in_queue)
        return None

    def augment(self, path):
        i = 0
        while i + 1 < len(path):
            u, v = path[i], path[i + 1]
            self.mate[u] = v
            self.mate[v] = u
            i += 2

    def maximum_matching(self):
        improved = True
        while improved:
            improved = False
            for v in range(self.n):
                if self.mate[v] == NIL:
                    path = self.find_augmenting_path(v)
                    if path is not None:
                        self.augment(path)
                        improved = True
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
    print("Edmonds' Blossom Algorithm (Simple) - Python Implementation")
    print("=============================================================")
    print()

    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <filename>")
        sys.exit(1)

    n, edges = load_graph(sys.argv[1])
    print(f"Graph: {n} vertices, {len(edges)} edges")

    t0 = time.time()
    eb = EdmondsBlossomSimple(n, edges)
    matching = eb.maximum_matching()
    t1 = time.time()

    validate_matching(n, eb.graph, matching)

    print(f"Matching size: {len(matching)}")
    print(f"Time: {int((t1 - t0) * 1000)} ms")


if __name__ == "__main__":
    main()

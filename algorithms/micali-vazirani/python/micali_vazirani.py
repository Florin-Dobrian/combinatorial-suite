"""
Micali-Vazirani Algorithm (Hybrid) - O(E√V) Maximum Matching

Hybrid approach:
- MV-style MIN phase (level building with even/odd tracking)
- Gabow-style MAX phase (simple path finding)

Python implementation — fully deterministic, no hash containers.
"""

import sys
import time
from bisect import bisect_left

NIL = -1
UNSET = float('inf')


class Node:
    __slots__ = ('preds', 'match', 'min_level', 'even_level', 'odd_level')

    def __init__(self):
        self.preds = []
        self.match = NIL
        self.min_level = UNSET
        self.even_level = UNSET
        self.odd_level = UNSET

    def set_min_level(self, level):
        self.min_level = level
        if level % 2:
            self.odd_level = level
        else:
            self.even_level = level

    def reset(self):
        self.preds.clear()
        self.min_level = UNSET
        self.even_level = UNSET
        self.odd_level = UNSET


class MicaliVazirani:
    def __init__(self, n, edges):
        self.n = n
        self.graph = [[] for _ in range(n)]
        self.nodes = [Node() for _ in range(n)]
        self.base = list(range(n))
        self.levels = []

        for u, v in edges:
            if u < n and v < n and u != v:
                self.graph[u].append(v)
                self.graph[v].append(u)
        for adj in self.graph:
            adj.sort()

    def find_base(self, v):
        if self.base[v] != v:
            self.base[v] = self.find_base(self.base[v])
        return self.base[v]

    def add_to_level(self, level, node):
        while len(self.levels) <= level:
            self.levels.append([])
        self.levels[level].append(node)

    def step_to(self, to, frm, level):
        next_level = level + 1
        tl = self.nodes[to].min_level
        if tl >= next_level:
            if tl != next_level:
                self.add_to_level(next_level, to)
                self.nodes[to].set_min_level(next_level)
            self.nodes[to].preds.append(frm)

    def phase_1(self):
        self.levels.clear()
        for i in range(self.n):
            self.base[i] = i
            self.nodes[i].reset()

        for i in range(self.n):
            if self.nodes[i].match == NIL:
                self.add_to_level(0, i)
                self.nodes[i].set_min_level(0)

        for i in range(self.n):
            if i >= len(self.levels) or not self.levels[i]:
                continue
            level_snap = list(self.levels[i])
            for cur in level_snap:
                if i % 2 == 0:
                    for nb in self.graph[cur]:
                        if nb != self.nodes[cur].match:
                            self.step_to(nb, cur, i)
                else:
                    if self.nodes[cur].match != NIL:
                        self.step_to(self.nodes[cur].match, cur, i)

    def phase_2(self):
        found = False

        for start in range(self.n):
            if self.nodes[start].match != NIL or self.nodes[start].min_level != 0:
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

                    if self.nodes[v].match == NIL and v != start:
                        pred[v] = u
                        endpoint = v
                        break

                    pred[v] = u
                    vis[bv] = True
                    mv = self.nodes[v].match
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
                    self.nodes[path[i]].match = path[i + 1]
                    self.nodes[path[i + 1]].match = path[i]
                    i += 2
                found = True

        return found

    def maximum_matching(self):
        while True:
            self.phase_1()
            if not self.phase_2():
                break

        matching = []
        for u in range(self.n):
            if self.nodes[u].match != NIL and self.nodes[u].match > u:
                matching.append((u, self.nodes[u].match))
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
    print("Micali-Vazirani Algorithm (Hybrid) - Python Implementation")
    print("=============================================================")
    print()

    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <filename>")
        sys.exit(1)

    n, edges = load_graph(sys.argv[1])
    print(f"Graph: {n} vertices, {len(edges)} edges")

    t0 = time.time()
    mv = MicaliVazirani(n, edges)
    matching = mv.maximum_matching()
    t1 = time.time()

    validate_matching(n, mv.graph, matching)

    print(f"Matching size: {len(matching)}")
    print(f"Time: {int((t1 - t0) * 1000)} ms")


if __name__ == "__main__":
    main()

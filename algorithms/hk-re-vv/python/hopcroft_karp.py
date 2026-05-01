"""
Hopcroft-Karp Algorithm - O(EâˆšV) Maximum Bipartite Matching

Python implementation â€” fully deterministic, no hash containers.
"""

import sys
import time
from bisect import bisect_left

NIL = -1
INF = float('inf')


class HopcroftKarp:
    def __init__(self, left_count, right_count, edges):
        self.left_count = left_count
        self.right_count = right_count
        self.graph = [[] for _ in range(left_count)]
        for u, v in edges:
            if 0 <= u < left_count and 0 <= v < right_count:
                self.graph[u].append(v)
        for adj in self.graph:
            adj.sort()

        self.pair_left = [NIL] * left_count
        self.pair_right = [NIL] * right_count
        self.dist = [0] * (left_count + 1)

    def bfs(self):
        queue = []
        qi = 0

        for u in range(self.left_count):
            if self.pair_left[u] == NIL:
                self.dist[u] = 0
                queue.append(u)
            else:
                self.dist[u] = INF

        self.dist[self.left_count] = INF

        while qi < len(queue):
            u = queue[qi]
            qi += 1
            if self.dist[u] < self.dist[self.left_count]:
                for v in self.graph[u]:
                    paired = self.left_count if self.pair_right[v] == NIL else self.pair_right[v]
                    if self.dist[paired] == INF:
                        self.dist[paired] = self.dist[u] + 1
                        if self.pair_right[v] != NIL:
                            queue.append(self.pair_right[v])

        return self.dist[self.left_count] != INF

    def dfs(self, u):
        if u == NIL:
            return True

        for v in self.graph[u]:
            paired = self.left_count if self.pair_right[v] == NIL else self.pair_right[v]
            if self.dist[paired] == self.dist[u] + 1:
                if self.dfs(self.pair_right[v]):
                    self.pair_right[v] = u
                    self.pair_left[u] = v
                    return True

        self.dist[u] = INF
        return False

    def maximum_matching(self, greedy_mode=0):
        if greedy_mode == 1:
            self.greedy_size = self._greedy_init()
        elif greedy_mode == 2:
            self.greedy_size = self._greedy_init_md()

        while self.bfs():
            for u in range(self.left_count):
                if self.pair_left[u] == NIL:
                    self.dfs(u)

        matching = []
        for u in range(self.left_count):
            if self.pair_left[u] != NIL:
                matching.append((u, self.pair_left[u]))
        matching.sort()
        return matching

    greedy_size = 0

    def _greedy_init(self):
        cnt = 0
        for u in range(self.left_count):
            if self.pair_left[u] != NIL:
                continue
            for v in self.graph[u]:
                if self.pair_right[v] == NIL:
                    self.pair_left[u] = v
                    self.pair_right[v] = u
                    cnt += 1
                    break
        return cnt

    def _greedy_init_md(self):
        cnt = 0
        lc = self.left_count
        rc = self.right_count
        # Compute right-side degrees
        rdeg = [0] * rc
        for u in range(lc):
            for v in self.graph[u]:
                rdeg[v] += 1
        # Compute left-side degrees
        ldeg = [len(self.graph[u]) for u in range(lc)]
        order = sorted(range(lc), key=lambda x: (ldeg[x], x))
        for u in order:
            if self.pair_left[u] != NIL:
                continue
            best = NIL
            best_deg = float('inf')
            for v in self.graph[u]:
                if self.pair_right[v] == NIL and rdeg[v] < best_deg:
                    best = v
                    best_deg = rdeg[v]
            if best >= 0:
                self.pair_left[u] = best
                self.pair_right[best] = u
                cnt += 1
        return cnt


def validate_matching(left_count, right_count, graph, matching):
    left_deg = [0] * left_count
    right_deg = [0] * right_count
    errors = 0

    for u, v in matching:
        pos = bisect_left(graph[u], v)
        if pos >= len(graph[u]) or graph[u][pos] != v:
            print(f"ERROR: Edge ({u}, {v}) not in graph!", file=sys.stderr)
            errors += 1
        left_deg[u] += 1
        right_deg[v] += 1

    for i in range(left_count):
        if left_deg[i] > 1:
            print(f"ERROR: Left vertex {i} in {left_deg[i]} edges!", file=sys.stderr)
            errors += 1
    for i in range(right_count):
        if right_deg[i] > 1:
            print(f"ERROR: Right vertex {i} in {right_deg[i]} edges!", file=sys.stderr)
            errors += 1

    matched_left = sum(1 for d in left_deg if d > 0)
    matched_right = sum(1 for d in right_deg if d > 0)

    print()
    print("=== Validation Report ===")
    print(f"Matching size: {len(matching)}")
    print(f"Matched vertices: {matched_left} left, {matched_right} right")
    print("VALIDATION FAILED" if errors > 0 else "VALIDATION PASSED")
    print("=========================")
    print()


def load_graph(filename):
    with open(filename, 'r') as f:
        parts = f.readline().split()
        left_count, right_count, m = int(parts[0]), int(parts[1]), int(parts[2])
        edges = []
        for line in f:
            p = line.split()
            if len(p) >= 2:
                edges.append((int(p[0]), int(p[1])))
    return left_count, right_count, edges


def main():
    print("Hopcroft-Karp Algorithm - Python Implementation")
    print("==================================================")
    print()

    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <filename> [--greedy|--greedy-md]")
        sys.exit(1)

    greedy_mode = 0
    for arg in sys.argv[2:]:
        if arg == "--greedy":
            greedy_mode = 1
        elif arg == "--greedy-md":
            greedy_mode = 2

    left_count, right_count, edges = load_graph(sys.argv[1])
    print(f"Graph: {left_count} left, {right_count} right, {len(edges)} edges")

    t0 = time.time()
    hk = HopcroftKarp(left_count, right_count, edges)
    matching = hk.maximum_matching(greedy_mode)
    t1 = time.time()

    validate_matching(left_count, right_count, hk.graph, matching)

    print(f"Matching size: {len(matching)}")
    if greedy_mode > 0:
        print(f"Greedy init size: {hk.greedy_size}")
        if matching:
            print(f"Greedy/Final: {100.0 * hk.greedy_size / len(matching):.2f}%")
        else:
            print("Greedy/Final: NA")
    print(f"Time: {int((t1 - t0) * 1000)} ms")


if __name__ == "__main__":
    main()

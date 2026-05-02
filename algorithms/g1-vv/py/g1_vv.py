"""
Gabow's Algorithm (Simple) – O(V * E) Maximum Matching

Faithful to Gabow 1976: forest BFS with blossom contraction via
union-find. No physical contraction – bases are tracked virtually.
Epoch-based interleaved LCA, path-only contraction, bridge recording
for augmentation through blossoms.

Forest search: each iteration labels ALL free vertices as EVEN roots
simultaneously and grows a search forest. An augmenting path is found
when two different trees meet (EVEN-EVEN edge across trees, detected
by find_lca returning NIL). One augmentation per iteration, then full
reset and repeat until no augmenting path exists.

Complexity: O(V * E) – each iteration does O(E) work, at most V/2
augmentations total.

All integers, no hash containers, fully deterministic.

Python implementation – faithful translation of the C++ version.
"""

import sys
import time

NIL = -1
UNLABELED = 0
EVEN = 1
ODD = 2


class GabowSimple:
    def __init__(self, n, edges):
        self.n = n
        self.graph = [[] for _ in range(n)]
        for u, v in edges:
            if 0 <= u < n and 0 <= v < n and u != v:
                self.graph[u].append(v)
                self.graph[v].append(u)
        for i in range(n):
            self.graph[i] = sorted(set(self.graph[i]))

        self.mate = [NIL] * n
        self.base = list(range(n))
        self.parent = [NIL] * n
        self.label = [UNLABELED] * n
        self.bridge_src = [NIL] * n
        self.bridge_tgt = [NIL] * n
        self.lca_tag1 = [0] * n
        self.lca_tag2 = [0] * n
        self.lca_epoch = 0
        self.greedy_size = 0

    # ---- Greedy initialization ----

    def greedy_init(self):
        cnt = 0
        for u in range(self.n):
            if self.mate[u] != NIL:
                continue
            for v in self.graph[u]:
                if self.mate[v] == NIL:
                    self.mate[u] = v
                    self.mate[v] = u
                    cnt += 1
                    break
        return cnt

    def greedy_init_md(self):
        cnt = 0
        n = self.n
        deg = [len(self.graph[u]) for u in range(n)]
        order = sorted(range(n), key=lambda x: (deg[x], x))
        for u in order:
            if self.mate[u] != NIL:
                continue
            best = NIL
            best_deg = float('inf')
            for v in self.graph[u]:
                if self.mate[v] == NIL and deg[v] < best_deg:
                    best = v
                    best_deg = deg[v]
            if best >= 0:
                self.mate[u] = best
                self.mate[best] = u
                cnt += 1
        return cnt

    # ---- Union-find base with path halving ----

    def find_base(self, v):
        while self.base[v] != v:
            self.base[v] = self.base[self.base[v]]
            v = self.base[v]
        return v

    # ---- Interleaved LCA using epoch tags ----

    def find_lca(self, u, v):
        self.lca_epoch += 1
        ep = self.lca_epoch
        hx = self.find_base(u)
        hy = self.find_base(v)
        self.lca_tag1[hx] = ep
        self.lca_tag2[hy] = ep
        while True:
            if self.lca_tag1[hy] == ep:
                return hy
            if self.lca_tag2[hx] == ep:
                return hx
            hxr = (self.mate[hx] == NIL)
            hyr = (self.mate[hy] == NIL)
            if hxr and hyr:
                return NIL  # different trees
            if not hxr:
                hx = self.find_base(self.parent[self.mate[hx]])
                self.lca_tag1[hx] = ep
            if not hyr:
                hy = self.find_base(self.parent[self.mate[hy]])
                self.lca_tag2[hy] = ep

    # ---- Path-only contraction ----

    def shrink_path(self, lca, x, y, queue, qtail):
        v = self.find_base(x)
        while v != lca:
            mv = self.mate[v]
            self.base[self.find_base(v)] = lca
            self.base[self.find_base(mv)] = lca
            self.base[lca] = lca
            self.bridge_src[mv] = x
            self.bridge_tgt[mv] = y
            if self.label[mv] != EVEN:
                self.label[mv] = EVEN
                queue[qtail[0]] = mv
                qtail[0] += 1
            v = self.find_base(self.parent[mv])
        return qtail

    # ---- Trace path for augmentation ----

    def trace_path(self, v, u, pairs):
        """Trace from vertex v to vertex u (or to root if u==NIL),
        collecting edge pairs for augmentation."""
        # Iterative with explicit stack
        # Frame: (v, u, phase, sb, tb)
        stack = [(v, u, 0, 0, 0)]

        while stack:
            fv, fu, phase, sb, tb = stack[-1]

            if fv == fu:
                stack.pop()
                continue

            if phase == 0:
                if self.bridge_src[fv] == NIL:
                    # Originally EVEN vertex (no bridge)
                    if self.mate[fv] == NIL:
                        # Root (free vertex) – done
                        stack.pop()
                        continue
                    mv = self.mate[fv]
                    pmv = self.parent[mv]
                    pairs.append((mv, pmv))
                    stack[-1] = (pmv, fu, 0, 0, 0)
                    continue
                # Has bridge – originally ODD, absorbed into blossom
                sb = self.bridge_src[fv]
                tb = self.bridge_tgt[fv]
                stack[-1] = (fv, fu, 1, sb, tb)
                stack.append((sb, self.mate[fv], 0, 0, 0))
                continue

            if phase == 1:
                pairs.append((sb, tb))
                stack[-1] = (fv, fu, 2, sb, tb)
                stack.append((tb, fu, 0, 0, 0))
                continue

            # phase == 2
            stack.pop()

    # ---- Augment along two-sided path ----

    def augment_two_sides(self, u, v):
        pairs = [(u, v)]
        self.trace_path(u, NIL, pairs)
        self.trace_path(v, NIL, pairs)
        for a, b in pairs:
            self.mate[a] = b
            self.mate[b] = a

    # ---- Main search iteration ----

    def find_and_augment(self):
        n = self.n
        # Reset per-iteration state
        for i in range(n):
            self.base[i] = i
            self.parent[i] = NIL
            self.label[i] = UNLABELED
            self.bridge_src[i] = NIL
            self.bridge_tgt[i] = NIL

        queue = [0] * n
        qhead = 0
        qtail = [0]  # use list for mutability in shrink_path

        # All free vertices become EVEN roots
        for v in range(n):
            if self.mate[v] == NIL:
                self.label[v] = EVEN
                queue[qtail[0]] = v
                qtail[0] += 1

        while qhead < qtail[0]:
            u = queue[qhead]
            qhead += 1
            # Check that u is still effectively EVEN
            if self.label[self.find_base(u)] != EVEN:
                continue

            for v in self.graph[u]:
                bu = self.find_base(u)
                bv = self.find_base(v)
                if bu == bv:
                    continue  # same blossom
                if v == self.mate[u]:
                    continue  # skip matching edge

                if self.label[bv] == UNLABELED:
                    # v is matched and unlabeled -> grow step
                    self.label[v] = ODD
                    self.parent[v] = u
                    w = self.mate[v]
                    self.label[w] = EVEN
                    queue[qtail[0]] = w
                    qtail[0] += 1

                elif self.label[bv] == EVEN:
                    # EVEN-EVEN edge: blossom or augmenting path
                    lca = self.find_lca(u, v)
                    if lca != NIL:
                        # Same tree -> blossom contraction
                        self.shrink_path(lca, u, v, queue, qtail)
                        self.shrink_path(lca, v, u, queue, qtail)
                    else:
                        # Different trees -> augmenting path!
                        self.augment_two_sides(u, v)
                        return True
                # label[bv] == ODD: ignore

        return False

    def maximum_matching(self, greedy_mode=0):
        if greedy_mode == 1:
            self.greedy_size = self.greedy_init()
        elif greedy_mode == 2:
            self.greedy_size = self.greedy_init_md()

        while self.find_and_augment():
            pass

        matching = []
        for u in range(self.n):
            if self.mate[u] != NIL and self.mate[u] > u:
                matching.append((u, self.mate[u]))
        matching.sort()
        return matching


# ---- Validation ----

def validate_matching(n, graph, matching):
    deg = [0] * n
    errors = 0
    for u, v in matching:
        if v not in graph[u]:
            print(f"ERROR: Edge ({u},{v}) not in graph!", file=sys.stderr)
            errors += 1
        deg[u] += 1
        deg[v] += 1
    for i in range(n):
        if deg[i] > 1:
            print(f"ERROR: Vertex {i} in {deg[i]} edges!", file=sys.stderr)
            errors += 1
    mt = sum(1 for d in deg if d > 0)
    status = "VALIDATION FAILED" if errors > 0 else "VALIDATION PASSED"
    print(f"\n=== Validation Report ===\n"
          f"Matching size: {len(matching)}\n"
          f"Matched vertices: {mt}\n"
          f"{status}\n"
          f"=========================\n")


# ---- Graph loading ----

def load_graph(filename):
    with open(filename) as f:
        header = f.readline().split()
        n, m = int(header[0]), int(header[1])
        edges = []
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                edges.append((int(parts[0]), int(parts[1])))
    return n, edges


# ---- Main ----

def main():
    print("Gabow's Algorithm (Simple) - Python Implementation")
    print("===================================================")
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

    n, edges = load_graph(sys.argv[1])
    print(f"Graph: {n} vertices, {len(edges)} edges")

    t0 = time.time()
    gabow = GabowSimple(n, edges)
    matching = gabow.maximum_matching(greedy_mode)
    t1 = time.time()

    validate_matching(n, gabow.graph, matching)

    print(f"Matching size: {len(matching)}")
    if greedy_mode > 0:
        gs = gabow.greedy_size
        fs = len(matching)
        print(f"Greedy init size: {gs}")
        if fs > 0:
            print(f"Greedy/Final: {100.0 * gs / fs:.2f}%")
        else:
            print("Greedy/Final: NA")
    print(f"Time: {int((t1 - t0) * 1000)} ms")


if __name__ == "__main__":
    main()

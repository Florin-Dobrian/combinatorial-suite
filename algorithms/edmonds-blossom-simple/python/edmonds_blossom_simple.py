"""
Edmonds' Blossom Algorithm (Simple) – Unweighted Maximum Cardinality Matching

Single-source BFS (tree, not forest). Each iteration grows an alternating
tree from one free vertex. Blossoms are shrunk into supernodes during the
search and expanded back to regular vertices after each search completes.

Blossom IDs are reset to n at the start of each BFS, so all indices fit
in plain Python ints.

Blossom data structure follows NetworkX: each blossom stores childs
(sub-blossom IDs in cycle order) and edges (connecting edge pairs).
augment_blossom recurses into nested sub-blossoms for correct path lifting.

Complexity: O(V^2 * E) worst case.

Python implementation – faithful translation of the C++ version.
"""

import sys
import time

NIL = -1


class Blos:
    __slots__ = ('childs', 'edges')
    def __init__(self):
        self.childs = []
        self.edges = []


class Solver:
    def __init__(self, n, edges):
        self.n = n
        self.adj = [[] for _ in range(n)]
        for u, v in edges:
            if u != v and 0 <= u < n and 0 <= v < n:
                self.adj[u].append(v)
                self.adj[v].append(u)
        for i in range(n):
            self.adj[i] = sorted(set(self.adj[i]))

        self.mate = [NIL] * n

        # Blossom storage. IDs 0..n-1 are trivial (one vertex each).
        # Non-trivial blossoms have id in [n, nblos).
        self.blos = [None] * n  # slots 0..n-1 unused; extended as needed
        self.nblos = n

        self.inblossom = list(range(n))
        self.blossomparent = [NIL] * n
        self.blossombase = list(range(n))

        # Per-search state (resized in ensure / reset_blossoms)
        self.label = []
        self.labeledge = []
        self.queue = []

        self.greedy_size = 0

    def ensure(self, b):
        if b < len(self.label):
            return
        old = len(self.label)
        self.label.extend([0] * (b + 1 - old))
        self.labeledge.extend([(NIL, NIL)] * (b + 1 - old))
        while len(self.blossomparent) <= b:
            self.blossomparent.append(NIL)
        while len(self.blossombase) <= b:
            self.blossombase.append(NIL)

    def is_blossom(self, b):
        return b >= self.n

    def leaves(self, b):
        """Yield all vertices (leaf nodes) inside blossom b."""
        if not self.is_blossom(b):
            return [b]
        result = []
        stack = [b]
        while stack:
            x = stack.pop()
            if not self.is_blossom(x):
                result.append(x)
            else:
                for c in self.blos[x].childs:
                    stack.append(c)
        return result

    # ---- Reset for a new BFS ----

    def reset_blossoms(self):
        n = self.n
        self.nblos = n
        self.blos = self.blos[:n]
        for i in range(n):
            self.inblossom[i] = i
            self.blossombase[i] = i
            self.blossomparent[i] = NIL
        self.label = [0] * n
        self.labeledge = [(NIL, NIL)] * n
        self.queue = []

    # ---- Tree building ----

    def assign_label(self, w, t, v):
        b = self.inblossom[w]
        self.ensure(b)
        self.label[b] = t
        self.label[w] = t
        if v != NIL:
            self.labeledge[w] = self.labeledge[b] = (v, w)
        else:
            self.labeledge[w] = self.labeledge[b] = (NIL, NIL)
        if t == 1:
            # S-blossom: add its leaves to the BFS queue
            for u in self.leaves(b):
                self.queue.append(u)
        elif t == 2:
            # T-blossom: label the mate of its base as S
            base = self.blossombase[b]
            self.assign_label(self.mate[base], 1, base)

    # ---- Blossom detection ----

    def scan_blossom(self, v, w):
        """Trace from two S-vertices to find their LCA (blossom base).
        Returns base vertex, or -2 if different trees."""
        path = []
        base = -2
        while v != -2 or w != -2:
            if v != -2:
                b = self.inblossom[v]
                if self.label[b] & 4:
                    base = self.blossombase[b]
                    break
                path.append(b)
                self.label[b] = 5  # breadcrumb
                le = self.labeledge[b]
                if le[0] == NIL:
                    v = -2  # reached root
                else:
                    v = le[0]
                    bt = self.inblossom[v]
                    v = self.labeledge[bt][0]
                if w != -2:
                    v, w = w, v
            else:
                v, w = w, v
        for b in path:
            self.label[b] = 1  # restore breadcrumbs
        return base

    # ---- Blossom contraction ----

    def add_blossom(self, base, v, w):
        bb = self.inblossom[base]
        bv = self.inblossom[v]
        bw = self.inblossom[w]

        bid = self.nblos
        self.nblos += 1
        if bid >= len(self.blos):
            self.blos.append(Blos())
        else:
            self.blos[bid] = Blos()
        self.ensure(bid)
        self.blossombase[bid] = base
        self.blossomparent[bid] = NIL
        self.blossomparent[bb] = bid

        bl = self.blos[bid]
        bl.edges.append((v, w))  # bridge edge

        # Trace from v back to base
        cv = v
        bcv = bv
        while bcv != bb:
            self.blossomparent[bcv] = bid
            bl.childs.append(bcv)
            bl.edges.append(self.labeledge[bcv])
            cv = self.labeledge[bcv][0]
            bcv = self.inblossom[cv]
        bl.childs.append(bb)
        bl.childs.reverse()
        bl.edges.reverse()

        # Trace from w back to base
        cw = w
        bcw = bw
        while bcw != bb:
            self.blossomparent[bcw] = bid
            bl.childs.append(bcw)
            le = self.labeledge[bcw]
            bl.edges.append((le[1], le[0]))  # reversed
            cw = self.labeledge[bcw][0]
            bcw = self.inblossom[cw]

        self.label[bid] = 1
        self.labeledge[bid] = self.labeledge[bb]

        # Relabel: T-vertices inside the blossom become S
        for u in self.leaves(bid):
            if self.label[self.inblossom[u]] == 2:
                self.queue.append(u)
            self.inblossom[u] = bid

    # ---- Blossom expansion ----

    def expand_blossom(self, b, endstage):
        stack = [(b, endstage, 0)]

        while stack:
            fb, fend, fidx = stack[-1]
            bl = self.blos[fb]

            if fidx < len(bl.childs):
                s = bl.childs[fidx]
                stack[-1] = (fb, fend, fidx + 1)
                self.blossomparent[s] = NIL
                if self.is_blossom(s):
                    if fend:
                        stack.append((s, True, 0))
                        continue
                    else:
                        for u in self.leaves(s):
                            self.inblossom[u] = s
                else:
                    self.inblossom[s] = s
            else:
                # All children processed
                if not fend and self.label[fb] == 2:
                    # Mid-stage T-blossom expansion: relabel children
                    bl2 = self.blos[fb]
                    entrychild = self.inblossom[self.labeledge[fb][1]]
                    k = len(bl2.childs)
                    j = 0
                    while j < k:
                        if bl2.childs[j] == entrychild:
                            break
                        j += 1
                    if j & 1:
                        j -= k
                        jstep = 1
                    else:
                        jstep = -1
                    lv_ = self.labeledge[fb][0]
                    lw_ = self.labeledge[fb][1]
                    while j != 0:
                        if jstep == 1:
                            pp = bl2.edges[j % k][0]
                            qq = bl2.edges[j % k][1]
                        else:
                            ei = (j - 1) % k
                            qq = bl2.edges[ei][0]
                            pp = bl2.edges[ei][1]
                        self.label[lw_] = 0
                        self.label[qq] = 0
                        self.assign_label(lw_, 2, lv_)
                        j += jstep
                        if jstep == 1:
                            lv_ = bl2.edges[j % k][0]
                            lw_ = bl2.edges[j % k][1]
                        else:
                            ei = (j - 1) % k
                            lw_ = bl2.edges[ei][0]
                            lv_ = bl2.edges[ei][1]
                        j += jstep

                    bwi = bl2.childs[j % k]
                    self.ensure(bwi)
                    self.label[lw_] = self.label[bwi] = 2
                    self.labeledge[lw_] = self.labeledge[bwi] = (lv_, lw_)
                    j += jstep
                    while bl2.childs[j % k] != entrychild:
                        bvi = bl2.childs[j % k]
                        self.ensure(bvi)
                        if self.label[bvi] == 1:
                            j += jstep
                            continue
                        found_v = NIL
                        if self.is_blossom(bvi):
                            for u in self.leaves(bvi):
                                if self.label[u]:
                                    found_v = u
                                    break
                        else:
                            found_v = bvi
                        if found_v != NIL and self.label[found_v]:
                            self.label[found_v] = 0
                            self.label[self.mate[self.blossombase[bvi]]] = 0
                            self.assign_label(found_v, 2, self.labeledge[found_v][0])
                        j += jstep

                self.label[fb] = 0
                bl.childs.clear()
                bl.edges.clear()
                stack.pop()

    # ---- Augmentation through blossoms ----

    def augment_blossom(self, b, v):
        # Iterative version using explicit stack
        # Frame: (b, v, phase, i, j, jstep)
        stack = [(b, v, 0, 0, 0, 0)]

        while stack:
            f = stack[-1]
            fb, fv, fphase, fi, fj, fjstep = f

            if fphase == 0:
                # Find sub-blossom containing v
                t = fv
                while self.blossomparent[t] != fb:
                    t = self.blossomparent[t]
                bl = self.blos[fb]
                k = len(bl.childs)
                i = 0
                while i < k:
                    if bl.childs[i] == t:
                        break
                    i += 1
                if self.is_blossom(t):
                    stack[-1] = (fb, fv, 1, i, 0, 0)
                    stack.append((t, fv, 0, 0, 0, 0))
                    continue
                if i & 1:
                    j = i - k
                    jstep = 1
                else:
                    j = i
                    jstep = -1
                stack[-1] = (fb, fv, 2, i, j, jstep)
                continue

            if fphase == 1:
                # After recursion into sub-blossom
                bl = self.blos[fb]
                k = len(bl.childs)
                if fi & 1:
                    j = fi - k
                    jstep = 1
                else:
                    j = fi
                    jstep = -1
                stack[-1] = (fb, fv, 2, fi, j, jstep)
                continue

            if fphase == 2:
                # Main loop: walk from position i toward position 0
                bl = self.blos[fb]
                k = len(bl.childs)
                if fj == 0:
                    # Done: rotate childs/edges so new base is first
                    if fi > 0:
                        bl.childs = bl.childs[fi:] + bl.childs[:fi]
                        bl.edges = bl.edges[fi:] + bl.edges[:fi]
                    self.blossombase[fb] = fv
                    stack.pop()
                    continue
                # Step to next pair of sub-blossoms
                fj += fjstep
                idx1 = fj % k
                c1 = bl.childs[idx1]
                if fjstep == 1:
                    ww = bl.edges[idx1][0]
                    xx = bl.edges[idx1][1]
                else:
                    ei = (fj - 1) % k
                    xx = bl.edges[ei][0]
                    ww = bl.edges[ei][1]
                if self.is_blossom(c1):
                    stack[-1] = (fb, fv, 3, fi, fj, fjstep)
                    stack.append((c1, ww, 0, 0, 0, 0))
                    continue
                stack[-1] = (fb, fv, 3, fi, fj, fjstep)

            if fphase == 3:
                # After optional recursion for c1, step to c2
                bl = self.blos[fb]
                k = len(bl.childs)
                idx1 = fj % k
                if fjstep == 1:
                    ww = bl.edges[idx1][0]
                    xx = bl.edges[idx1][1]
                else:
                    ei = (fj - 1) % k
                    xx = bl.edges[ei][0]
                    ww = bl.edges[ei][1]
                fj += fjstep
                idx2 = fj % k
                c2 = bl.childs[idx2]
                if self.is_blossom(c2):
                    stack[-1] = (fb, fv, 4, fi, fj, fjstep)
                    stack.append((c2, xx, 0, 0, 0, 0))
                    continue
                stack[-1] = (fb, fv, 4, fi, fj, fjstep)

            if fphase == 4:
                # After optional recursion for c2, set mate pair
                bl = self.blos[fb]
                k = len(bl.childs)
                prev_j = fj - fjstep
                idx1 = prev_j % k
                if fjstep == 1:
                    ww = bl.edges[idx1][0]
                    xx = bl.edges[idx1][1]
                else:
                    ei = (prev_j - 1) % k
                    xx = bl.edges[ei][0]
                    ww = bl.edges[ei][1]
                self.mate[ww] = xx
                self.mate[xx] = ww
                stack[-1] = (fb, fv, 2, fi, fj, fjstep)  # continue loop

    # ---- Augmenting path: trace from v back to root ----

    def augment_path(self, v, w):
        s, j = v, w
        while True:
            bs = self.inblossom[s]
            if self.is_blossom(bs):
                self.augment_blossom(bs, s)
            self.mate[s] = j
            le = self.labeledge[bs]
            if le[0] == NIL:
                break  # root
            t = le[0]  # T-vertex
            bt = self.inblossom[t]
            le2 = self.labeledge[bt]
            s = le2[0]
            j = le2[1]
            if self.is_blossom(bt):
                self.augment_blossom(bt, j)
            self.mate[j] = s
        self.mate[w] = v

    # ---- Greedy initialization ----

    def greedy_init(self):
        cnt = 0
        for u in range(self.n):
            if self.mate[u] != NIL:
                continue
            for v in self.adj[u]:
                if self.mate[v] == NIL:
                    self.mate[u] = v
                    self.mate[v] = u
                    cnt += 1
                    break
        return cnt

    def greedy_init_md(self):
        cnt = 0
        n = self.n
        deg = [0] * n
        for u in range(n):
            deg[u] = len(self.adj[u])
        order = sorted(range(n), key=lambda x: (deg[x], x))
        for u in order:
            if self.mate[u] != NIL:
                continue
            best = NIL
            bd = float('inf')
            for v in self.adj[u]:
                if self.mate[v] == NIL and deg[v] < bd:
                    best = v
                    bd = deg[v]
            if best >= 0:
                self.mate[u] = best
                self.mate[best] = u
                cnt += 1
        return cnt

    # ---- Main solver ----

    def solve(self, greedy_mode=0):
        if greedy_mode == 1:
            self.greedy_size = self.greedy_init()
        elif greedy_mode == 2:
            self.greedy_size = self.greedy_init_md()

        improved = True
        while improved:
            improved = False
            for root in range(self.n):
                if self.mate[root] != NIL:
                    continue

                # Fresh search from this root
                self.reset_blossoms()
                self.assign_label(root, 1, NIL)

                augmented = False
                while self.queue and not augmented:
                    v = self.queue.pop()
                    if self.label[self.inblossom[v]] != 1:
                        continue  # stale
                    for w in self.adj[v]:
                        bv = self.inblossom[v]
                        bw = self.inblossom[w]
                        if bv == bw:
                            continue
                        self.ensure(bw)
                        if self.label[bw] == 0:
                            if self.mate[w] == NIL:
                                self.augment_path(v, w)
                                augmented = True
                                break
                            self.assign_label(w, 2, v)
                        elif self.label[bw] == 1:
                            base = self.scan_blossom(v, w)
                            if base >= 0:
                                self.add_blossom(base, v, w)

                # Expand all remaining blossoms (endstage)
                for b in range(self.n, self.nblos):
                    if (b < len(self.blos) and self.blos[b] is not None
                            and self.blos[b].childs and self.blossomparent[b] == NIL):
                        self.expand_blossom(b, True)

                if augmented:
                    improved = True
                    break

        result = []
        for u in range(self.n):
            if self.mate[u] > u:
                result.append((u, self.mate[u]))
        result.sort()
        return result


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
    print("Edmonds' Blossom Algorithm (Simple) - Python Implementation")
    print("=============================================================")
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
    sol = Solver(n, edges)
    matching = sol.solve(greedy_mode)
    t1 = time.time()

    validate_matching(n, sol.adj, matching)

    print(f"Matching size: {len(matching)}")
    if greedy_mode > 0:
        print(f"Greedy init size: {sol.greedy_size}")
        if matching:
            print(f"Greedy/Final: {100.0 * sol.greedy_size / len(matching):.2f}%")
        else:
            print("Greedy/Final: NA")
    print(f"Time: {int((t1 - t0) * 1000)} ms")


if __name__ == "__main__":
    main()

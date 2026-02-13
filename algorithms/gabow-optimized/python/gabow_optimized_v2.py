"""
Gabow's Scaling Algorithm (Optimized) - O(E√V) Maximum Matching

Version 2: precomputed H-adjacency, no contracted_into.

Pure cardinality (unweighted) — integer weights conceptually all 1.
Phase 1: BFS by levels (Delta), detect blossoms. Build contracted graph H
         with precomputed adjacency lists.
Phase 2: Find all shortest augmenting paths in H (iterative DFS with blossom
         contraction) using H-adjacency. Unfold to G via bridges.

Based on LEDA-7's mc_matching_gabow architecture, adapted for pure cardinality.
All integers, fully deterministic.
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
            if 0 <= u < n and 0 <= v < n and u != v:
                self.graph[u].append(v)
                self.graph[v].append(u)
        for i in range(n):
            self.graph[i] = sorted(set(self.graph[i]))

        self.mate = [NIL] * n
        self.label = [UNLABELED] * n
        self.parent = [NIL] * n
        self.source_bridge = [NIL] * n
        self.target_bridge = [NIL] * n
        self.base_par = list(range(n))
        self.dbase_par = list(range(n))
        self.level_queue = [[] for _ in range(n + 2)]
        self.lca_tag1 = [0] * n
        self.lca_tag2 = [0] * n
        self.lca_epoch = 0
        self.in_tree = [False] * n
        self.tree_nodes = []
        self.delta = 0

        self.h_adj = [[] for _ in range(n)]  # precomputed H-adjacency

        self.rep = [0] * n
        self.mate_h = [NIL] * n
        self.label_h = [UNLABELED] * n
        self.parent_h_src = [NIL] * n
        self.parent_h_tgt = [NIL] * n
        self.bridge_h_src = [NIL] * n
        self.bridge_h_tgt = [NIL] * n
        self.dir_h = [0] * n
        self.even_time_h = [0] * n
        self.t_h = 0
        self.db2_par = list(range(n))

    # ---- union-find: base ----
    def find_base(self, v):
        while self.base_par[v] != v:
            self.base_par[v] = self.base_par[self.base_par[v]]
            v = self.base_par[v]
        return v

    def union_base(self, a, b, r):
        a = self.find_base(a)
        b = self.find_base(b)
        self.base_par[a] = r
        self.base_par[b] = r

    # ---- union-find: dbase ----
    def find_dbase(self, v):
        while self.dbase_par[v] != v:
            self.dbase_par[v] = self.dbase_par[self.dbase_par[v]]
            v = self.dbase_par[v]
        return v

    def union_dbase(self, a, b):
        a = self.find_dbase(a)
        b = self.find_dbase(b)
        if a != b:
            self.dbase_par[a] = b

    def make_rep_dbase(self, v):
        r = self.find_dbase(v)
        if r != v:
            self.dbase_par[r] = v
            self.dbase_par[v] = v

    # ---- union-find: dbase2 ----
    def find_db2(self, v):
        while self.db2_par[v] != v:
            self.db2_par[v] = self.db2_par[self.db2_par[v]]
            v = self.db2_par[v]
        return v

    def union_db2(self, a, b):
        a = self.find_db2(a)
        b = self.find_db2(b)
        if a != b:
            self.db2_par[a] = b

    def make_rep_db2(self, v):
        r = self.find_db2(v)
        if r != v:
            self.db2_par[r] = v
            self.db2_par[v] = v

    # ---- interleaved LCA ----
    def find_lca(self, u, v):
        self.lca_epoch += 1
        ep = self.lca_epoch
        hx = self.find_base(u)
        hy = self.find_base(v)
        self.lca_tag1[hx] = ep
        self.lca_tag2[hy] = ep
        while True:
            if self.lca_tag1[hy] == ep: return hy
            if self.lca_tag2[hx] == ep: return hx
            hxr = (self.mate[hx] == NIL or self.parent[self.mate[hx]] == NIL)
            hyr = (self.mate[hy] == NIL or self.parent[self.mate[hy]] == NIL)
            if hxr and hyr: return NIL
            if not hxr:
                hx = self.find_base(self.parent[self.mate[hx]])
                self.lca_tag1[hx] = ep
            if not hyr:
                hy = self.find_base(self.parent[self.mate[hy]])
                self.lca_tag2[hy] = ep

    # ---- shrink_path ----
    def shrink_path(self, b, x, y, dunions):
        v = self.find_base(x)
        while v != b:
            self.union_base(v, b, b)
            dunions.append((v, b))
            mv = self.mate[v]
            self.union_base(mv, b, b)
            dunions.append((mv, b))
            self.base_par[b] = b
            self.source_bridge[mv] = x
            self.target_bridge[mv] = y
            d = self.delta
            for w in self.graph[mv]:
                if w == self.mate[mv]: continue
                bw = self.find_base(w)
                if self.label[bw] == ODD: continue
                if self.label[bw] == UNLABELED:
                    self.level_queue[d + 1].append((mv, w))
                elif self.label[bw] == EVEN:
                    self.level_queue[d].append((mv, w))
            v = self.find_base(self.parent[mv])
        dunions.append((b, b))

    # ================================================================
    #                          PHASE 1
    # ================================================================
    def phase_1(self):
        self.delta = 0
        self.tree_nodes = []
        for q in self.level_queue:
            q.clear()
        dunions = []

        for i in range(self.n):
            self.base_par[i] = i
            self.dbase_par[i] = i
            self.label[i] = UNLABELED
            self.parent[i] = NIL
            self.source_bridge[i] = NIL
            self.target_bridge[i] = NIL
            self.in_tree[i] = False

        for v in range(self.n):
            if self.mate[v] == NIL:
                self.label[v] = EVEN
                self.in_tree[v] = True
                self.tree_nodes.append(v)
                for u in self.graph[v]:
                    if u == self.mate[v]: continue
                    bu = self.find_base(u)
                    if self.label[bu] == ODD: continue
                    if self.label[bu] == UNLABELED:
                        self.level_queue[1].append((v, u))
                    elif self.label[bu] == EVEN:
                        self.level_queue[0].append((v, u))

        found_sap = False

        while self.delta <= self.n:
            d = self.delta
            while self.level_queue[d]:
                z, u = self.level_queue[d].pop()
                bz = self.find_base(z)
                bu = self.find_base(u)
                if self.label[bz] != EVEN:
                    z, u = u, z
                    bz, bu = bu, bz
                if bz == bu or self.label[bz] != EVEN: continue
                if u == self.mate[z] or self.label[bu] == ODD: continue

                if self.label[bu] == UNLABELED:
                    mv = self.mate[u]
                    if mv == NIL: continue
                    self.parent[u] = z
                    self.parent[mv] = u
                    self.label[u] = ODD
                    self.label[mv] = EVEN
                    self.in_tree[u] = True
                    self.in_tree[mv] = True
                    self.tree_nodes.append(u)
                    self.tree_nodes.append(mv)
                    for w in self.graph[mv]:
                        if w == self.mate[mv]: continue
                        bw = self.find_base(w)
                        if self.label[bw] == ODD: continue
                        if self.label[bw] == UNLABELED:
                            self.level_queue[d + 1].append((mv, w))
                        elif self.label[bw] == EVEN:
                            self.level_queue[d].append((mv, w))

                elif self.label[bu] == EVEN:
                    lca = self.find_lca(z, u)
                    if lca != NIL:
                        self.shrink_path(lca, z, u, dunions)
                        self.shrink_path(lca, u, z, dunions)
                    else:
                        found_sap = True

            if found_sap:
                # Build H: mateH and h_adj
                for u in self.tree_nodes:
                    self.mate_h[u] = NIL
                    uh = self.find_dbase(u)
                    mv = self.mate[u]
                    if mv != NIL and self.in_tree[mv]:
                        vh = self.find_dbase(mv)
                        if uh != vh:
                            self.mate_h[uh] = vh
                            self.mate_h[vh] = uh
                # Build h_adj: non-matching edges between different dbase components
                for u in self.tree_nodes:
                    uh = self.find_dbase(u)
                    for w in self.graph[u]:
                        if not self.in_tree[w]: continue
                        if self.mate[u] == w: continue
                        wh = self.find_dbase(w)
                        if uh == wh: continue
                        self.h_adj[uh].append((u, w))
                return True

            for a, b in dunions:
                if a == b:
                    self.make_rep_dbase(a)
                else:
                    self.union_dbase(a, b)
            dunions.clear()
            self.delta += 1

        return False

    # ================================================================
    #                          PHASE 2
    # ================================================================

    def find_ap_hg(self, root_vh):
        """Iterative DFS in H using precomputed h_adj."""
        # Stack: [vh, edge_idx]
        stk = [[root_vh, 0]]

        while stk:
            f = stk[-1]
            vh = f[0]
            adj = self.h_adj[vh]

            found_next = False
            while f[1] < len(adj):
                v, w = adj[f[1]]
                f[1] += 1

                uh = self.find_db2(self.rep[w])
                if uh == self.find_db2(vh): continue
                if self.mate_h[vh] == uh: continue
                if self.label_h[uh] == ODD: continue

                if self.label_h[uh] == UNLABELED:
                    muh = self.mate_h[uh]
                    if muh == NIL:
                        self.label_h[uh] = ODD
                        self.parent_h_src[uh] = w
                        self.parent_h_tgt[uh] = v
                        return uh
                    self.label_h[uh] = ODD
                    self.parent_h_src[uh] = w
                    self.parent_h_tgt[uh] = v
                    self.label_h[muh] = EVEN
                    self.even_time_h[muh] = self.t_h
                    self.t_h += 1
                    stk.append([muh, 0])
                    found_next = True
                    break

                elif self.label_h[uh] == EVEN:
                    bh = self.find_db2(vh)
                    zh = self.find_db2(uh)
                    if self.even_time_h[bh] < self.even_time_h[zh]:
                        tmp = []
                        endpoints = []
                        cur = zh
                        while cur != bh:
                            endpoints.append(cur)
                            mc = self.mate_h[cur]
                            endpoints.append(mc)
                            tmp.append(mc)
                            ps = self.parent_h_src[mc]
                            pt = self.parent_h_tgt[mc]
                            nxt = self.rep[pt] if self.rep[ps] == mc else self.rep[ps]
                            cur = self.find_db2(nxt)
                        for nd in endpoints:
                            self.union_db2(nd, bh)
                        self.make_rep_db2(bh)
                        for mc in tmp:
                            self.bridge_h_src[mc] = v
                            self.bridge_h_tgt[mc] = w
                            self.dir_h[mc] = -1
                        for i in range(len(tmp) - 1, -1, -1):
                            stk.append([tmp[i], 0])
                        found_next = True
                        break

            if not found_next:
                stk.pop()

        return NIL

    def trace_h_path(self, vh, uh, edges_out):
        """Iterative trace from vh to uh in H, collecting non-matching G-edges."""
        stk = [[vh, uh, 0, 0, 0, 0, 0]]
        while stk:
            f = stk[-1]
            if f[0] == f[1]:
                stk.pop()
                continue
            if self.label_h[f[0]] == EVEN:
                mvh = self.mate_h[f[0]]
                ps = self.parent_h_src[mvh]
                pt = self.parent_h_tgt[mvh]
                edges_out.append((ps, pt))
                f[0] = self.rep[pt] if self.rep[ps] == mvh else self.rep[ps]
                continue
            if f[2] == 0:
                bs = self.bridge_h_src[f[0]]
                bt = self.bridge_h_tgt[f[0]]
                f[3] = bs
                f[4] = bt
                if self.dir_h[f[0]] == 1:
                    f[5] = self.rep[bs]
                    f[6] = self.rep[bt]
                else:
                    f[5] = self.rep[bt]
                    f[6] = self.rep[bs]
                f[2] = 1
                mt = self.rep[self.mate_h[f[0]]] if self.mate_h[f[0]] != NIL else f[0]
                stk.append([f[5], mt, 0, 0, 0, 0, 0])
                continue
            if f[2] == 1:
                edges_out.append((f[3], f[4]))
                f[2] = 2
                stk.append([f[6], f[1], 0, 0, 0, 0, 0])
                continue
            stk.pop()

    def find_path_in_g(self, v, u, pairs):
        """Iterative unfold within single H-node."""
        stk = [[v, u, 0, 0, 0]]
        while stk:
            f = stk[-1]
            if f[0] == f[1]:
                stk.pop()
                continue
            if f[2] == 0:
                if self.label[f[0]] == EVEN:
                    mv = self.mate[f[0]]
                    pmv = self.parent[mv]
                    pairs.append((mv, pmv))
                    f[0] = pmv
                    continue
                f[3] = self.source_bridge[f[0]]
                f[4] = self.target_bridge[f[0]]
                f[2] = 1
                stk.append([f[3], self.mate[f[0]], 0, 0, 0])
                continue
            if f[2] == 1:
                pairs.append((f[3], f[4]))
                f[2] = 2
                stk.append([f[4], f[1], 0, 0, 0])
                continue
            stk.pop()

    def augment_g(self, h_edges):
        """Unfold H-edges to G and augment."""
        pairs = []
        for u, v in h_edges:
            pairs.append((u, v))
            self.find_path_in_g(u, self.rep[u], pairs)
            self.find_path_in_g(v, self.rep[v], pairs)
        for a, b in pairs:
            self.mate[a] = b
            self.mate[b] = a

    def phase_2(self):
        """Find all SAPs in H, unfold and augment."""
        for v in self.tree_nodes:
            self.rep[v] = self.find_dbase(v)
            self.label_h[v] = UNLABELED
            self.parent_h_src[v] = NIL
            self.parent_h_tgt[v] = NIL
            self.bridge_h_src[v] = NIL
            self.bridge_h_tgt[v] = NIL
            self.dir_h[v] = 0
            self.even_time_h[v] = 0
            self.db2_par[v] = v
        self.t_h = 0

        all_paths = []
        for vh in self.tree_nodes:
            if vh != self.rep[vh]: continue
            if self.label_h[vh] != UNLABELED or self.mate_h[vh] != NIL: continue
            self.label_h[vh] = EVEN
            self.even_time_h[vh] = self.t_h
            self.t_h += 1

            free_node = self.find_ap_hg(vh)
            if free_node != NIL:
                h_nm = []
                ps = self.parent_h_src[free_node]
                pt = self.parent_h_tgt[free_node]
                h_nm.append((ps, pt))
                nxt = self.rep[pt] if self.rep[ps] == free_node else self.rep[ps]
                self.trace_h_path(nxt, vh, h_nm)
                all_paths.append(h_nm)

        for he in all_paths:
            self.augment_g(he)

        # Clean up
        for v in self.tree_nodes:
            db = self.find_dbase(v)
            self.h_adj[db].clear()
            self.h_adj[v].clear()
            self.mate_h[v] = NIL

    # ================================================================
    #                      MAIN ENTRY POINT
    # ================================================================
    def maximum_matching(self):
        # Greedy init
        for u in range(self.n):
            if self.mate[u] != NIL: continue
            for v in self.graph[u]:
                if self.mate[v] == NIL:
                    self.mate[u] = v
                    self.mate[v] = u
                    break
        while self.phase_1():
            self.phase_2()

        result = []
        for u in range(self.n):
            if self.mate[u] != NIL and self.mate[u] > u:
                result.append((u, self.mate[u]))
        result.sort()
        return result


# ================================================================
#                    VALIDATION AND MAIN
# ================================================================

def validate_matching(n, graph, matching):
    deg = [0] * n
    errors = 0
    for u, v in matching:
        idx = bisect_left(graph[u], v)
        if idx >= len(graph[u]) or graph[u][idx] != v:
            print(f"ERROR: Edge ({u}, {v}) not in graph!", file=sys.stderr)
            errors += 1
        deg[u] += 1
        deg[v] += 1
    for i in range(n):
        if deg[i] > 1:
            print(f"ERROR: Vertex {i} in {deg[i]} edges!", file=sys.stderr)
            errors += 1
    matched = sum(1 for d in deg if d > 0)
    print(f"\n=== Validation Report ===")
    print(f"Matching size: {len(matching)}")
    print(f"Matched vertices: {matched}")
    print("VALIDATION PASSED" if errors == 0 else "VALIDATION FAILED")
    print(f"=========================\n")


def load_graph(filename):
    with open(filename) as f:
        n, m = map(int, f.readline().split())
        edges = []
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                edges.append((int(parts[0]), int(parts[1])))
    return n, edges


def main():
    print("Gabow's Scaling Algorithm (Optimized V2) - Python Implementation")
    print("==================================================================\n")

    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <filename>")
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

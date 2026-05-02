"""
Micali-Vazirani Pure Algorithm - O(EâˆšV) Maximum Matching

True MV with DDFS, tenacity, regular + hanging bridges, petal contraction.
Faithful port of C++ micali_vazirani_pure.cpp.

All integers, no hash containers, fully deterministic.
"""

import sys
import time
from bisect import bisect_left

NIL = -1

# DDFS result codes
DDFS_EMPTY = 0
DDFS_PETAL = 1
DDFS_PATH  = 2


# =========================================================================
# Node
# =========================================================================
class Node:
    __slots__ = (
        'preds', 'pred_to', 'hanging_bridges',
        'min_level', 'max_level', 'even_level', 'odd_level',
        'match', 'bud', 'above', 'below',
        'ddfs_green', 'ddfs_red', 'number_preds',
        'deleted', 'visited',
    )

    def __init__(self):
        self.preds = []
        self.pred_to = []          # list of (target, index in target's preds)
        self.hanging_bridges = []
        self.min_level = NIL
        self.max_level = NIL
        self.even_level = NIL
        self.odd_level = NIL
        self.match = NIL
        self.bud = NIL
        self.above = NIL
        self.below = NIL
        self.ddfs_green = NIL
        self.ddfs_red = NIL
        self.number_preds = 0
        self.deleted = False
        self.visited = False

    def set_min_level(self, level):
        self.min_level = level
        if level % 2:
            self.odd_level = level
        else:
            self.even_level = level

    def set_max_level(self, level):
        self.max_level = level
        if level % 2:
            self.odd_level = level
        else:
            self.even_level = level

    def outer(self):
        return (self.even_level != NIL and
                (self.odd_level == NIL or self.even_level < self.odd_level))

    def inner(self):
        return not self.outer()

    def reset(self):
        self.preds.clear()
        self.pred_to.clear()
        self.hanging_bridges.clear()
        self.min_level = NIL
        self.max_level = NIL
        self.even_level = NIL
        self.odd_level = NIL
        self.bud = NIL
        self.above = NIL
        self.below = NIL
        self.ddfs_green = NIL
        self.ddfs_red = NIL
        self.number_preds = 0
        self.deleted = False
        self.visited = False


# =========================================================================
# DDFSResult
# =========================================================================
class DDFSResult:
    __slots__ = ('nodes_seen', 'bottleneck')

    def __init__(self):
        self.nodes_seen = []
        self.bottleneck = NIL


# =========================================================================
# MVGraph â€” the full algorithm
# =========================================================================
class MVGraph:
    def __init__(self):
        self.nodes = []
        self.edges = []          # flat adjacency (CSR values)
        self.adj_start = []
        self.deg = []
        self.levels = []
        self.bridges = []        # bridges by tenacity bucket
        self.green_stack = []
        self.red_stack = []
        self.path_found = []
        self.last_ddfs = DDFSResult()
        self.matchnum = 0
        self.bridgenum = 0
        self.todonum = 0

    # ---- construction ----
    def build(self, n, edge_list):
        self.nodes = [Node() for _ in range(n)]
        adj = [[] for _ in range(n)]
        for u, v in edge_list:
            if 0 <= u < n and 0 <= v < n and u != v:
                adj[u].append(v)
                adj[v].append(u)
        for i in range(n):
            adj[i].sort()

        self.adj_start = [0] * n
        self.deg = [0] * n
        self.edges = []
        for i in range(n):
            self.adj_start[i] = len(self.edges)
            self.deg[i] = len(adj[i])
            self.edges.extend(adj[i])

    # ---- greedy initialization ----
    greedy_size = 0

    def greedy_init(self):
        n = len(self.nodes)
        cnt = 0
        for j in range(n):
            if self.nodes[j].match == NIL:
                for k in range(self.deg[j]):
                    i = self.edges[self.adj_start[j] + k]
                    if self.nodes[i].match == NIL:
                        self.nodes[j].match = i
                        self.nodes[i].match = j
                        self.matchnum += 1
                        cnt += 1
                        break
        return cnt

    def greedy_init_md(self):
        n = len(self.nodes)
        cnt = 0
        order = sorted(range(n), key=lambda x: (self.deg[x], x))
        for j in order:
            if self.nodes[j].match == NIL:
                best = NIL
                best_deg = float('inf')
                for k in range(self.deg[j]):
                    i = self.edges[self.adj_start[j] + k]
                    if self.nodes[i].match == NIL and self.deg[i] < best_deg:
                        best = i
                        best_deg = self.deg[i]
                if best != NIL:
                    self.nodes[j].match = best
                    self.nodes[best].match = j
                    self.matchnum += 1
                    cnt += 1
        return cnt

    # ---- helpers ----
    def add_to_level(self, level, node):
        while level >= len(self.levels):
            self.levels.append([])
        self.levels[level].append(node)
        self.todonum += 1

    def add_to_bridges(self, level, n1, n2):
        while level >= len(self.bridges):
            self.bridges.append([])
        self.bridges[level].append((n1, n2))
        self.bridgenum += 1

    def tenacity(self, n1, n2):
        if self.nodes[n1].match == n2:  # matched bridge
            if self.nodes[n1].odd_level != NIL and self.nodes[n2].odd_level != NIL:
                return self.nodes[n1].odd_level + self.nodes[n2].odd_level + 1
        else:  # unmatched bridge
            if self.nodes[n1].even_level != NIL and self.nodes[n2].even_level != NIL:
                return self.nodes[n1].even_level + self.nodes[n2].even_level + 1
        return NIL

    def bud_star(self, c):
        b = self.nodes[c].bud
        if b == NIL:
            return c
        return self.bud_star(b)

    def bud_star_includes(self, c, goal):
        if c == goal:
            return True
        b = self.nodes[c].bud
        if b == NIL:
            return False
        return self.bud_star_includes(b, goal)

    # ---- reset between phases ----
    def phase_reset(self):
        for v in self.levels:
            v.clear()
        for v in self.bridges:
            v.clear()
        self.bridgenum = 0
        self.todonum = 0
        n = len(self.nodes)
        for i in range(n):
            self.nodes[i].reset()
            if self.nodes[i].match == NIL:
                self.add_to_level(0, i)
                self.nodes[i].set_min_level(0)

    # ---- step_to: core level-building step ----
    def step_to(self, to, frm, level):
        level += 1
        tl = self.nodes[to].min_level
        if tl == NIL or tl >= level:
            if tl != level:
                self.add_to_level(level, to)
                self.nodes[to].set_min_level(level)
            self.nodes[to].preds.append(frm)
            self.nodes[to].number_preds += 1
            self.nodes[frm].pred_to.append((to, len(self.nodes[to].preds) - 1))
        else:
            # found a bridge
            ten = self.tenacity(to, frm)
            if ten == NIL:
                self.nodes[to].hanging_bridges.append(frm)
                self.nodes[frm].hanging_bridges.append(to)
            else:
                self.add_to_bridges((ten - 1) // 2, to, frm)

    # ---- MIN phase ----
    def MIN(self, i):
        if i >= len(self.levels):
            return
        k = 0
        while k < len(self.levels[i]):
            current = self.levels[i][k]
            k += 1
            self.todonum -= 1
            nd = self.nodes[current]
            if i % 2 == 0:
                for j in range(self.deg[current]):
                    edge = self.edges[self.adj_start[current] + j]
                    if edge != nd.match:
                        self.step_to(edge, current, i)
            else:
                if nd.match != NIL:
                    self.step_to(nd.match, current, i)

    # ---- MAX phase ----
    def MAX(self, i):
        found = False
        if i >= len(self.bridges):
            return False

        j = 0
        while j < len(self.bridges[i]):
            n1, n2 = self.bridges[i][j]
            j += 1
            self.bridgenum -= 1
            if self.nodes[n1].deleted or self.nodes[n2].deleted:
                continue

            result = self.DDFS(n1, n2)
            if result == DDFS_EMPTY:
                continue

            if result == DDFS_PATH:
                self.find_path(n1, n2)
                self.augment_path()
                if len(self.nodes) // 2 <= self.matchnum:
                    return True
                self.remove_path()
                found = True
            else:  # DDFS_PETAL
                b = self.last_ddfs.bottleneck
                current_ten = i * 2 + 1
                for itt in list(self.last_ddfs.nodes_seen):
                    self.nodes[itt].bud = b
                    self.nodes[itt].set_max_level(current_ten - self.nodes[itt].min_level)
                    self.add_to_level(self.nodes[itt].max_level, itt)
                    for hanging in self.nodes[itt].hanging_bridges:
                        hanging_ten = self.tenacity(itt, hanging)
                        if hanging_ten != NIL:
                            self.add_to_bridges((hanging_ten - 1) // 2, itt, hanging)
        return found

    # ==================================================================
    # DDFS â€” Double Depth-First Search
    # ==================================================================

    def add_pred_to_stack(self, cur, stack):
        for pred in self.nodes[cur].preds:
            if pred != NIL:
                stack.append((cur, pred))

    def prepare_next(self, nx):
        """nx is [first, second]; mutated in place."""
        if nx[0] != NIL:
            self.nodes[nx[0]].below = nx[1]
        nx[1] = self.bud_star(nx[1])

    @staticmethod
    def edge_valid(e):
        return not (e[0] == NIL and e[1] == NIL)

    @staticmethod
    def node_from_stack(e, S):
        """Pop from S into e (mutate list e in place)."""
        if S:
            top = S.pop()
            e[0] = top[0]
            e[1] = top[1]
        else:
            e[0] = NIL
            e[1] = NIL

    def L(self, e):
        n = self.bud_star(e[1])
        return self.nodes[n].min_level

    def step_into(self, C_ref, nx, S, green_top, red_top):
        """C_ref is [value]; nx is [first, second]. Mutated in place."""
        self.prepare_next(nx)
        if not self.nodes[nx[1]].visited:
            self.nodes[nx[1]].above = nx[0]
            C_ref[0] = nx[1]
            nd = self.nodes[C_ref[0]]
            nd.visited = True
            nd.ddfs_green = green_top
            nd.ddfs_red = red_top
            self.last_ddfs.nodes_seen.append(C_ref[0])
            self.add_pred_to_stack(C_ref[0], S)
        self.node_from_stack(nx, S)

    def DDFS(self, green_top, red_top):
        self.last_ddfs.nodes_seen.clear()
        self.last_ddfs.bottleneck = NIL

        Sr = self.red_stack
        Sg = self.green_stack
        Sr.clear()
        Sg.clear()

        G = [NIL]
        R = [NIL]

        if self.bud_star(red_top) == self.bud_star(green_top):
            return DDFS_EMPTY
        if self.nodes[green_top].min_level == 0 and self.nodes[red_top].min_level == 0:
            return DDFS_PATH

        Ng = [NIL, green_top]
        Nr = [NIL, red_top]
        red_before = [NIL, NIL]
        green_before = [NIL, NIL]

        while (R[0] == NIL or G[0] == NIL or
               self.nodes[R[0]].min_level > 0 or self.nodes[G[0]].min_level > 0):

            while self.edge_valid(Nr) and self.edge_valid(Ng) and self.L(Nr) != self.L(Ng):

                while self.edge_valid(Nr) and self.L(Nr) > self.L(Ng):
                    self.step_into(R, Nr, Sr, green_top, red_top)

                if not self.edge_valid(Nr):
                    Nr[0] = red_before[0]
                    Nr[1] = red_before[1]
                    tmp = red_before[0]
                    while self.nodes[tmp].above != NIL:
                        rc = self.nodes[tmp].above
                        for ri in self.nodes[rc].preds:
                            if ri == NIL:
                                continue
                            if self.bud_star(ri) == tmp:
                                self.nodes[rc].below = ri
                                break
                        tmp = self.nodes[tmp].above

                while self.edge_valid(Ng) and self.L(Nr) < self.L(Ng):
                    self.step_into(G, Ng, Sg, green_top, red_top)

                if not self.edge_valid(Ng):
                    Ng[0] = green_before[0]
                    Ng[1] = green_before[1]
                    tmp = green_before[0]
                    while self.nodes[tmp].above != NIL:
                        rc = self.nodes[tmp].above
                        for ri in self.nodes[rc].preds:
                            if ri == NIL:
                                continue
                            if self.bud_star(ri) == tmp:
                                self.nodes[rc].below = ri
                                break
                        tmp = self.nodes[tmp].above

            if self.bud_star(Nr[1]) == self.bud_star(Ng[1]):
                if Sr:
                    red_before[0] = Nr[0]
                    red_before[1] = Nr[1]
                    self.prepare_next(Nr)
                    self.node_from_stack(Nr, Sr)
                    if self.edge_valid(Nr):
                        R[0] = Nr[0]
                    else:
                        Nr[0] = red_before[0]
                        Nr[1] = red_before[1]
                elif Sg:
                    green_before[0] = Ng[0]
                    green_before[1] = Ng[1]
                    self.prepare_next(Ng)
                    self.node_from_stack(Ng, Sg)
                    if self.edge_valid(Ng):
                        G[0] = Ng[0]
                    else:
                        Ng[0] = green_before[0]
                        Ng[1] = green_before[1]
                else:
                    self.prepare_next(Nr)
                    self.prepare_next(Ng)
                    self.last_ddfs.bottleneck = Nr[1]
                    return DDFS_PETAL
            else:
                self.step_into(R, Nr, Sr, green_top, red_top)
                self.step_into(G, Ng, Sg, green_top, red_top)

        return DDFS_PATH

    # ==================================================================
    # Path finding and augmentation
    # ==================================================================

    def find_path(self, n1, n2):
        self.path_found.clear()
        self.walk_down_path(n1)
        self.path_found.reverse()
        self.walk_down_path(n2)

    def walk_down_path(self, start):
        cur = start
        while cur != NIL:
            if self.nodes[cur].bud != NIL:
                cur = self.walk_blossom(cur)
            else:
                self.path_found.append(cur)
                cur = self.nodes[cur].below

    def jump_bridge(self, cur):
        if self.nodes[cur].ddfs_green == cur:
            return self.nodes[cur].ddfs_red
        if self.nodes[cur].ddfs_red == cur:
            return self.nodes[cur].ddfs_green
        if self.bud_star_includes(self.nodes[cur].ddfs_green, cur):
            before = len(self.path_found)
            b = self.nodes[cur].ddfs_green
            while b != cur:
                b = self.walk_blossom(b)
            self.path_found[before:] = self.path_found[before:][::-1]
            return self.nodes[cur].ddfs_red
        else:
            before = len(self.path_found)
            b = self.nodes[cur].ddfs_red
            while b != cur:
                b = self.walk_blossom(b)
            self.path_found[before:] = self.path_found[before:][::-1]
            return self.nodes[cur].ddfs_green

    def walk_blossom(self, cur):
        if self.nodes[cur].outer():
            cur = self.walk_blossom_down(cur, NIL)
        else:
            cur = self.walk_blossom_up(cur)
            before = cur
            cur = self.jump_bridge(cur)
            cur = self.walk_blossom_down(cur, before)
        return cur

    def walk_blossom_down(self, cur, before):
        if before == NIL:
            before = cur
        b = self.nodes[cur].bud
        while cur != NIL and cur != b:
            if (self.nodes[cur].ddfs_green != self.nodes[before].ddfs_green or
                    self.nodes[cur].ddfs_red != self.nodes[before].ddfs_red):
                cur = self.walk_blossom(cur)
            else:
                self.path_found.append(cur)
                cur = self.nodes[cur].below
        return cur

    def walk_blossom_up(self, cur):
        while True:
            self.path_found.append(cur)
            if self.nodes[cur].above == NIL:
                break
            b = self.nodes[self.nodes[cur].above].below
            if b != cur and self.bud_star_includes(b, cur):
                before = len(self.path_found)
                while b != cur:
                    b = self.walk_blossom(b)
                self.path_found[before:] = self.path_found[before:][::-1]
            cur = self.nodes[cur].above
        return cur

    def augment_path(self):
        i = 0
        while i + 1 < len(self.path_found):
            n1 = self.path_found[i]
            n2 = self.path_found[i + 1]
            self.nodes[n1].match = n2
            self.nodes[n2].match = n1
            i += 2
        self.matchnum += 1

    def remove_path(self):
        while self.path_found:
            current = self.path_found.pop()
            if not self.nodes[current].deleted:
                self.nodes[current].deleted = True
                for tgt, idx in self.nodes[current].pred_to:
                    nd = self.nodes[tgt]
                    if not nd.deleted:
                        nd.preds[idx] = NIL
                        nd.number_preds -= 1
                        if nd.number_preds <= 0:
                            self.path_found.append(tgt)

    # ---- main matching driver ----
    def max_match(self):
        n = len(self.nodes)
        for i in range(n):
            if self.nodes[i].match == NIL:
                self.add_to_level(0, i)
                self.nodes[i].set_min_level(0)

        found = self.max_match_phase()
        while n // 2 > self.matchnum and found:
            self.phase_reset()
            found = self.max_match_phase()

    def max_match_phase(self):
        n = len(self.nodes)
        found = False
        for i in range(n // 2 + 1):
            if found:
                break
            if self.todonum <= 0 and self.bridgenum <= 0:
                return False
            self.MIN(i)
            found = self.MAX(i)
        return found

    def get_matching(self):
        result = []
        n = len(self.nodes)
        for i in range(n):
            if self.nodes[i].match != NIL and self.nodes[i].match > i:
                result.append((i, self.nodes[i].match))
        return result


# =========================================================================
# Validation
# =========================================================================
def validate_matching(n, edges_flat, adj_start, deg, matching):
    vdeg = [0] * n
    errors = 0

    for u, v in matching:
        # binary search in sorted adjacency
        start = adj_start[u]
        end = start + deg[u]
        pos = bisect_left(edges_flat, v, start, end)
        if pos >= end or edges_flat[pos] != v:
            print(f"ERROR: Edge ({u}, {v}) not in graph!", file=sys.stderr)
            errors += 1
        vdeg[u] += 1
        vdeg[v] += 1

    for i in range(n):
        if vdeg[i] > 1:
            print(f"ERROR: Vertex {i} in {vdeg[i]} edges!", file=sys.stderr)
            errors += 1

    matched = sum(1 for d in vdeg if d > 0)

    print()
    print("=== Validation Report ===")
    print(f"Matching size: {len(matching)}")
    print(f"Matched vertices: {matched}")
    print("VALIDATION FAILED" if errors > 0 else "VALIDATION PASSED")
    print("=========================")
    print()


# =========================================================================
# File I/O and main
# =========================================================================
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
    print("Micali-Vazirani Pure Algorithm - Python Implementation")
    print("=======================================================")
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

    n, edge_list = load_graph(sys.argv[1])
    print(f"Graph: {n} vertices, {len(edge_list)} edges")

    t0 = time.time()
    mv = MVGraph()
    mv.build(n, edge_list)
    if greedy_mode == 1:
        mv.greedy_size = mv.greedy_init()
    elif greedy_mode == 2:
        mv.greedy_size = mv.greedy_init_md()
    mv.max_match()
    t1 = time.time()

    matching = mv.get_matching()

    validate_matching(n, mv.edges, mv.adj_start, mv.deg, matching)

    print(f"Matching size: {len(matching)}")
    if greedy_mode > 0:
        print(f"Greedy init size: {mv.greedy_size}")
        if matching:
            print(f"Greedy/Final: {100.0 * mv.greedy_size / len(matching):.2f}%")
        else:
            print("Greedy/Final: NA")
    print(f"Time: {int((t1 - t0) * 1000)} ms")


if __name__ == "__main__":
    main()

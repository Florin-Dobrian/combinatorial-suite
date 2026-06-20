/*
 * Hopcroft-Karp Bidirectional, Iterative, with DDFS (HKB2) - O(E sqrt(V))
 *
 * Variant 2 ("BMV"): bidirectional BFS interleaved with DDFS at the level
 * where bridges first appear. The DDFS coordination (coordinated cursors,
 * Sr/Sg descent stacks, level-tie tie-break via L(), step_into, above/below
 * pointers) is borrowed from MV. In bipartite there are no blossoms, so the
 * contracted graph H of MV equals G: bud_star is identity, only min_level
 * matters, hanging bridges / petals never fire, and a same-level bridge at
 * common level i has tenacity 2i+1 -> bucket index i.
 *
 * BFS predecessors are recovered directly from the CSR via the parity rule
 * (even level: unique pred = mate; odd level: non-matching neighbors at
 * level-1), so only num_preds (a counter driving cascade-delete) is kept per
 * node. The algorithm runs in a unified vertex namespace [0, sN+tN): S in
 * [0, sN), T in [sN, sN+tN). Top-level returns num_phases.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;

const DDFS_EMPTY: i32 = 0;
const DDFS_PATH: i32 = 2;
/* DDFS_PETAL omitted -- never returned in bipartite. */

/* ---------- Input: BipartiteGraph ---------- */

#[allow(dead_code)]
struct BipartiteGraph {
    s_num_vtxs: usize,
    t_num_vtxs: usize,
    num_edgs: usize,
    s_idx: Vec<usize>,
    t_idx: Vec<usize>,
    s_adj: Vec<i32>,
    t_adj: Vec<i32>,
}

/* ---------- Output: BipartiteMatching ---------- */

#[allow(dead_code)]
struct BipartiteMatching {
    s_num_vtxs: usize,
    t_num_vtxs: usize,
    num_edgs: usize,
    s_mate: Vec<i32>,
    t_mate: Vec<i32>,
}

/* ---------- Per-vertex node ---------- */

#[allow(dead_code)]
struct Hkb2Node {
    min_level: i32,
    even_level: i32,
    odd_level: i32,
    mate: i32,        // 'match' in C++ (reserved word in Rust)
    above: i32,
    below: i32,
    ddfs_green: i32,
    ddfs_red: i32,
    num_preds: i32,   // count of still-alive BFS predecessors
    deleted: u8,
    visited: u8,
}

impl Hkb2Node {
    fn new() -> Self {
        Hkb2Node {
            min_level: NIL, even_level: NIL, odd_level: NIL, mate: NIL,
            above: NIL, below: NIL, ddfs_green: NIL, ddfs_red: NIL,
            num_preds: 0, deleted: 0, visited: 0,
        }
    }
    fn set_min_level(&mut self, lvl: i32) {
        self.min_level = lvl;
        if lvl & 1 != 0 { self.odd_level = lvl; } else { self.even_level = lvl; }
    }
    fn reset(&mut self) {
        self.min_level = NIL; self.even_level = NIL; self.odd_level = NIL;
        self.above = NIL; self.below = NIL; self.ddfs_green = NIL; self.ddfs_red = NIL;
        self.num_preds = 0; self.deleted = 0; self.visited = 0;
        /* mate intentionally NOT reset -- persists across phases. */
    }
}

/* ---------- Internal state ---------- */

struct Hkb2State {
    n_vtxs: usize,
    s_num_vtxs: usize,

    v_idx: Vec<usize>,        // unified CSR, both directions
    v_adj: Vec<i32>,

    nodes: Vec<Hkb2Node>,

    levels: Vec<Vec<i32>>,
    br_buckets: Vec<Vec<(i32, i32)>>,

    s_stk_r: Vec<(i32, i32)>, // red descent stack (reused across DDFS calls)
    s_stk_g: Vec<(i32, i32)>, // green descent stack
    path: Vec<i32>,           // reconstructed AP

    todo_num: i32,
    bridge_num: i32,
}

/* ---------- Free helpers for DDFS ---------- */

fn edge_valid(e: (i32, i32)) -> bool {
    !(e.0 == NIL && e.1 == NIL)
}

fn node_from_stack(e: &mut (i32, i32), s: &mut Vec<(i32, i32)>) {
    match s.pop() {
        Some(x) => *e = x,
        None => *e = (NIL, NIL),
    }
}

/* ---------- HKB2 algorithm (operates purely on Hkb2State) ---------- */

impl Hkb2State {
    fn add_to_level(&mut self, lvl: i32, v: i32) {
        if lvl as usize >= self.levels.len() {
            self.levels.resize_with((lvl + 1) as usize, Vec::new);
        }
        self.levels[lvl as usize].push(v);
        self.todo_num += 1;
    }

    fn add_to_bridges(&mut self, bucket: i32, u: i32, v: i32) {
        if bucket as usize >= self.br_buckets.len() {
            self.br_buckets.resize_with((bucket + 1) as usize, Vec::new);
        }
        self.br_buckets[bucket as usize].push((u, v));
        self.bridge_num += 1;
    }

    /* step_to: prop step or same-level bridge. Bucket index = common level. */
    fn step_to(&mut self, to: i32, from: i32, level: i32) {
        let new_level = level + 1;
        let tl = self.nodes[to as usize].min_level;
        if tl == NIL || tl >= new_level {
            if tl != new_level {
                self.add_to_level(new_level, to);
                self.nodes[to as usize].set_min_level(new_level);
            }
            self.nodes[to as usize].num_preds += 1;
        } else if tl == level {
            self.add_to_bridges(level, to, from);
        }
        /* else: back-edge (tl < level), silently dropped */
    }

    fn min_pass(&mut self, i: i32) {
        if i as usize >= self.levels.len() { return; }
        let sz = self.levels[i as usize].len();
        for k in 0..sz {
            let cur = self.levels[i as usize][k];
            self.todo_num -= 1;
            if i & 1 == 0 {
                let cur_mate = self.nodes[cur as usize].mate;
                let b = self.v_idx[cur as usize];
                let e = self.v_idx[cur as usize + 1];
                for j in b..e {
                    let edge = self.v_adj[j];
                    if edge != cur_mate { self.step_to(edge, cur, i); }
                }
            } else {
                let cur_mate = self.nodes[cur as usize].mate;
                if cur_mate != NIL { self.step_to(cur_mate, cur, i); }
            }
        }
    }

    /* ---- DDFS ---- */

    /* Enumerate live BFS predecessors of `cur`, pushing (cur, pred) pairs. */
    fn add_preds_to_stack(&self, cur: i32, stk: &mut Vec<(i32, i32)>) {
        let n = &self.nodes[cur as usize];
        if n.num_preds <= 0 { return; }
        let target_lvl = n.min_level - 1;
        if n.min_level & 1 == 0 {
            let m = n.mate;
            if m != NIL && self.nodes[m as usize].deleted == 0
                && self.nodes[m as usize].min_level == target_lvl {
                stk.push((cur, m));
            }
        } else {
            let mate = n.mate;
            let b = self.v_idx[cur as usize];
            let e = self.v_idx[cur as usize + 1];
            for j in b..e {
                let nb = self.v_adj[j];
                if nb == mate { continue; }
                if self.nodes[nb as usize].deleted != 0 { continue; }
                if self.nodes[nb as usize].min_level != target_lvl { continue; }
                stk.push((cur, nb));
            }
        }
    }

    fn prepare_next(&mut self, nx: (i32, i32)) {
        if nx.0 != NIL { self.nodes[nx.0 as usize].below = nx.1; }
    }

    fn level_of(&self, e: (i32, i32)) -> i32 {
        self.nodes[e.1 as usize].min_level
    }

    fn step_into(&mut self, c: &mut i32, nx: &mut (i32, i32),
                 s: &mut Vec<(i32, i32)>, green_top: i32, red_top: i32) {
        self.prepare_next(*nx);
        if self.nodes[nx.1 as usize].visited == 0 {
            let cc = nx.1;
            self.nodes[nx.1 as usize].above = nx.0;
            *c = cc;
            {
                let n = &mut self.nodes[cc as usize];
                n.visited = 1;
                n.ddfs_green = green_top;
                n.ddfs_red = red_top;
            }
            self.add_preds_to_stack(cc, s);
        }
        node_from_stack(nx, s);
    }

    fn ddfs(&mut self, green_top: i32, red_top: i32) -> i32 {
        let mut sr = std::mem::take(&mut self.s_stk_r);
        let mut sg = std::mem::take(&mut self.s_stk_g);
        sr.clear();
        sg.clear();
        let r = self.ddfs_inner(&mut sr, &mut sg, green_top, red_top);
        self.s_stk_r = sr;
        self.s_stk_g = sg;
        r
    }

    fn ddfs_inner(&mut self, sr: &mut Vec<(i32, i32)>, sg: &mut Vec<(i32, i32)>,
                  green_top: i32, red_top: i32) -> i32 {
        let mut r_cur: i32 = NIL;
        let mut g_cur: i32 = NIL;

        if red_top == green_top { return DDFS_EMPTY; }
        if self.nodes[green_top as usize].min_level == 0
            && self.nodes[red_top as usize].min_level == 0 {
            return DDFS_PATH;
        }

        let mut ng: (i32, i32) = (NIL, green_top);
        let mut nr: (i32, i32) = (NIL, red_top);
        let mut red_before: (i32, i32) = (NIL, NIL);
        let mut green_before: (i32, i32) = (NIL, NIL);

        while r_cur == NIL || g_cur == NIL
            || self.nodes[r_cur as usize].min_level > 0
            || self.nodes[g_cur as usize].min_level > 0 {

            while edge_valid(nr) && edge_valid(ng) && self.level_of(nr) != self.level_of(ng) {

                while edge_valid(nr) && self.level_of(nr) > self.level_of(ng) {
                    self.step_into(&mut r_cur, &mut nr, sr, green_top, red_top);
                }

                if !edge_valid(nr) {
                    nr = red_before;
                    let mut tmp = red_before.0;
                    while tmp != NIL && self.nodes[tmp as usize].above != NIL {
                        let rc = self.nodes[tmp as usize].above;
                        if self.nodes[tmp as usize].deleted == 0 {
                            self.nodes[rc as usize].below = tmp;
                        }
                        tmp = self.nodes[tmp as usize].above;
                    }
                }

                while edge_valid(ng) && self.level_of(nr) < self.level_of(ng) {
                    self.step_into(&mut g_cur, &mut ng, sg, green_top, red_top);
                }

                if !edge_valid(ng) {
                    ng = green_before;
                    let mut tmp = green_before.0;
                    while tmp != NIL && self.nodes[tmp as usize].above != NIL {
                        let rc = self.nodes[tmp as usize].above;
                        if self.nodes[tmp as usize].deleted == 0 {
                            self.nodes[rc as usize].below = tmp;
                        }
                        tmp = self.nodes[tmp as usize].above;
                    }
                }
            }

            if nr.1 == ng.1 {
                /* Cursors collide. Backtrack red first, then green. */
                if !sr.is_empty() {
                    red_before = nr;
                    self.prepare_next(nr);
                    node_from_stack(&mut nr, sr);
                    if edge_valid(nr) { r_cur = nr.0; } else { nr = red_before; }
                } else if !sg.is_empty() {
                    green_before = ng;
                    self.prepare_next(ng);
                    node_from_stack(&mut ng, sg);
                    if edge_valid(ng) { g_cur = ng.0; } else { ng = green_before; }
                } else {
                    /* Bipartite-dead branch (would be DDFS_PETAL in MV). */
                    self.prepare_next(nr);
                    self.prepare_next(ng);
                    return DDFS_EMPTY;
                }
            } else {
                self.step_into(&mut r_cur, &mut nr, sr, green_top, red_top);
                self.step_into(&mut g_cur, &mut ng, sg, green_top, red_top);
            }
        }
        DDFS_PATH
    }

    /* find_path / walk_down_path: flat below-pointer descent (bud is NIL). */
    fn find_path(&mut self, n1: i32, n2: i32) {
        self.path.clear();
        self.walk_down_path(n1);
        self.path.reverse();
        self.walk_down_path(n2);
    }

    fn walk_down_path(&mut self, start: i32) {
        let mut cur = start;
        while cur != NIL {
            self.path.push(cur);
            cur = self.nodes[cur as usize].below;
        }
    }

    fn augment_path(&mut self) {
        let mut i = 0;
        while i + 1 < self.path.len() {
            let n1 = self.path[i];
            let n2 = self.path[i + 1];
            self.nodes[n1 as usize].mate = n2;
            self.nodes[n2 as usize].mate = n1;
            i += 2;
        }
    }

    /* removePath: mark AP vertices deleted, cascade-delete preds reaching 0. */
    fn remove_path(&mut self) {
        let mut stk = std::mem::take(&mut self.path);
        while let Some(cur) = stk.pop() {
            if self.nodes[cur as usize].deleted != 0 { continue; }
            self.nodes[cur as usize].deleted = 1;
            let cur_lvl = self.nodes[cur as usize].min_level;
            if cur_lvl == NIL { continue; }
            let child_lvl = cur_lvl + 1;

            if cur_lvl & 1 == 1 {
                /* odd -> unique matching child at even level */
                let m = self.nodes[cur as usize].mate;
                if m != NIL && self.nodes[m as usize].deleted == 0
                    && self.nodes[m as usize].min_level == child_lvl {
                    self.nodes[m as usize].num_preds -= 1;
                    if self.nodes[m as usize].num_preds <= 0 { stk.push(m); }
                }
            } else {
                /* even -> non-matching children at odd level */
                let mate = self.nodes[cur as usize].mate;
                let b = self.v_idx[cur as usize];
                let e = self.v_idx[cur as usize + 1];
                for j in b..e {
                    let nb = self.v_adj[j];
                    if nb == mate { continue; }
                    if self.nodes[nb as usize].deleted != 0 { continue; }
                    if self.nodes[nb as usize].min_level != child_lvl { continue; }
                    self.nodes[nb as usize].num_preds -= 1;
                    if self.nodes[nb as usize].num_preds <= 0 { stk.push(nb); }
                }
            }
        }
        self.path = stk; /* restore buffer (logically empty), reuse capacity */
    }

    fn max_pass(&mut self, i: i32) -> bool {
        let mut found = false;
        if i as usize >= self.br_buckets.len() { return false; }
        let len = self.br_buckets[i as usize].len();
        for j in 0..len {
            let (n1, n2) = self.br_buckets[i as usize][j];
            self.bridge_num -= 1;
            if self.nodes[n1 as usize].deleted != 0 || self.nodes[n2 as usize].deleted != 0 {
                continue;
            }
            let result = self.ddfs(n1, n2);
            if result == DDFS_EMPTY { continue; }
            /* DDFS_PATH */
            self.find_path(n1, n2);
            self.augment_path();
            self.remove_path();
            found = true;
        }
        found
    }

    fn reset_phase(&mut self) {
        for v in &mut self.levels { v.clear(); }
        for v in &mut self.br_buckets { v.clear(); }
        self.bridge_num = 0;
        self.todo_num = 0;
        for i in 0..self.n_vtxs {
            self.nodes[i].reset();
            if self.nodes[i].mate == NIL {
                self.add_to_level(0, i as i32);
                self.nodes[i].set_min_level(0);
            }
        }
    }

    fn run_phase(&mut self) -> bool {
        let mut found = false;
        let limit = self.n_vtxs as i32 / 2 + 1;
        let mut i = 0;
        while i < limit && !found {
            if self.todo_num <= 0 && self.bridge_num <= 0 { return false; }
            self.min_pass(i);
            found = self.max_pass(i);
            i += 1;
        }
        found
    }
}

/* ---------- BipartiteGraph / BipartiteMatching construction ---------- */

fn build_bipartite_graph(s_num_vtxs: usize, t_num_vtxs: usize,
                         edges: &[(usize, usize)]) -> BipartiteGraph {
    let mut s_tmp: Vec<Vec<i32>> = vec![Vec::new(); s_num_vtxs];
    let mut t_tmp: Vec<Vec<i32>> = vec![Vec::new(); t_num_vtxs];
    for &(u, v) in edges {
        if u < s_num_vtxs && v < t_num_vtxs {
            s_tmp[u].push(v as i32);
            t_tmp[v].push(u as i32);
        }
    }
    for adj in &mut s_tmp { adj.sort_unstable(); adj.dedup(); }
    for adj in &mut t_tmp { adj.sort_unstable(); adj.dedup(); }

    let mut s_idx = vec![0usize; s_num_vtxs + 1];
    for u in 0..s_num_vtxs { s_idx[u + 1] = s_idx[u] + s_tmp[u].len(); }
    let mut s_adj: Vec<i32> = Vec::with_capacity(s_idx[s_num_vtxs]);
    for u in 0..s_num_vtxs { s_adj.extend_from_slice(&s_tmp[u]); }

    let mut t_idx = vec![0usize; t_num_vtxs + 1];
    for v in 0..t_num_vtxs { t_idx[v + 1] = t_idx[v] + t_tmp[v].len(); }
    let mut t_adj: Vec<i32> = Vec::with_capacity(t_idx[t_num_vtxs]);
    for v in 0..t_num_vtxs { t_adj.extend_from_slice(&t_tmp[v]); }

    let num_edgs = s_idx[s_num_vtxs];
    BipartiteGraph { s_num_vtxs, t_num_vtxs, num_edgs, s_idx, t_idx, s_adj, t_adj }
}

fn empty_bipartite_matching(g: &BipartiteGraph) -> BipartiteMatching {
    BipartiteMatching {
        s_num_vtxs: g.s_num_vtxs,
        t_num_vtxs: g.t_num_vtxs,
        num_edgs: 0,
        s_mate: vec![NIL; g.s_num_vtxs],
        t_mate: vec![NIL; g.t_num_vtxs],
    }
}

/* ---------- Greedy initial matching ---------- */

fn greedy_init(g: &BipartiteGraph, m: &mut BipartiteMatching) -> usize {
    let mut cnt: usize = 0;
    for u in 0..g.s_num_vtxs {
        if m.s_mate[u] != NIL { continue; }
        let s = g.s_idx[u];
        let e = g.s_idx[u + 1];
        for j in s..e {
            let v = g.s_adj[j];
            if m.t_mate[v as usize] == NIL {
                m.s_mate[u] = v;
                m.t_mate[v as usize] = u as i32;
                cnt += 1;
                break;
            }
        }
    }
    m.num_edgs += cnt;
    cnt
}

fn greedy_init_md(g: &BipartiteGraph, m: &mut BipartiteMatching) -> usize {
    let mut cnt: usize = 0;
    let mut deg = vec![0usize; g.t_num_vtxs];
    for u in 0..g.s_num_vtxs {
        let s = g.s_idx[u];
        let e = g.s_idx[u + 1];
        for j in s..e { deg[g.s_adj[j] as usize] += 1; }
    }
    let mut order: Vec<usize> = (0..g.s_num_vtxs).collect();
    order.sort_unstable_by(|&a, &b| {
        let da = g.s_idx[a + 1] - g.s_idx[a];
        let db = g.s_idx[b + 1] - g.s_idx[b];
        da.cmp(&db).then(a.cmp(&b))
    });
    for u in order {
        if m.s_mate[u] != NIL { continue; }
        let mut best: i32 = NIL;
        let mut best_deg = usize::MAX;
        let s = g.s_idx[u];
        let e = g.s_idx[u + 1];
        for j in s..e {
            let v = g.s_adj[j];
            if m.t_mate[v as usize] == NIL && deg[v as usize] < best_deg {
                best = v;
                best_deg = deg[v as usize];
            }
        }
        if best != NIL {
            m.s_mate[u] = best;
            m.t_mate[best as usize] = u as i32;
            cnt += 1;
        }
    }
    m.num_edgs += cnt;
    cnt
}

/* ---------- State construction / writeback ---------- */

fn build_state(g: &BipartiteGraph, m: &BipartiteMatching) -> Hkb2State {
    let s_num = g.s_num_vtxs;
    let n_vtxs = g.s_num_vtxs + g.t_num_vtxs;

    let mut v_idx = vec![0usize; n_vtxs + 1];
    for v in 0..g.s_num_vtxs {
        v_idx[v + 1] = v_idx[v] + (g.s_idx[v + 1] - g.s_idx[v]);
    }
    for t in 0..g.t_num_vtxs {
        let v = g.s_num_vtxs + t;
        v_idx[v + 1] = v_idx[v] + (g.t_idx[t + 1] - g.t_idx[t]);
    }
    let mut v_adj = vec![0i32; v_idx[n_vtxs]];
    for v in 0..g.s_num_vtxs {
        let src_begin = g.s_idx[v];
        let src_end = g.s_idx[v + 1];
        let dst_begin = v_idx[v];
        for k in src_begin..src_end {
            v_adj[dst_begin + (k - src_begin)] = g.s_num_vtxs as i32 + g.s_adj[k];
        }
    }
    for t in 0..g.t_num_vtxs {
        let v = g.s_num_vtxs + t;
        let src_begin = g.t_idx[t];
        let src_end = g.t_idx[t + 1];
        let dst_begin = v_idx[v];
        for k in src_begin..src_end {
            v_adj[dst_begin + (k - src_begin)] = g.t_adj[k];
        }
    }

    let mut nodes: Vec<Hkb2Node> = (0..n_vtxs).map(|_| Hkb2Node::new()).collect();
    for s in 0..g.s_num_vtxs {
        if m.s_mate[s] != NIL {
            let t = m.s_mate[s];
            nodes[s].mate = g.s_num_vtxs as i32 + t;
        }
    }
    for t in 0..g.t_num_vtxs {
        if m.t_mate[t] != NIL {
            let s = m.t_mate[t];
            let v = g.s_num_vtxs + t;
            nodes[v].mate = s;
        }
    }

    Hkb2State {
        n_vtxs,
        s_num_vtxs: s_num,
        v_idx,
        v_adj,
        nodes,
        levels: Vec::with_capacity(n_vtxs / 2 + 1),
        br_buckets: Vec::with_capacity(n_vtxs / 2 + 1),
        s_stk_r: Vec::new(),
        s_stk_g: Vec::new(),
        path: Vec::new(),
        todo_num: 0,
        bridge_num: 0,
    }
}

fn writeback_matching(state: &Hkb2State, m: &mut BipartiteMatching) {
    let mut count: i32 = 0;
    for x in m.s_mate.iter_mut() { *x = NIL; }
    for x in m.t_mate.iter_mut() { *x = NIL; }
    for s in 0..state.s_num_vtxs {
        let mt = state.nodes[s].mate;
        if mt != NIL {
            let t = mt - state.s_num_vtxs as i32;
            m.s_mate[s] = t;
            m.t_mate[t as usize] = s as i32;
            count += 1;
        }
    }
    m.num_edgs = count as usize;
}

/* ---------- Top-level HKB2 ---------- */

fn hkb2_iterative_mcm(g: &BipartiteGraph, m: &mut BipartiteMatching) -> i32 {
    let mut state = build_state(g, m);

    let mut num_phases: i32 = 0;
    loop {
        state.reset_phase();
        if !state.run_phase() { break; }
        num_phases += 1;
    }

    writeback_matching(&state, m);
    num_phases
}

/* ---------- Validation ---------- */

fn validate_bipartite_matching(g: &BipartiteGraph, m: &BipartiteMatching) {
    let mut errors = 0;
    let mut matched_s = 0;
    let mut matched_t = 0;

    for u in 0..g.s_num_vtxs {
        if m.s_mate[u] != NIL {
            matched_s += 1;
            let v = m.s_mate[u];
            if v < 0 || (v as usize) >= g.t_num_vtxs {
                eprintln!("ERROR: s_mate[{}] = {} out of range", u, v);
                errors += 1;
            } else if m.t_mate[v as usize] != u as i32 {
                eprintln!("ERROR: s_mate[{}]={} but t_mate[{}]={}", u, v, v, m.t_mate[v as usize]);
                errors += 1;
            } else {
                let st = g.s_idx[u];
                let en = g.s_idx[u + 1];
                if g.s_adj[st..en].binary_search(&v).is_err() {
                    eprintln!("ERROR: edge ({},{}) not in graph", u, v);
                    errors += 1;
                }
            }
        }
    }
    for v in 0..g.t_num_vtxs {
        if m.t_mate[v] != NIL { matched_t += 1; }
    }
    if matched_s != matched_t {
        eprintln!("ERROR: S-matched ({}) != T-matched ({})", matched_s, matched_t);
        errors += 1;
    }

    println!("\n=== Validation Report ===");
    println!("Matching size: {}", m.num_edgs);
    println!("S matched: {}, T matched: {}", matched_s, matched_t);
    println!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
    println!("=========================\n");
}

/* ---------- Graph loader ---------- */

fn load_graph(filename: &str)
    -> Result<(usize, usize, Vec<(usize, usize)>), Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let first = lines.next().ok_or("Empty file")??;
    let parts: Vec<&str> = first.split_whitespace().collect();
    if parts.len() != 3 {
        return Err("First line must have 3 numbers".into());
    }
    let s_num_vtxs: usize = parts[0].parse()?;
    let t_num_vtxs: usize = parts[1].parse()?;
    let m: usize = parts[2].parse()?;

    let mut edges = Vec::with_capacity(m);
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let u: usize = parts[0].parse()?;
            let v: usize = parts[1].parse()?;
            edges.push((u, v));
        }
    }
    Ok((s_num_vtxs, t_num_vtxs, edges))
}

/* ---------- Main ---------- */

fn main() {
    println!("Hopcroft-Karp Bidirectional MV-style (HKB2) - Rust Implementation (CSR)");
    println!("=========================================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename> [--greedy|--greedy-md]", args[0]);
        std::process::exit(1);
    }
    let greedy_mode: i32 = if args.iter().any(|a| a == "--greedy-md") { 2 }
                          else if args.iter().any(|a| a == "--greedy") { 1 }
                          else { 0 };

    match load_graph(&args[1]) {
        Ok((s_num_vtxs, t_num_vtxs, edges)) => {
            println!("Graph: {} s-vertices, {} t-vertices, {} edges",
                     s_num_vtxs, t_num_vtxs, edges.len());

            let g = build_bipartite_graph(s_num_vtxs, t_num_vtxs, &edges);
            let mut m = empty_bipartite_matching(&g);

            let start = Instant::now();

            let greedy_size: usize = match greedy_mode {
                1 => greedy_init(&g, &mut m),
                2 => greedy_init_md(&g, &mut m),
                _ => 0,
            };

            let num_phases = hkb2_iterative_mcm(&g, &mut m);

            let duration = start.elapsed();

            validate_bipartite_matching(&g, &m);

            println!("Phases: {}", num_phases);
            println!("Matching size: {}", m.num_edgs);
            if greedy_mode > 0 {
                println!("Greedy init size: {}", greedy_size);
                if m.num_edgs > 0 {
                    println!("Greedy/Final: {:.2}%", 100.0 * greedy_size as f64 / m.num_edgs as f64);
                } else {
                    println!("Greedy/Final: NA");
                }
            }
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

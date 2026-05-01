/*
 * Hopcroft-Karp Pure Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * CSR adjacency: contiguous flat arrays.
 * Matchbox-style implementation:
 *   - Iterative stack-based DFS (no recursion)
 *   - Edge index array for O(E) per-phase DFS
 *   - Selective cleanup (reset only visited vertices)
 *   - Circular queue containers
 *   - Bidirectional search (from smaller partition)
 *   - Lookahead variant (optional, separate top-level function)
 *   - Greedy and greedy-md initialization
 *   - Comprehensive per-phase statistics (feature = "stats")
 *
 * Refactored to separate input (BipartiteGraph), output (BipartiteMatching),
 * and algorithm state (HKPState).
 *
 * State-side note on direction:
 *   When rvrs == true, HKPState's "s" fields refer to the input graph's
 *   t-vertices (we search starting from t). The output BipartiteMatching
 *   is always in original (s, t) coordinates regardless of rvrs.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;
const INF_LEVEL: i32 = i32::MAX;

/* ---------- Vertex status during BFS/DFS ---------- */

#[repr(u8)]
#[derive(Copy, Clone, PartialEq, Eq)]
enum Stt {
    Idle = 0,
    BfsQueued,
    BfsDone,
    DfsActive,
    DfsDone,
    Last,
}

/* ---------- Helper container types ---------- */

struct CircQueue {
    buf: Vec<i32>,
    cap: usize,
    sz: usize,
    head: usize,
    tail: usize,
}

impl CircQueue {
    fn new() -> Self { CircQueue { buf: Vec::new(), cap: 0, sz: 0, head: 0, tail: 0 } }
    fn init(&mut self, c: usize) {
        self.cap = c; self.sz = 0; self.head = 0; self.tail = 0;
        self.buf = vec![0i32; c];
    }
    fn is_empty(&self) -> bool { self.sz == 0 }
    fn push(&mut self, v: i32) {
        debug_assert!(self.sz < self.cap);
        self.buf[self.tail] = v;
        self.tail = if self.tail + 1 < self.cap { self.tail + 1 } else { 0 };
        self.sz += 1;
    }
    fn front(&self) -> i32 { debug_assert!(self.sz > 0); self.buf[self.head] }
    fn pop(&mut self) {
        debug_assert!(self.sz > 0);
        self.head = if self.head + 1 < self.cap { self.head + 1 } else { 0 };
        self.sz -= 1;
    }
}

struct Stack {
    buf: Vec<i32>,
    sz: usize,
}

impl Stack {
    fn new() -> Self { Stack { buf: Vec::new(), sz: 0 } }
    fn init(&mut self, c: usize) { self.buf = vec![0i32; c]; self.sz = 0; }
    fn is_empty(&self) -> bool { self.sz == 0 }
    fn push(&mut self, v: i32) { self.buf[self.sz] = v; self.sz += 1; }
    fn top(&self) -> i32 { self.buf[self.sz - 1] }
    fn pop(&mut self) { self.sz -= 1; }
}

struct IdxQueue {
    nxt: Vec<i32>,
    prv: Vec<i32>,
    head: i32,
    sz: usize,
    #[allow(dead_code)]
    cap: usize,
}

impl IdxQueue {
    fn new() -> Self { IdxQueue { nxt: Vec::new(), prv: Vec::new(), head: NIL, sz: 0, cap: 0 } }
    fn init(&mut self, c: usize) {
        self.cap = c; self.sz = 0; self.head = NIL;
        self.nxt = vec![NIL; c];
        self.prv = vec![NIL; c];
    }
    fn push(&mut self, v: i32) {
        let vu = v as usize;
        self.nxt[vu] = self.head;
        self.prv[vu] = NIL;
        if self.head != NIL { self.prv[self.head as usize] = v; }
        self.head = v;
        self.sz += 1;
    }
    fn erase(&mut self, v: i32) {
        let vu = v as usize;
        if self.prv[vu] != NIL { self.nxt[self.prv[vu] as usize] = self.nxt[vu]; }
        else { self.head = self.nxt[vu]; }
        if self.nxt[vu] != NIL { self.prv[self.nxt[vu] as usize] = self.prv[vu]; }
        self.nxt[vu] = NIL;
        self.prv[vu] = NIL;
        self.sz -= 1;
    }
    fn first(&self) -> i32 { self.head }
    fn next(&self, v: i32) -> i32 { self.nxt[v as usize] }
}

/* ---------- Statistics (compiled only with --features stats) ---------- */

#[cfg(feature = "stats")]
#[derive(Clone, Copy)]
struct PhaseStats {
    bfs_vtx: i64,
    bfs_edg: i64,
    dfs_vtx: i64,
    dfs_edg: i64,
    num_augmentations: i32,
    shortest_path_len: i32,
    agg_aug_path_len: i64,
    min_aug_path_len: i32,
    max_aug_path_len: i32,
}

#[cfg(feature = "stats")]
impl PhaseStats {
    fn new() -> Self {
        PhaseStats {
            bfs_vtx: 0, bfs_edg: 0, dfs_vtx: 0, dfs_edg: 0,
            num_augmentations: 0,
            shortest_path_len: 0,
            agg_aug_path_len: 0,
            min_aug_path_len: INF_LEVEL,
            max_aug_path_len: 0,
        }
    }
}

#[cfg(feature = "stats")]
struct AlgoStats {
    phases: Vec<PhaseStats>,
    greedy_card: i32,
    reversed: bool,
}

#[cfg(feature = "stats")]
impl AlgoStats {
    fn new() -> Self { AlgoStats { phases: Vec::new(), greedy_card: 0, reversed: false } }
    fn reset(&mut self) { self.phases.clear(); self.greedy_card = 0; self.reversed = false; }
    fn print(&self) {
        let mut tb_v = 0i64; let mut tb_e = 0i64;
        let mut td_v = 0i64; let mut td_e = 0i64;
        let mut t_aug = 0i32; let mut t_apl = 0i64;
        for p in &self.phases {
            tb_v += p.bfs_vtx; tb_e += p.bfs_edg;
            td_v += p.dfs_vtx; td_e += p.dfs_edg;
            t_aug += p.num_augmentations;
            t_apl += p.agg_aug_path_len;
        }
        println!("\n=== Algorithm Statistics ===");
        println!("Search direction: {}", if self.reversed { "T->S (reversed)" } else { "S->T (normal)" });
        println!("Greedy init card: {}", self.greedy_card);
        println!("Phases: {}", self.phases.len());
        println!("Total BFS vertices: {}, edges: {}", tb_v, tb_e);
        println!("Total DFS vertices: {}, edges: {}", td_v, td_e);
        println!("Total augmentations: {}", t_aug);
        println!("Total aug path length: {}", t_apl);
        if t_aug > 0 {
            println!("Avg aug path length: {:.2}", t_apl as f64 / t_aug as f64);
        }
        println!("\nPer-phase breakdown:");
        for (i, p) in self.phases.iter().enumerate() {
            let avg = if p.num_augmentations > 0 {
                p.agg_aug_path_len as f64 / p.num_augmentations as f64
            } else { 0.0 };
            println!("  Phase {}: shortest={} augs={} bfs_v={} bfs_e={} dfs_v={} dfs_e={} avg_len={:.1}",
                     i + 1, p.shortest_path_len, p.num_augmentations,
                     p.bfs_vtx, p.bfs_edg, p.dfs_vtx, p.dfs_edg, avg);
        }
        println!("===========================\n");
    }
}

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

/* ---------- State: HKPState ---------- */

struct HKPState {
    rvrs: bool,
    s_count: usize,
    #[allow(dead_code)]
    t_count: usize,

    s_level: Vec<i32>,
    t_level: Vec<i32>,
    s_stt: Vec<Stt>,
    t_stt: Vec<Stt>,
    s_ptr: Vec<i32>,
    s_edge_idx: Vec<usize>,
    s_lkhd_idx: Vec<usize>,

    s_bfs_queue: CircQueue,
    s_done_queue: CircQueue,
    t_done_queue: CircQueue,
    s_dfs_stack: Stack,
    s_last_queue: CircQueue,
    t_last_queue: CircQueue,
    s_exposed: IdxQueue,

    #[cfg(feature = "stats")]
    stats: AlgoStats,
}

/* ---------- BipartiteGraph construction ---------- */

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
    for u in 0..s_num_vtxs {
        s_adj.extend_from_slice(&s_tmp[u]);
    }

    let mut t_idx = vec![0usize; t_num_vtxs + 1];
    for v in 0..t_num_vtxs { t_idx[v + 1] = t_idx[v] + t_tmp[v].len(); }
    let mut t_adj: Vec<i32> = Vec::with_capacity(t_idx[t_num_vtxs]);
    for v in 0..t_num_vtxs {
        t_adj.extend_from_slice(&t_tmp[v]);
    }

    let num_edgs = s_idx[s_num_vtxs];

    BipartiteGraph {
        s_num_vtxs,
        t_num_vtxs,
        num_edgs,
        s_idx,
        t_idx,
        s_adj,
        t_adj,
    }
}

/* ---------- BipartiteMatching construction ---------- */

fn empty_bipartite_matching(g: &BipartiteGraph) -> BipartiteMatching {
    BipartiteMatching {
        s_num_vtxs: g.s_num_vtxs,
        t_num_vtxs: g.t_num_vtxs,
        num_edgs: 0,
        s_mate: vec![NIL; g.s_num_vtxs],
        t_mate: vec![NIL; g.t_num_vtxs],
    }
}

/* ---------- HKPState construction ---------- */

fn init_hkp_state(g: &BipartiteGraph, rvrs: bool) -> HKPState {
    let s_count = if rvrs { g.t_num_vtxs } else { g.s_num_vtxs };
    let t_count = if rvrs { g.s_num_vtxs } else { g.t_num_vtxs };

    let mut s = HKPState {
        rvrs,
        s_count,
        t_count,
        s_level: vec![INF_LEVEL; s_count],
        t_level: vec![INF_LEVEL; t_count],
        s_stt: vec![Stt::Idle; s_count],
        t_stt: vec![Stt::Idle; t_count],
        s_ptr: vec![NIL; s_count],
        s_edge_idx: vec![0usize; s_count],
        s_lkhd_idx: vec![0usize; s_count],
        s_bfs_queue: CircQueue::new(),
        s_done_queue: CircQueue::new(),
        t_done_queue: CircQueue::new(),
        s_dfs_stack: Stack::new(),
        s_last_queue: CircQueue::new(),
        t_last_queue: CircQueue::new(),
        s_exposed: IdxQueue::new(),
        #[cfg(feature = "stats")]
        stats: AlgoStats::new(),
    };

    let sc = s_count.max(1);
    let tc = t_count.max(1);
    s.s_bfs_queue.init(sc);
    s.s_done_queue.init(sc);
    s.t_done_queue.init(tc);
    s.s_dfs_stack.init(sc);
    let mp = s_count.min(t_count).max(1);
    s.s_last_queue.init(mp);
    s.t_last_queue.init(mp);
    s.s_exposed.init(s_count);

    #[cfg(feature = "stats")]
    {
        s.stats.reset();
        s.stats.reversed = rvrs;
    }

    s
}

/* ---------- CSR adjacency / mate accessors honoring rvrs ---------- */

#[inline]
fn src_idx_begin(g: &BipartiteGraph, s: &HKPState, u: i32) -> usize {
    if s.rvrs { g.t_idx[u as usize] } else { g.s_idx[u as usize] }
}
#[inline]
fn src_idx_end(g: &BipartiteGraph, s: &HKPState, u: i32) -> usize {
    if s.rvrs { g.t_idx[u as usize + 1] } else { g.s_idx[u as usize + 1] }
}
#[inline]
fn src_adj_at<'a>(g: &'a BipartiteGraph, s: &HKPState, j: usize) -> i32 {
    if s.rvrs { g.t_adj[j] } else { g.s_adj[j] }
}

#[inline]
fn src_mate(m: &BipartiteMatching, s: &HKPState, u: i32) -> i32 {
    if s.rvrs { m.t_mate[u as usize] } else { m.s_mate[u as usize] }
}

#[inline]
fn tgt_mate(m: &BipartiteMatching, s: &HKPState, v: i32) -> i32 {
    if s.rvrs { m.s_mate[v as usize] } else { m.t_mate[v as usize] }
}

#[inline]
fn set_mate(m: &mut BipartiteMatching, s: &HKPState, u: i32, v: i32) {
    if s.rvrs {
        m.t_mate[u as usize] = v;
        m.s_mate[v as usize] = u;
    } else {
        m.s_mate[u as usize] = v;
        m.t_mate[v as usize] = u;
    }
}

/* ---------- Greedy initial matching: simple ---------- */

fn greedy_init(g: &BipartiteGraph, m: &mut BipartiteMatching) -> usize {
    let mut cnt: usize = 0;
    for u in 0..g.s_num_vtxs {
        if m.s_mate[u] != NIL { continue; }
        let st = g.s_idx[u];
        let en = g.s_idx[u + 1];
        for j in st..en {
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

/* ---------- Greedy initial matching: min-degree ---------- */

fn greedy_init_md(g: &BipartiteGraph, m: &mut BipartiteMatching) -> usize {
    let mut cnt: usize = 0;
    let mut deg = vec![0usize; g.t_num_vtxs];
    for u in 0..g.s_num_vtxs {
        let st = g.s_idx[u];
        let en = g.s_idx[u + 1];
        for j in st..en { deg[g.s_adj[j] as usize] += 1; }
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
        let st = g.s_idx[u];
        let en = g.s_idx[u + 1];
        for j in st..en {
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

/* ---------- BFS: compute shortest augmenting path length ---------- */

fn bfs_shortest_level(
    g: &BipartiteGraph,
    m: &BipartiteMatching,
    s: &mut HKPState,
    #[cfg(feature = "stats")] ps: &mut PhaseStats,
) -> i32 {
    let mut shortest = INF_LEVEL;

    let mut sf = s.s_exposed.first();
    while sf != NIL {
        #[cfg(feature = "stats")] { ps.bfs_vtx += 1; }
        let sfu = sf as usize;
        s.s_level[sfu] = 0;
        s.s_bfs_queue.push(sf);
        s.s_stt[sfu] = Stt::BfsQueued;
        sf = s.s_exposed.next(sf);
    }

    let mut cur_layer = 0;
    while !s.s_bfs_queue.is_empty() {
        let src = s.s_bfs_queue.front();
        s.s_bfs_queue.pop();
        s.s_done_queue.push(src);
        let su = src as usize;
        s.s_stt[su] = Stt::BfsDone;

        if cur_layer < s.s_level[su] / 2 {
            if shortest != INF_LEVEL { break; }
            cur_layer += 1;
        }

        let a_beg = src_idx_begin(g, s, src);
        let a_end = src_idx_end(g, s, src);
        for j in a_beg..a_end {
            let t = src_adj_at(g, s, j);
            #[cfg(feature = "stats")] { ps.bfs_edg += 1; }
            let tu = t as usize;
            if s.t_stt[tu] == Stt::BfsDone { continue; }
            #[cfg(feature = "stats")] { ps.bfs_vtx += 1; }
            s.t_level[tu] = s.s_level[su] + 1;
            s.t_done_queue.push(t);
            s.t_stt[tu] = Stt::BfsDone;

            let ss = tgt_mate(m, s, t);
            if ss == NIL {
                shortest = s.t_level[tu];
            } else {
                #[cfg(feature = "stats")] { ps.bfs_edg += 1; ps.bfs_vtx += 1; }
                let ssu = ss as usize;
                s.s_level[ssu] = s.t_level[tu] + 1;
                s.s_bfs_queue.push(ss);
                s.s_stt[ssu] = Stt::BfsQueued;
            }
        }
    }
    shortest
}

/* ---------- DFS: HK-style ---------- */

fn dfs_find_paths_hk(
    g: &BipartiteGraph,
    m: &BipartiteMatching,
    s: &mut HKPState,
    shortest: i32,
    #[cfg(feature = "stats")] ps: &mut PhaseStats,
) {
    let mut sf = s.s_exposed.first();
    while sf != NIL {
        let sf_next = s.s_exposed.next(sf);
        #[cfg(feature = "stats")] { ps.dfs_vtx += 1; }
        s.s_dfs_stack.push(sf);
        s.s_stt[sf as usize] = Stt::DfsActive;

        while !s.s_dfs_stack.is_empty() {
            let src = s.s_dfs_stack.top();
            let su = src as usize;
            let a_beg = src_idx_begin(g, s, src);
            let a_end = src_idx_end(g, s, src);
            let num_edg = a_end - a_beg;

            while s.s_edge_idx[su] < num_edg {
                let t = src_adj_at(g, s, a_beg + s.s_edge_idx[su]);
                #[cfg(feature = "stats")] { ps.dfs_edg += 1; }
                let tu = t as usize;
                let ss = tgt_mate(m, s, t);
                let skip = (s.t_level[tu] != s.s_level[su] + 1) ||
                           ((ss == NIL) && (s.t_stt[tu] == Stt::Last)) ||
                           ((ss != NIL) &&
                              ((s.s_level[su] + 1 == shortest) ||
                               (s.s_stt[ss as usize] == Stt::DfsActive) ||
                               (s.s_stt[ss as usize] == Stt::DfsDone)));
                if skip { s.s_edge_idx[su] += 1; } else { break; }
            }

            if s.s_edge_idx[su] < num_edg {
                let t = src_adj_at(g, s, a_beg + s.s_edge_idx[su]);
                let tu = t as usize;
                let ss = tgt_mate(m, s, t);
                #[cfg(feature = "stats")] { ps.dfs_vtx += 1; }
                if ss == NIL {
                    s.s_last_queue.push(src);
                    s.t_last_queue.push(t);
                    s.t_stt[tu] = Stt::Last;
                    while !s.s_dfs_stack.is_empty() {
                        let top = s.s_dfs_stack.top();
                        s.s_stt[top as usize] = Stt::DfsDone;
                        s.s_dfs_stack.pop();
                    }
                    break;
                }
                #[cfg(feature = "stats")] { ps.dfs_edg += 1; ps.dfs_vtx += 1; }
                let ssu = ss as usize;
                s.s_ptr[ssu] = src;
                s.s_dfs_stack.push(ss);
                s.s_stt[ssu] = Stt::DfsActive;
                s.s_edge_idx[su] += 1;
            } else {
                s.s_dfs_stack.pop();
                s.s_stt[su] = Stt::DfsDone;
            }
        }
        sf = sf_next;
    }
}

/* ---------- DFS: lookahead-style ---------- */

fn dfs_find_paths_lkhd(
    g: &BipartiteGraph,
    m: &BipartiteMatching,
    s: &mut HKPState,
    #[cfg(feature = "stats")] ps: &mut PhaseStats,
) {
    let mut sf = s.s_exposed.first();
    while sf != NIL {
        let sf_next = s.s_exposed.next(sf);
        #[cfg(feature = "stats")] { ps.dfs_vtx += 1; }
        s.s_dfs_stack.push(sf);
        s.s_stt[sf as usize] = Stt::DfsActive;

        while !s.s_dfs_stack.is_empty() {
            let src = s.s_dfs_stack.top();
            let su = src as usize;
            let a_beg = src_idx_begin(g, s, src);
            let a_end = src_idx_end(g, s, src);
            let num_edg = a_end - a_beg;

            /* Lookahead: check for direct free t-neighbor */
            while s.s_lkhd_idx[su] < num_edg {
                let t = src_adj_at(g, s, a_beg + s.s_lkhd_idx[su]);
                #[cfg(feature = "stats")] { ps.bfs_edg += 1; }
                let tu = t as usize;
                let ss = tgt_mate(m, s, t);
                let skip = (ss == NIL && s.t_stt[tu] != Stt::Idle) || ss != NIL;
                if skip { s.s_lkhd_idx[su] += 1; } else { break; }
            }

            if s.s_lkhd_idx[su] < num_edg {
                let t = src_adj_at(g, s, a_beg + s.s_lkhd_idx[su]);
                let tu = t as usize;
                #[cfg(feature = "stats")] { ps.bfs_vtx += 1; }
                s.s_last_queue.push(src);
                s.t_last_queue.push(t);
                s.t_stt[tu] = Stt::Last;
                while !s.s_dfs_stack.is_empty() {
                    let sss = s.s_dfs_stack.top();
                    s.s_dfs_stack.pop();
                    s.s_done_queue.push(sss);
                    s.s_stt[sss as usize] = Stt::DfsDone;
                }
                s.s_lkhd_idx[su] += 1;
                break;
            }

            /* Regular DFS deeper */
            while s.s_edge_idx[su] < num_edg {
                let t = src_adj_at(g, s, a_beg + s.s_edge_idx[su]);
                #[cfg(feature = "stats")] { ps.dfs_edg += 1; }
                let ss = tgt_mate(m, s, t);
                let skip = ss == NIL || (ss != NIL && s.s_stt[ss as usize] != Stt::Idle);
                if skip { s.s_edge_idx[su] += 1; } else { break; }
            }

            if s.s_edge_idx[su] < num_edg {
                let t = src_adj_at(g, s, a_beg + s.s_edge_idx[su]);
                let ss = tgt_mate(m, s, t);
                #[cfg(feature = "stats")] { ps.dfs_vtx += 1; ps.dfs_edg += 1; ps.dfs_vtx += 1; }
                let ssu = ss as usize;
                s.s_ptr[ssu] = src;
                s.s_dfs_stack.push(ss);
                s.s_stt[ssu] = Stt::DfsActive;
                s.s_edge_idx[su] += 1;
            } else {
                s.s_dfs_stack.pop();
                s.s_done_queue.push(src);
                s.s_stt[su] = Stt::DfsDone;
            }
        }
        sf = sf_next;
    }
}

/* ---------- Augment along path from sLast to root via sPtr ---------- */

fn augment(m: &mut BipartiteMatching, s: &mut HKPState, s_last: i32, t_last: i32) -> i32 {
    let mut src = s_last;
    let mut t = t_last;
    let mut k = 0i32;
    while src != NIL {
        let su = src as usize;
        let tt = src_mate(m, s, src);
        set_mate(m, s, src, t);
        if s.s_ptr[su] == NIL { s.s_exposed.erase(src); }
        src = s.s_ptr[su];
        t = tt;
        k += 1;
    }
    2 * k - 1
}

/* ---------- Cleanup after HK phase ---------- */

fn cleanup_hk(s: &mut HKPState) {
    while !s.s_bfs_queue.is_empty() {
        let u = s.s_bfs_queue.front(); s.s_bfs_queue.pop();
        let uu = u as usize;
        s.s_ptr[uu] = NIL; s.s_level[uu] = INF_LEVEL;
        s.s_stt[uu] = Stt::Idle;
        s.s_edge_idx[uu] = 0; s.s_lkhd_idx[uu] = 0;
    }
    while !s.s_done_queue.is_empty() {
        let u = s.s_done_queue.front(); s.s_done_queue.pop();
        let uu = u as usize;
        s.s_ptr[uu] = NIL; s.s_level[uu] = INF_LEVEL;
        s.s_stt[uu] = Stt::Idle;
        s.s_edge_idx[uu] = 0; s.s_lkhd_idx[uu] = 0;
    }
    while !s.t_done_queue.is_empty() {
        let t = s.t_done_queue.front(); s.t_done_queue.pop();
        let tu = t as usize;
        s.t_level[tu] = INF_LEVEL;
        s.t_stt[tu] = Stt::Idle;
    }
}

/* ---------- Top-level Hopcroft-Karp Pure (HK mode) ---------- */

fn hopcroft_karp_pure(g: &BipartiteGraph, m: &mut BipartiteMatching, rvrs: bool) {
    let mut s = init_hkp_state(g, rvrs);

    #[cfg(feature = "stats")]
    {
        s.stats.greedy_card = m.num_edgs as i32;
    }

    /* Push exposed src vertices into the IdxQueue */
    for u in 0..s.s_count {
        if src_mate(m, &s, u as i32) == NIL { s.s_exposed.push(u as i32); }
    }

    let mut phase_count: i32 = 0;
    let mut new_edgs: usize = 0;
    loop {
        #[cfg(feature = "stats")]
        let mut ps = PhaseStats::new();

        let shortest = bfs_shortest_level(
            g, m, &mut s,
            #[cfg(feature = "stats")] &mut ps,
        );

        if shortest == INF_LEVEL {
            cleanup_hk(&mut s);
            break;
        }

        phase_count += 1;
        #[cfg(feature = "stats")] { ps.shortest_path_len = shortest; }

        dfs_find_paths_hk(
            g, m, &mut s, shortest,
            #[cfg(feature = "stats")] &mut ps,
        );

        while !s.t_last_queue.is_empty() {
            let sl = s.s_last_queue.front(); s.s_last_queue.pop();
            let tl = s.t_last_queue.front(); s.t_last_queue.pop();
            let _path_len = augment(m, &mut s, sl, tl);
            new_edgs += 1;
            #[cfg(feature = "stats")] {
                ps.num_augmentations += 1;
                ps.agg_aug_path_len += _path_len as i64;
                if _path_len < ps.min_aug_path_len { ps.min_aug_path_len = _path_len; }
                if _path_len > ps.max_aug_path_len { ps.max_aug_path_len = _path_len; }
            }
        }

        cleanup_hk(&mut s);

        #[cfg(feature = "stats")] { s.stats.phases.push(ps); }
    }

    m.num_edgs += new_edgs;

    #[cfg(feature = "stats")] { s.stats.print(); }
    println!("Phases: {}", phase_count);
}

/* ---------- Top-level Hopcroft-Karp Pure (lookahead mode) ---------- */

fn hopcroft_karp_pure_lkhd(g: &BipartiteGraph, m: &mut BipartiteMatching, rvrs: bool) {
    let mut s = init_hkp_state(g, rvrs);

    #[cfg(feature = "stats")]
    {
        s.stats.greedy_card = m.num_edgs as i32;
    }

    for u in 0..s.s_count {
        if src_mate(m, &s, u as i32) == NIL { s.s_exposed.push(u as i32); }
    }

    let mut pass_count: i32 = 0;
    let mut new_edgs: usize = 0;
    loop {
        #[cfg(feature = "stats")]
        let mut ps = PhaseStats::new();

        dfs_find_paths_lkhd(
            g, m, &mut s,
            #[cfg(feature = "stats")] &mut ps,
        );

        if s.t_last_queue.is_empty() { break; }

        pass_count += 1;
        while !s.t_last_queue.is_empty() {
            let sl = s.s_last_queue.front(); s.s_last_queue.pop();
            let tl = s.t_last_queue.front(); s.t_last_queue.pop();
            s.t_stt[tl as usize] = Stt::Idle;
            let _path_len = augment(m, &mut s, sl, tl);
            new_edgs += 1;
            #[cfg(feature = "stats")] {
                ps.num_augmentations += 1;
                ps.agg_aug_path_len += _path_len as i64;
            }
        }

        /* Cleanup: reset only visited vertices */
        while !s.s_done_queue.is_empty() {
            let u = s.s_done_queue.front(); s.s_done_queue.pop();
            let uu = u as usize;
            s.s_ptr[uu] = NIL;
            s.s_stt[uu] = Stt::Idle;
            s.s_edge_idx[uu] = 0; s.s_lkhd_idx[uu] = 0;
        }

        #[cfg(feature = "stats")] { s.stats.phases.push(ps); }
    }

    m.num_edgs += new_edgs;

    #[cfg(feature = "stats")] { s.stats.print(); }
    println!("Phases: {}", pass_count);
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
                if !g.s_adj[st..en].binary_search(&v).is_ok() {
                    eprintln!("ERROR: edge ({},{}) not in graph", u, v);
                    errors += 1;
                }
            }
        }
    }
    for v in 0..g.t_num_vtxs {
        if m.t_mate[v] != NIL { matched_t += 1; }
    }

    println!("\n=== Validation Report ===");
    println!("Matching size: {}", m.num_edgs);
    println!("S matched: {}, T matched: {}", matched_s, matched_t);
    println!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
    println!("=========================\n");
}

/* ---------- Graph loader ---------- */

fn load_graph(filename: &str) -> Result<(usize, usize, Vec<(usize, usize)>), Box<dyn std::error::Error>> {
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
    println!("Hopcroft-Karp Pure Algorithm - Rust Implementation (CSR)");
    println!("==========================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename> [--greedy|--greedy-md] [--lkhd] [--rvrs|--fwd]", args[0]);
        std::process::exit(1);
    }

    let greedy_mode: i32 = if args.iter().any(|a| a == "--greedy-md") { 2 }
                          else if args.iter().any(|a| a == "--greedy") { 1 }
                          else { 0 };
    let use_lkhd = args.iter().any(|a| a == "--lkhd");
    let mut force_rvrs = false;
    let mut auto_dir = true;
    for a in &args[2..] {
        if a == "--rvrs" { force_rvrs = true; auto_dir = false; }
        else if a == "--fwd" { force_rvrs = false; auto_dir = false; }
    }

    match load_graph(&args[1]) {
        Ok((s_num_vtxs, t_num_vtxs, edges)) => {
            println!("Graph: {} s-vertices, {} t-vertices, {} edges",
                     s_num_vtxs, t_num_vtxs, edges.len());

            let use_rvrs = if auto_dir { t_num_vtxs < s_num_vtxs } else { force_rvrs };

            println!("Mode: {}", if use_lkhd { "lookahead DFS" } else { "HK (BFS+DFS)" });
            println!("Direction: {}{}",
                     if use_rvrs { "T->S (reversed)" } else { "S->T (normal)" },
                     if auto_dir { " (auto)" } else { " (forced)" });

            let bipartite_graph = build_bipartite_graph(s_num_vtxs, t_num_vtxs, &edges);
            let mut bipartite_matching = empty_bipartite_matching(&bipartite_graph);

            let start = Instant::now();

            let greedy_size: usize = match greedy_mode {
                1 => greedy_init(&bipartite_graph, &mut bipartite_matching),
                2 => greedy_init_md(&bipartite_graph, &mut bipartite_matching),
                _ => 0,
            };

            if use_lkhd {
                hopcroft_karp_pure_lkhd(&bipartite_graph, &mut bipartite_matching, use_rvrs);
            } else {
                hopcroft_karp_pure(&bipartite_graph, &mut bipartite_matching, use_rvrs);
            }

            let duration = start.elapsed();

            validate_bipartite_matching(&bipartite_graph, &bipartite_matching);

            println!("Matching size: {}", bipartite_matching.num_edgs);
            if greedy_mode > 0 {
                let fs = bipartite_matching.num_edgs;
                println!("Greedy init size: {}", greedy_size);
                if fs > 0 {
                    println!("Greedy/Final: {:.2}%", 100.0 * greedy_size as f64 / fs as f64);
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

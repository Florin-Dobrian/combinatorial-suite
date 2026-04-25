/*
 * Hopcroft-Karp Pure Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * Incorporates all optimizations from the Matchbox library:
 * - Iterative stack-based DFS (no recursion)
 * - Edge index array (s_idx) for O(E) per-phase DFS
 * - Selective cleanup (reset only visited vertices)
 * - Circular queue containers (no overflow on reuse)
 * - Bidirectional search (from smaller partition)
 * - Lookahead variant (eMplPathDfsLkhd style, optional, via --lkhd)
 * - Greedy and greedy-md initialization
 * - Comprehensive per-phase statistics (feature = "stats")
 *
 * Input format: "L R M" header then M lines of "u v" edges (0-indexed).
 *
 * Rust port. All integers, no hash containers, fully deterministic.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

/* =========================================================================
 * Constants
 * ========================================================================= */
const NIL: i32 = -1;
const INF_LEVEL: i32 = i32::MAX;

/* Vertex status during BFS/DFS */
#[repr(u8)]
#[derive(Copy, Clone, PartialEq, Eq)]
enum Stt {
    Idle = 0,
    BfsQueued,   /* in BFS queue, not yet processed */
    BfsDone,     /* BFS processed */
    DfsActive,   /* on DFS stack */
    DfsDone,     /* DFS finished (found path or dead end) */
    Last,        /* T-vertex: endpoint of found augmenting path */
}

/* =========================================================================
 * Circular Queue (matches Matchbox VecItmQue)
 *
 * Fixed-capacity circular buffer. Capacity = max simultaneous items.
 * Push/Pop wrap around, so total pushes can exceed capacity safely.
 * ========================================================================= */
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
    #[allow(dead_code)]
    fn clear(&mut self) { self.sz = 0; self.head = 0; self.tail = 0; }
}

/* =========================================================================
 * Stack (matches Matchbox VecItmStk)
 * ========================================================================= */
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
    #[allow(dead_code)]
    fn clear(&mut self) { self.sz = 0; }
}

/* =========================================================================
 * Indexed Queue (matches Matchbox LstItmIdxdQue but array-based)
 *
 * Supports O(1) Push, Erase, First, Next for the exposed-vertex set.
 * Doubly-linked list embedded in arrays.
 * ========================================================================= */
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
    #[allow(dead_code)]
    fn is_empty(&self) -> bool { self.sz == 0 }
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

/* =========================================================================
 * Statistics (compiled only with --features stats)
 * ========================================================================= */
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

/* =========================================================================
 * HopcroftKarpPure
 * ========================================================================= */
struct HopcroftKarpPure {
    s_count: usize,
    t_count: usize,
    src_adj: Vec<Vec<usize>>,

    s_mate: Vec<i32>,
    t_mate: Vec<i32>,
    s_level: Vec<i32>,
    t_level: Vec<i32>,
    s_stt: Vec<Stt>,
    t_stt: Vec<Stt>,
    s_ptr: Vec<i32>,
    s_idx: Vec<usize>,
    s_lkhd: Vec<usize>,

    s_bfs_queue: CircQueue,
    s_done_queue: CircQueue,
    t_done_queue: CircQueue,
    s_dfs_stack: Stack,
    s_last_queue: CircQueue,
    t_last_queue: CircQueue,
    s_exposed: IdxQueue,

    greedy_size: usize,
    rvrs: bool,

    #[cfg(feature = "stats")]
    stats: AlgoStats,
}

impl HopcroftKarpPure {
    fn new() -> Self {
        HopcroftKarpPure {
            s_count: 0, t_count: 0,
            src_adj: Vec::new(),
            s_mate: Vec::new(), t_mate: Vec::new(),
            s_level: Vec::new(), t_level: Vec::new(),
            s_stt: Vec::new(), t_stt: Vec::new(),
            s_ptr: Vec::new(),
            s_idx: Vec::new(),
            s_lkhd: Vec::new(),
            s_bfs_queue: CircQueue::new(),
            s_done_queue: CircQueue::new(),
            t_done_queue: CircQueue::new(),
            s_dfs_stack: Stack::new(),
            s_last_queue: CircQueue::new(),
            t_last_queue: CircQueue::new(),
            s_exposed: IdxQueue::new(),
            greedy_size: 0,
            rvrs: false,
            #[cfg(feature = "stats")]
            stats: AlgoStats::new(),
        }
    }

    fn build(&mut self, left_count: usize, right_count: usize,
             edges: &[(usize, usize)], force_rvrs: bool) {
        self.rvrs = force_rvrs;
        if self.rvrs { self.s_count = right_count; self.t_count = left_count; }
        else         { self.s_count = left_count;  self.t_count = right_count; }

        self.src_adj = vec![Vec::<usize>::new(); self.s_count];
        for &(u, v) in edges {
            if self.rvrs {
                if v < self.s_count && u < self.t_count { self.src_adj[v].push(u); }
            } else {
                if u < self.s_count && v < self.t_count { self.src_adj[u].push(v); }
            }
        }
        for adj in &mut self.src_adj { adj.sort_unstable(); adj.dedup(); }

        self.s_mate  = vec![NIL; self.s_count];
        self.t_mate  = vec![NIL; self.t_count];
        self.s_level = vec![INF_LEVEL; self.s_count];
        self.t_level = vec![INF_LEVEL; self.t_count];
        self.s_stt   = vec![Stt::Idle; self.s_count];
        self.t_stt   = vec![Stt::Idle; self.t_count];
        self.s_ptr   = vec![NIL; self.s_count];
        self.s_idx   = vec![0usize; self.s_count];
        self.s_lkhd  = vec![0usize; self.s_count];

        let sc = self.s_count.max(1);
        let tc = self.t_count.max(1);
        self.s_bfs_queue.init(sc);
        self.s_done_queue.init(sc);
        self.t_done_queue.init(tc);
        self.s_dfs_stack.init(sc);
        let mp = self.s_count.min(self.t_count).max(1);
        self.s_last_queue.init(mp);
        self.t_last_queue.init(mp);
        self.s_exposed.init(self.s_count);
        self.greedy_size = 0;
    }

    /* Greedy: simple sequential */
    fn greedy_init(&mut self) -> usize {
        let mut cnt = 0usize;
        for s in 0..self.s_count {
            if self.s_mate[s] != NIL { continue; }
            let n = self.src_adj[s].len();
            for i in 0..n {
                let t = self.src_adj[s][i];
                if self.t_mate[t] == NIL {
                    self.s_mate[s] = t as i32;
                    self.t_mate[t] = s as i32;
                    cnt += 1;
                    break;
                }
            }
        }
        cnt
    }

    /* Greedy: min-degree heuristic */
    fn greedy_init_md(&mut self) -> usize {
        let mut cnt = 0usize;
        let mut t_deg = vec![0usize; self.t_count];
        for s in 0..self.s_count {
            for &t in &self.src_adj[s] { t_deg[t] += 1; }
        }
        let mut order: Vec<usize> = (0..self.s_count).collect();
        order.sort_unstable_by(|&a, &b| {
            self.src_adj[a].len().cmp(&self.src_adj[b].len()).then(a.cmp(&b))
        });
        for s in order {
            if self.s_mate[s] != NIL { continue; }
            let mut best: i32 = NIL;
            let mut best_deg = usize::MAX;
            let n = self.src_adj[s].len();
            for i in 0..n {
                let t = self.src_adj[s][i];
                if self.t_mate[t] == NIL && t_deg[t] < best_deg {
                    best = t as i32;
                    best_deg = t_deg[t];
                }
            }
            if best != NIL {
                self.s_mate[s] = best;
                self.t_mate[best as usize] = s as i32;
                cnt += 1;
            }
        }
        cnt
    }

    /* BFS: compute shortest augmenting path length */
    fn bfs_shortest_level(
        &mut self,
        #[cfg(feature = "stats")] ps: &mut PhaseStats,
    ) -> i32 {
        let mut shortest = INF_LEVEL;

        let mut sf = self.s_exposed.first();
        while sf != NIL {
            #[cfg(feature = "stats")] { ps.bfs_vtx += 1; }
            let sfu = sf as usize;
            self.s_level[sfu] = 0;
            self.s_bfs_queue.push(sf);
            self.s_stt[sfu] = Stt::BfsQueued;
            sf = self.s_exposed.next(sf);
        }

        let mut cur_layer = 0;
        while !self.s_bfs_queue.is_empty() {
            let s = self.s_bfs_queue.front();
            self.s_bfs_queue.pop();
            self.s_done_queue.push(s);
            let su = s as usize;
            self.s_stt[su] = Stt::BfsDone;

            if cur_layer < self.s_level[su] / 2 {
                if shortest != INF_LEVEL { break; }
                cur_layer += 1;
            }

            let num_edg = self.src_adj[su].len();
            for i in 0..num_edg {
                let t = self.src_adj[su][i];
                #[cfg(feature = "stats")] { ps.bfs_edg += 1; }
                if self.t_stt[t] == Stt::BfsDone { continue; }
                #[cfg(feature = "stats")] { ps.bfs_vtx += 1; }
                self.t_level[t] = self.s_level[su] + 1;
                self.t_done_queue.push(t as i32);
                self.t_stt[t] = Stt::BfsDone;

                let ss = self.t_mate[t];
                if ss == NIL {
                    shortest = self.t_level[t];
                } else {
                    #[cfg(feature = "stats")] { ps.bfs_edg += 1; ps.bfs_vtx += 1; }
                    let ssu = ss as usize;
                    self.s_level[ssu] = self.t_level[t] + 1;
                    self.s_bfs_queue.push(ss);
                    self.s_stt[ssu] = Stt::BfsQueued;
                }
            }
        }
        shortest
    }

    /* DFS: find maximal set of shortest augmenting paths (HK) */
    fn dfs_find_paths_hk(
        &mut self,
        shortest: i32,
        #[cfg(feature = "stats")] ps: &mut PhaseStats,
    ) {
        let mut sf = self.s_exposed.first();
        while sf != NIL {
            let sf_next = self.s_exposed.next(sf);
            #[cfg(feature = "stats")] { ps.dfs_vtx += 1; }
            self.s_dfs_stack.push(sf);
            self.s_stt[sf as usize] = Stt::DfsActive;

            while !self.s_dfs_stack.is_empty() {
                let s = self.s_dfs_stack.top();
                let su = s as usize;
                let num_edg = self.src_adj[su].len();

                while self.s_idx[su] < num_edg {
                    let t = self.src_adj[su][self.s_idx[su]];
                    #[cfg(feature = "stats")] { ps.dfs_edg += 1; }
                    let ss = self.t_mate[t];
                    let skip = (self.t_level[t] != self.s_level[su] + 1) ||
                               ((ss == NIL) && (self.t_stt[t] == Stt::Last)) ||
                               ((ss != NIL) &&
                                  ((self.s_level[su] + 1 == shortest) ||
                                   (self.s_stt[ss as usize] == Stt::DfsActive) ||
                                   (self.s_stt[ss as usize] == Stt::DfsDone)));
                    if skip { self.s_idx[su] += 1; } else { break; }
                }

                if self.s_idx[su] < num_edg {
                    let t = self.src_adj[su][self.s_idx[su]];
                    let ss = self.t_mate[t];
                    #[cfg(feature = "stats")] { ps.dfs_vtx += 1; }
                    if ss == NIL {
                        self.s_last_queue.push(s);
                        self.t_last_queue.push(t as i32);
                        self.t_stt[t] = Stt::Last;
                        while !self.s_dfs_stack.is_empty() {
                            let top = self.s_dfs_stack.top();
                            self.s_stt[top as usize] = Stt::DfsDone;
                            self.s_dfs_stack.pop();
                        }
                        break;
                    }
                    #[cfg(feature = "stats")] { ps.dfs_edg += 1; ps.dfs_vtx += 1; }
                    let ssu = ss as usize;
                    self.s_ptr[ssu] = s;
                    self.s_dfs_stack.push(ss);
                    self.s_stt[ssu] = Stt::DfsActive;
                    self.s_idx[su] += 1;
                } else {
                    self.s_dfs_stack.pop();
                    self.s_stt[su] = Stt::DfsDone;
                }
            }
            sf = sf_next;
        }
    }

    /* DFS: find augmenting paths with lookahead (no BFS layering, matchbox-style) */
    fn dfs_find_paths_lkhd(
        &mut self,
        #[cfg(feature = "stats")] ps: &mut PhaseStats,
    ) {
        let mut sf = self.s_exposed.first();
        while sf != NIL {
            let sf_next = self.s_exposed.next(sf);
            #[cfg(feature = "stats")] { ps.dfs_vtx += 1; }
            self.s_dfs_stack.push(sf);
            self.s_stt[sf as usize] = Stt::DfsActive;

            while !self.s_dfs_stack.is_empty() {
                let s = self.s_dfs_stack.top();
                let su = s as usize;
                let num_edg = self.src_adj[su].len();

                /* Lookahead: check for direct free t-neighbor */
                while self.s_lkhd[su] < num_edg {
                    let t = self.src_adj[su][self.s_lkhd[su]];
                    #[cfg(feature = "stats")] { ps.bfs_edg += 1; }
                    let ss = self.t_mate[t];
                    let skip = (ss == NIL && self.t_stt[t] != Stt::Idle) || ss != NIL;
                    if skip { self.s_lkhd[su] += 1; } else { break; }
                }

                if self.s_lkhd[su] < num_edg {
                    let t = self.src_adj[su][self.s_lkhd[su]];
                    #[cfg(feature = "stats")] { ps.bfs_vtx += 1; }
                    self.s_last_queue.push(s);
                    self.t_last_queue.push(t as i32);
                    self.t_stt[t] = Stt::Last;
                    while !self.s_dfs_stack.is_empty() {
                        let sss = self.s_dfs_stack.top();
                        self.s_dfs_stack.pop();
                        self.s_done_queue.push(sss);
                        self.s_stt[sss as usize] = Stt::DfsDone;
                    }
                    self.s_lkhd[su] += 1;
                    break;
                }

                /* Regular DFS deeper */
                while self.s_idx[su] < num_edg {
                    let t = self.src_adj[su][self.s_idx[su]];
                    #[cfg(feature = "stats")] { ps.dfs_edg += 1; }
                    let ss = self.t_mate[t];
                    let skip = ss == NIL || (ss != NIL && self.s_stt[ss as usize] != Stt::Idle);
                    if skip { self.s_idx[su] += 1; } else { break; }
                }

                if self.s_idx[su] < num_edg {
                    let t = self.src_adj[su][self.s_idx[su]];
                    let ss = self.t_mate[t];
                    #[cfg(feature = "stats")] { ps.dfs_vtx += 1; ps.dfs_edg += 1; ps.dfs_vtx += 1; }
                    let ssu = ss as usize;
                    self.s_ptr[ssu] = s;
                    self.s_dfs_stack.push(ss);
                    self.s_stt[ssu] = Stt::DfsActive;
                    self.s_idx[su] += 1;
                } else {
                    self.s_dfs_stack.pop();
                    self.s_done_queue.push(s);
                    self.s_stt[su] = Stt::DfsDone;
                }
            }
            sf = sf_next;
        }
    }

    /* Augment along path from s_last to root via s_ptr */
    fn augment(&mut self, s_last: i32, t_last: i32) -> i32 {
        let mut s = s_last;
        let mut t = t_last;
        let mut k = 0i32;
        while s != NIL {
            let su = s as usize;
            let tt = self.s_mate[su];
            self.s_mate[su] = t;
            self.t_mate[t as usize] = s;
            if self.s_ptr[su] == NIL { self.s_exposed.erase(s); }
            s = self.s_ptr[su];
            t = tt;
            k += 1;
        }
        2 * k - 1
    }

    /* Cleanup after HK phase (also resets lookahead indices) */
    fn cleanup_hk(&mut self) {
        while !self.s_bfs_queue.is_empty() {
            let s = self.s_bfs_queue.front();
            self.s_bfs_queue.pop();
            let su = s as usize;
            self.s_ptr[su] = NIL; self.s_level[su] = INF_LEVEL;
            self.s_stt[su] = Stt::Idle;
            self.s_idx[su] = 0; self.s_lkhd[su] = 0;
        }
        while !self.s_done_queue.is_empty() {
            let s = self.s_done_queue.front();
            self.s_done_queue.pop();
            let su = s as usize;
            self.s_ptr[su] = NIL; self.s_level[su] = INF_LEVEL;
            self.s_stt[su] = Stt::Idle;
            self.s_idx[su] = 0; self.s_lkhd[su] = 0;
        }
        while !self.t_done_queue.is_empty() {
            let t = self.t_done_queue.front();
            self.t_done_queue.pop();
            let tu = t as usize;
            self.t_level[tu] = INF_LEVEL;
            self.t_stt[tu] = Stt::Idle;
        }
    }

    /* Main solver: HK mode */
    fn solve_hk(&mut self, greedy_mode: i32) -> i32 {
        self.greedy_size = match greedy_mode {
            1 => self.greedy_init(),
            2 => self.greedy_init_md(),
            _ => 0,
        };
        let mut card = self.greedy_size as i32;

        #[cfg(feature = "stats")]
        {
            self.stats.reset();
            self.stats.greedy_card = self.greedy_size as i32;
            self.stats.reversed = self.rvrs;
        }

        for s in 0..self.s_count {
            if self.s_mate[s] == NIL { self.s_exposed.push(s as i32); }
        }

        let mut phase_count = 0;
        loop {
            #[cfg(feature = "stats")]
            let mut ps = PhaseStats::new();

            let shortest = self.bfs_shortest_level(
                #[cfg(feature = "stats")] &mut ps,
            );

            if shortest == INF_LEVEL {
                self.cleanup_hk();
                break;
            }

            phase_count += 1;
            #[cfg(feature = "stats")] { ps.shortest_path_len = shortest; }

            self.dfs_find_paths_hk(shortest,
                #[cfg(feature = "stats")] &mut ps,
            );

            while !self.t_last_queue.is_empty() {
                let sl = self.s_last_queue.front(); self.s_last_queue.pop();
                let tl = self.t_last_queue.front(); self.t_last_queue.pop();
                let _path_len = self.augment(sl, tl);
                card += 1;
                #[cfg(feature = "stats")] {
                    ps.num_augmentations += 1;
                    ps.agg_aug_path_len += _path_len as i64;
                    if _path_len < ps.min_aug_path_len { ps.min_aug_path_len = _path_len; }
                    if _path_len > ps.max_aug_path_len { ps.max_aug_path_len = _path_len; }
                }
            }

            self.cleanup_hk();

            #[cfg(feature = "stats")] { self.stats.phases.push(ps); }
        }

        #[cfg(feature = "stats")] { self.stats.print(); }
        println!("Phases: {}", phase_count);
        card
    }

    /* Main solver: lookahead mode (DFS-only, matchbox-style) */
    fn solve_lkhd(&mut self, greedy_mode: i32) -> i32 {
        self.greedy_size = match greedy_mode {
            1 => self.greedy_init(),
            2 => self.greedy_init_md(),
            _ => 0,
        };
        let mut card = self.greedy_size as i32;

        #[cfg(feature = "stats")]
        {
            self.stats.reset();
            self.stats.greedy_card = self.greedy_size as i32;
            self.stats.reversed = self.rvrs;
        }

        for s in 0..self.s_count {
            if self.s_mate[s] == NIL { self.s_exposed.push(s as i32); }
        }

        let mut pass_count = 0;
        loop {
            #[cfg(feature = "stats")]
            let mut ps = PhaseStats::new();

            self.dfs_find_paths_lkhd(
                #[cfg(feature = "stats")] &mut ps,
            );

            if self.t_last_queue.is_empty() { break; }

            pass_count += 1;
            while !self.t_last_queue.is_empty() {
                let sl = self.s_last_queue.front(); self.s_last_queue.pop();
                let tl = self.t_last_queue.front(); self.t_last_queue.pop();
                self.t_stt[tl as usize] = Stt::Idle;
                let _path_len = self.augment(sl, tl);
                card += 1;
                #[cfg(feature = "stats")] {
                    ps.num_augmentations += 1;
                    ps.agg_aug_path_len += _path_len as i64;
                }
            }

            /* Cleanup: reset only visited vertices */
            while !self.s_done_queue.is_empty() {
                let s = self.s_done_queue.front();
                self.s_done_queue.pop();
                let su = s as usize;
                self.s_ptr[su] = NIL;
                self.s_stt[su] = Stt::Idle;
                self.s_idx[su] = 0; self.s_lkhd[su] = 0;
            }

            #[cfg(feature = "stats")] { self.stats.phases.push(ps); }
        }

        #[cfg(feature = "stats")] { self.stats.print(); }
        println!("Phases: {}", pass_count);
        card
    }

    /* Get matching in original (left, right) coordinates */
    fn get_matching(&self) -> Vec<(usize, usize)> {
        let mut result: Vec<(usize, usize)> = Vec::new();
        for s in 0..self.s_count {
            if self.s_mate[s] != NIL {
                if self.rvrs { result.push((self.s_mate[s] as usize, s)); }
                else         { result.push((s, self.s_mate[s] as usize)); }
            }
        }
        result.sort_unstable();
        result
    }
}

/* =========================================================================
 * Validation
 * ========================================================================= */
fn validate_matching(
    left_count: usize, right_count: usize,
    graph: &[Vec<usize>], matching: &[(usize, usize)], rvrs: bool,
) {
    let mut ldeg = vec![0i32; left_count];
    let mut rdeg = vec![0i32; right_count];
    let mut errors = 0;

    for &(u, v) in matching {
        let adj = if rvrs { &graph[v] } else { &graph[u] };
        let target = if rvrs { u } else { v };
        if !adj.binary_search(&target).is_ok() {
            eprintln!("ERROR: Edge ({}, {}) not in graph!", u, v);
            errors += 1;
        }
        ldeg[u] += 1;
        rdeg[v] += 1;
    }
    for i in 0..left_count {
        if ldeg[i] > 1 { eprintln!("ERROR: Left {} in {} edges!", i, ldeg[i]); errors += 1; }
    }
    for i in 0..right_count {
        if rdeg[i] > 1 { eprintln!("ERROR: Right {} in {} edges!", i, rdeg[i]); errors += 1; }
    }
    let ml = ldeg.iter().filter(|&&d| d > 0).count();
    let mr = rdeg.iter().filter(|&&d| d > 0).count();

    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.len());
    println!("Left matched: {}, Right matched: {}", ml, mr);
    println!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
    println!("=========================\n");
}

fn load_graph(filename: &str) -> Result<(usize, usize, Vec<(usize, usize)>), Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let first = lines.next().ok_or("Empty file")??;
    let parts: Vec<&str> = first.split_whitespace().collect();
    if parts.len() != 3 {
        return Err("First line must have 3 numbers".into());
    }
    let left_count: usize = parts[0].parse()?;
    let right_count: usize = parts[1].parse()?;
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
    Ok((left_count, right_count, edges))
}

fn main() {
    println!("Hopcroft-Karp Pure Algorithm - Rust Implementation");
    println!("====================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename> [--greedy|--greedy-md] [--lkhd] [--rvrs|--fwd]", args[0]);
        std::process::exit(1);
    }

    let greedy_mode: i32 =
        if args.iter().any(|a| a == "--greedy-md") { 2 }
        else if args.iter().any(|a| a == "--greedy") { 1 }
        else { 0 };
    let use_lkhd = args.iter().any(|a| a == "--lkhd");
    let forced_rvrs = args.iter().any(|a| a == "--rvrs");
    let forced_fwd  = args.iter().any(|a| a == "--fwd");
    let auto_dir = !forced_rvrs && !forced_fwd;
    let mut use_rvrs = forced_rvrs;

    match load_graph(&args[1]) {
        Ok((lc, rc, edges)) => {
            println!("Graph: {} left, {} right, {} edges", lc, rc, edges.len());
            if auto_dir { use_rvrs = rc < lc; }

            println!("Mode: {}", if use_lkhd { "lookahead DFS" } else { "HK (BFS+DFS)" });
            println!("Direction: {}{}",
                     if use_rvrs { "T->S (reversed)" } else { "S->T (normal)" },
                     if auto_dir { " (auto)" } else { " (forced)" });

            let start = Instant::now();
            let mut hk = HopcroftKarpPure::new();
            hk.build(lc, rc, &edges, use_rvrs);
            let card = if use_lkhd { hk.solve_lkhd(greedy_mode) }
                       else        { hk.solve_hk(greedy_mode) };
            let duration = start.elapsed();

            let matching = hk.get_matching();
            validate_matching(lc, rc, &hk.src_adj, &matching, use_rvrs);

            println!("Matching size: {}", card);
            if greedy_mode > 0 {
                println!("Greedy init size: {}", hk.greedy_size);
                if card > 0 { println!("Greedy/Final: {:.2}%", 100.0 * hk.greedy_size as f64 / card as f64); }
            }
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

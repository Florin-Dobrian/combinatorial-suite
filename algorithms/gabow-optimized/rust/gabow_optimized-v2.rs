/*
 * Gabow's Scaling Algorithm (Optimized) - O(E√V) Maximum Matching
 *
 * Pure cardinality (unweighted) version — integer weights conceptually all 1.
 *
 * Phase 1: BFS by levels (Delta), detect blossoms.
 *          Build contracted graph H: edges connecting different dbase
 *          components that were processed during BFS.
 * Phase 2: Find all shortest augmenting paths in H (iterative DFS
 *          with blossom contraction), unfold to G via bridges.
 *
 * Based on LEDA-7's mc_matching_gabow architecture, stripped of weighted
 * dual machinery. Rust implementation — fully deterministic, no hash containers.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;
const UNLABELED: i32 = 0;
const EVEN: i32 = 1;
const ODD: i32 = 2;

struct GabowOptimized {
    n: usize,
    graph: Vec<Vec<usize>>,
    mate: Vec<i32>,

    label: Vec<i32>,
    parent: Vec<i32>,
    source_bridge: Vec<i32>,
    target_bridge: Vec<i32>,

    base_par: Vec<usize>,
    dbase_par: Vec<usize>,

    level_queue: Vec<Vec<(usize, usize)>>,

    lca_tag1: Vec<usize>,
    lca_tag2: Vec<usize>,
    lca_epoch: usize,

    in_tree: Vec<bool>,
    tree_nodes: Vec<usize>,
    delta: i32,

    rep: Vec<usize>,
    mate_h: Vec<i32>,
    label_h: Vec<i32>,
    parent_h_src: Vec<i32>,
    parent_h_tgt: Vec<i32>,
    bridge_h_src: Vec<i32>,
    bridge_h_tgt: Vec<i32>,
    dir_h: Vec<i32>,
    even_time_h: Vec<i32>,
    t_h: i32,
    db2_par: Vec<usize>,
    contracted_into: Vec<Vec<usize>>,
}

impl GabowOptimized {
    fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut graph = vec![Vec::new(); n];
        for &(u, v) in edges {
            if u < n && v < n && u != v {
                graph[u].push(v);
                graph[v].push(u);
            }
        }
        for adj in &mut graph { adj.sort_unstable(); adj.dedup(); }

        GabowOptimized {
            n, graph,
            mate: vec![NIL; n],
            label: vec![UNLABELED; n],
            parent: vec![NIL; n],
            source_bridge: vec![NIL; n],
            target_bridge: vec![NIL; n],
            base_par: (0..n).collect(),
            dbase_par: (0..n).collect(),
            level_queue: vec![Vec::new(); n + 2],
            lca_tag1: vec![0; n],
            lca_tag2: vec![0; n],
            lca_epoch: 0,
            in_tree: vec![false; n],
            tree_nodes: Vec::new(),
            delta: 0,
            rep: vec![0; n],
            mate_h: vec![NIL; n],
            label_h: vec![UNLABELED; n],
            parent_h_src: vec![NIL; n],
            parent_h_tgt: vec![NIL; n],
            bridge_h_src: vec![NIL; n],
            bridge_h_tgt: vec![NIL; n],
            dir_h: vec![0; n],
            even_time_h: vec![0; n],
            t_h: 0,
            db2_par: (0..n).collect(),
            contracted_into: vec![Vec::new(); n],
        }
    }

    /* ---- union-find: base ---- */
    fn find_base(&mut self, mut v: usize) -> usize {
        while self.base_par[v] != v {
            self.base_par[v] = self.base_par[self.base_par[v]];
            v = self.base_par[v];
        }
        v
    }
    fn union_base(&mut self, a: usize, b: usize, r: usize) {
        let fa = self.find_base(a);
        let fb = self.find_base(b);
        self.base_par[fa] = r;
        self.base_par[fb] = r;
    }

    /* ---- union-find: dbase ---- */
    fn find_dbase(&mut self, mut v: usize) -> usize {
        while self.dbase_par[v] != v {
            self.dbase_par[v] = self.dbase_par[self.dbase_par[v]];
            v = self.dbase_par[v];
        }
        v
    }
    fn union_dbase(&mut self, a: usize, b: usize) {
        let fa = self.find_dbase(a);
        let fb = self.find_dbase(b);
        if fa != fb { self.dbase_par[fa] = fb; }
    }
    fn make_rep_dbase(&mut self, v: usize) {
        let r = self.find_dbase(v);
        if r != v { self.dbase_par[r] = v; self.dbase_par[v] = v; }
    }

    /* ---- union-find: dbase2 ---- */
    fn find_db2(&mut self, mut v: usize) -> usize {
        while self.db2_par[v] != v {
            self.db2_par[v] = self.db2_par[self.db2_par[v]];
            v = self.db2_par[v];
        }
        v
    }
    fn union_db2(&mut self, a: usize, b: usize) {
        let fa = self.find_db2(a);
        let fb = self.find_db2(b);
        if fa != fb { self.db2_par[fa] = fb; }
    }
    fn make_rep_db2(&mut self, v: usize) {
        let r = self.find_db2(v);
        if r != v { self.db2_par[r] = v; self.db2_par[v] = v; }
    }

    /* ---- interleaved LCA ---- */
    fn find_lca(&mut self, u: usize, v: usize) -> i32 {
        self.lca_epoch += 1;
        let ep = self.lca_epoch;
        let mut hx = self.find_base(u);
        let mut hy = self.find_base(v);
        self.lca_tag1[hx] = ep;
        self.lca_tag2[hy] = ep;
        loop {
            if self.lca_tag1[hy] == ep { return hy as i32; }
            if self.lca_tag2[hx] == ep { return hx as i32; }
            let hxr = self.mate[hx] == NIL || self.parent[self.mate[hx] as usize] == NIL;
            let hyr = self.mate[hy] == NIL || self.parent[self.mate[hy] as usize] == NIL;
            if hxr && hyr { return NIL; }
            if !hxr {
                hx = self.find_base(self.parent[self.mate[hx] as usize] as usize);
                self.lca_tag1[hx] = ep;
            }
            if !hyr {
                hy = self.find_base(self.parent[self.mate[hy] as usize] as usize);
                self.lca_tag2[hy] = ep;
            }
        }
    }

    /* ---- shrink_path ---- */
    fn shrink_path(&mut self, b: usize, x: usize, y: usize,
                   dunions: &mut Vec<(usize, usize)>) {
        let mut v = self.find_base(x);
        while v != b {
            self.union_base(v, b, b);
            dunions.push((v, b));
            let mv = self.mate[v] as usize;
            self.union_base(mv, b, b);
            dunions.push((mv, b));
            self.base_par[b] = b;
            self.source_bridge[mv] = x as i32;
            self.target_bridge[mv] = y as i32;
            let d = self.delta;
            let neighbors: Vec<usize> = self.graph[mv].clone();
            for w in neighbors {
                if w as i32 == self.mate[mv] { continue; }
                let bw = self.find_base(w);
                if self.label[bw] == ODD { continue; }
                if self.label[bw] == UNLABELED {
                    self.level_queue[(d + 1) as usize].push((mv, w));
                } else if self.label[bw] == EVEN {
                    self.level_queue[d as usize].push((mv, w));
                }
            }
            v = self.find_base(self.parent[mv] as usize);
        }
        dunions.push((b, b));
    }

    /* ================================================================ */
    /*                          PHASE 1                                 */
    /* ================================================================ */
    fn phase_1(&mut self) -> bool {
        self.delta = 0;
        self.tree_nodes.clear();
        for q in &mut self.level_queue { q.clear(); }
        let mut dunions: Vec<(usize, usize)> = Vec::new();

        for i in 0..self.n {
            self.base_par[i] = i;
            self.dbase_par[i] = i;
            self.label[i] = UNLABELED;
            self.parent[i] = NIL;
            self.source_bridge[i] = NIL;
            self.target_bridge[i] = NIL;
            self.in_tree[i] = false;
        }

        /* Free vertices are EVEN roots at Delta=0 */
        for v in 0..self.n {
            if self.mate[v] == NIL {
                self.label[v] = EVEN;
                self.in_tree[v] = true;
                self.tree_nodes.push(v);
                let neighbors: Vec<usize> = self.graph[v].clone();
                for u in neighbors {
                    if u as i32 == self.mate[v] { continue; }
                    let bu = self.find_base(u);
                    if self.label[bu] == ODD { continue; }
                    if self.label[bu] == UNLABELED {
                        self.level_queue[1].push((v, u));
                    } else if self.label[bu] == EVEN {
                        self.level_queue[0].push((v, u));
                    }
                }
            }
        }

        let mut found_sap = false;

        while self.delta <= self.n as i32 {
            let d = self.delta as usize;
            while !self.level_queue[d].is_empty() {
                let (mut z, mut u) = self.level_queue[d].pop().unwrap();
                let mut bz = self.find_base(z);
                let mut bu = self.find_base(u);
                if self.label[bz] != EVEN {
                    std::mem::swap(&mut z, &mut u);
                    std::mem::swap(&mut bz, &mut bu);
                }
                if bz == bu || self.label[bz] != EVEN { continue; }
                if u as i32 == self.mate[z] || self.label[bu] == ODD { continue; }

                if self.label[bu] == UNLABELED {
                    let mv = self.mate[u];
                    if mv == NIL { continue; }
                    let mv = mv as usize;
                    self.parent[u] = z as i32;
                    self.parent[mv] = u as i32;
                    self.label[u] = ODD;
                    self.label[mv] = EVEN;
                    self.in_tree[u] = true;
                    self.in_tree[mv] = true;
                    self.tree_nodes.push(u);
                    self.tree_nodes.push(mv);
                    let neighbors: Vec<usize> = self.graph[mv].clone();
                    let delta = self.delta;
                    for w in neighbors {
                        if w as i32 == self.mate[mv] { continue; }
                        let bw = self.find_base(w);
                        if self.label[bw] == ODD { continue; }
                        if self.label[bw] == UNLABELED {
                            self.level_queue[(delta + 1) as usize].push((mv, w));
                        } else if self.label[bw] == EVEN {
                            self.level_queue[delta as usize].push((mv, w));
                        }
                    }
                } else if self.label[bu] == EVEN {
                    let lca = self.find_lca(z, u);
                    if lca != NIL {
                        let lca = lca as usize;
                        self.shrink_path(lca, z, u, &mut dunions);
                        self.shrink_path(lca, u, z, &mut dunions);
                    } else {
                        found_sap = true;
                    }
                }
            }

            if found_sap {
                /* Build H: contracted_into and mateH */
                let tn: Vec<usize> = self.tree_nodes.clone();
                for &v in &tn {
                    let db = self.find_dbase(v);
                    self.contracted_into[db].push(v);
                    self.mate_h[v] = NIL;
                }
                for &u in &tn {
                    let uh = self.find_dbase(u);
                    let mv = self.mate[u];
                    if mv != NIL && self.in_tree[mv as usize] {
                        let vh = self.find_dbase(mv as usize);
                        if uh != vh {
                            self.mate_h[uh] = vh as i32;
                            self.mate_h[vh] = uh as i32;
                        }
                    }
                }
                return true;
            }

            for (a, b) in dunions.drain(..) {
                if a == b { self.make_rep_dbase(a); }
                else { self.union_dbase(a, b); }
            }
            self.delta += 1;
        }
        false
    }

    /* ================================================================ */
    /*                          PHASE 2                                 */
    /* ================================================================ */

    /* find_apHG: ITERATIVE DFS in H.
     * Scans graph[v] for each G-vertex in contracted_into[vh].
     * Returns the free H-node found, or NIL. */
    fn find_ap_hg(&mut self, root_vh: usize) -> i32 {
        struct Frame { vh: usize, ci_idx: usize, adj_idx: usize }
        let mut stk: Vec<Frame> = vec![Frame { vh: root_vh, ci_idx: 0, adj_idx: 0 }];

        'outer: while let Some(f) = stk.last_mut() {
            let vh = f.vh;
            while f.ci_idx < self.contracted_into[vh].len() {
                let v = self.contracted_into[vh][f.ci_idx];
                while f.adj_idx < self.graph[v].len() {
                    let w = self.graph[v][f.adj_idx];
                    f.adj_idx += 1;

                    if !self.in_tree[w] { continue; }
                    if self.mate[v] == w as i32 { continue; }
                    let dv = self.find_dbase(v);
                    let dw = self.find_dbase(w);
                    if dv == dw { continue; }
                    let uh = self.find_db2(self.rep[w]);
                    if self.mate_h[vh] != NIL && self.mate_h[vh] as usize == uh { continue; }
                    if self.label_h[uh] == ODD { continue; }

                    if self.label_h[uh] == UNLABELED {
                        let muh = self.mate_h[uh];
                        if muh == NIL {
                            self.label_h[uh] = ODD;
                            self.parent_h_src[uh] = w as i32;
                            self.parent_h_tgt[uh] = v as i32;
                            return uh as i32;
                        }
                        self.label_h[uh] = ODD;
                        self.parent_h_src[uh] = w as i32;
                        self.parent_h_tgt[uh] = v as i32;
                        let muh = muh as usize;
                        self.label_h[muh] = EVEN;
                        self.even_time_h[muh] = self.t_h;
                        self.t_h += 1;
                        stk.push(Frame { vh: muh, ci_idx: 0, adj_idx: 0 });
                        continue 'outer;
                    } else if self.label_h[uh] == EVEN {
                        let bh = self.find_db2(vh);
                        let zh = self.find_db2(uh);
                        if self.even_time_h[bh] < self.even_time_h[zh] {
                            let mut tmp: Vec<usize> = Vec::new();
                            let mut endpoints: Vec<usize> = Vec::new();
                            let mut cur = zh;
                            while cur != bh {
                                endpoints.push(cur);
                                let mc = self.mate_h[cur] as usize;
                                endpoints.push(mc);
                                tmp.push(mc);
                                let ps = self.parent_h_src[mc] as usize;
                                let pt = self.parent_h_tgt[mc] as usize;
                                let next = if self.rep[ps] == mc {
                                    self.rep[pt]
                                } else {
                                    self.rep[ps]
                                };
                                cur = self.find_db2(next);
                            }
                            for &nd in &endpoints { self.union_db2(nd, bh); }
                            self.make_rep_db2(bh);
                            for &mc in &tmp {
                                self.bridge_h_src[mc] = v as i32;
                                self.bridge_h_tgt[mc] = w as i32;
                                self.dir_h[mc] = -1;
                            }
                            for i in (0..tmp.len()).rev() {
                                stk.push(Frame { vh: tmp[i], ci_idx: 0, adj_idx: 0 });
                            }
                            continue 'outer;
                        }
                    }
                }
                f.ci_idx += 1;
                f.adj_idx = 0;
            }
            stk.pop();
        }
        NIL
    }

    /* trace_h_path: iterative, collects non-matching G-edges along H-path */
    fn trace_h_path(&mut self, start_vh: usize, start_uh: usize,
                    edges_out: &mut Vec<(usize, usize)>) {
        struct Frame { vh: usize, uh: usize, phase: i32, bs: usize, bt: usize,
                       side_a: usize, side_b: usize }
        let mut stk = vec![Frame { vh: start_vh, uh: start_uh, phase: 0,
                                   bs: 0, bt: 0, side_a: 0, side_b: 0 }];
        while let Some(f) = stk.last_mut() {
            if f.vh == f.uh { stk.pop(); continue; }
            if self.label_h[f.vh] == EVEN {
                let mvh = self.mate_h[f.vh] as usize;
                let ps = self.parent_h_src[mvh] as usize;
                let pt = self.parent_h_tgt[mvh] as usize;
                edges_out.push((ps, pt));
                f.vh = if self.rep[ps] == mvh { self.rep[pt] } else { self.rep[ps] };
                continue;
            }
            if f.phase == 0 {
                let bs = self.bridge_h_src[f.vh] as usize;
                let bt = self.bridge_h_tgt[f.vh] as usize;
                f.bs = bs; f.bt = bt;
                if self.dir_h[f.vh] == 1 {
                    f.side_a = self.rep[bs]; f.side_b = self.rep[bt];
                } else {
                    f.side_a = self.rep[bt]; f.side_b = self.rep[bs];
                }
                f.phase = 1;
                let mt = if self.mate_h[f.vh] != NIL {
                    self.rep[self.mate_h[f.vh] as usize]
                } else { f.vh };
                let sa = f.side_a;
                stk.push(Frame { vh: sa, uh: mt, phase: 0,
                                 bs: 0, bt: 0, side_a: 0, side_b: 0 });
                continue;
            }
            if f.phase == 1 {
                edges_out.push((f.bs, f.bt));
                f.phase = 2;
                let sb = f.side_b; let uh = f.uh;
                stk.push(Frame { vh: sb, uh, phase: 0,
                                 bs: 0, bt: 0, side_a: 0, side_b: 0 });
                continue;
            }
            stk.pop();
        }
    }

    /* find_path_in_g: iterative unfold within single H-node */
    fn find_path_in_g(&mut self, start_v: usize, start_u: usize,
                      pairs: &mut Vec<(usize, usize)>) {
        struct Frame { v: usize, u: usize, phase: i32, sb: usize, tb: usize }
        let mut stk = vec![Frame { v: start_v, u: start_u, phase: 0, sb: 0, tb: 0 }];
        while let Some(f) = stk.last_mut() {
            if f.v == f.u { stk.pop(); continue; }
            if f.phase == 0 {
                if self.label[f.v] == EVEN {
                    let mv = self.mate[f.v] as usize;
                    let pmv = self.parent[mv] as usize;
                    pairs.push((mv, pmv));
                    f.v = pmv;
                    continue;
                }
                let sb = self.source_bridge[f.v] as usize;
                let tb = self.target_bridge[f.v] as usize;
                let mv = self.mate[f.v] as usize;
                f.sb = sb; f.tb = tb; f.phase = 1;
                stk.push(Frame { v: sb, u: mv, phase: 0, sb: 0, tb: 0 });
                continue;
            }
            if f.phase == 1 {
                pairs.push((f.sb, f.tb));
                f.phase = 2;
                let tb = f.tb; let u = f.u;
                stk.push(Frame { v: tb, u, phase: 0, sb: 0, tb: 0 });
                continue;
            }
            stk.pop();
        }
    }

    /* augment_g: unfold H-edges to G and augment */
    fn augment_g(&mut self, h_edges: &[(usize, usize)]) {
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        for &(u, v) in h_edges {
            pairs.push((u, v));
            let ru = self.rep[u];
            self.find_path_in_g(u, ru, &mut pairs);
            let rv = self.rep[v];
            self.find_path_in_g(v, rv, &mut pairs);
        }
        for &(a, b) in &pairs {
            self.mate[a] = b as i32;
            self.mate[b] = a as i32;
        }
    }

    /* phase_2: find all SAPs in H, unfold and augment */
    fn phase_2(&mut self) {
        let tn: Vec<usize> = self.tree_nodes.clone();
        for &v in &tn {
            let db = self.find_dbase(v);
            self.rep[v] = db;
            self.label_h[v] = UNLABELED;
            self.parent_h_src[v] = NIL; self.parent_h_tgt[v] = NIL;
            self.bridge_h_src[v] = NIL; self.bridge_h_tgt[v] = NIL;
            self.dir_h[v] = 0;
            self.even_time_h[v] = 0;
            self.db2_par[v] = v;
        }
        self.t_h = 0;

        let mut all_paths: Vec<Vec<(usize, usize)>> = Vec::new();
        for &vh in &tn {
            if vh != self.rep[vh] { continue; }
            if self.label_h[vh] != UNLABELED || self.mate_h[vh] != NIL { continue; }
            self.label_h[vh] = EVEN;
            self.even_time_h[vh] = self.t_h;
            self.t_h += 1;

            let free_node = self.find_ap_hg(vh);
            if free_node != NIL {
                let free_node = free_node as usize;
                let mut h_nm: Vec<(usize, usize)> = Vec::new();
                let ps = self.parent_h_src[free_node] as usize;
                let pt = self.parent_h_tgt[free_node] as usize;
                h_nm.push((ps, pt));
                let next = if self.rep[ps] == free_node {
                    self.rep[pt]
                } else {
                    self.rep[ps]
                };
                self.trace_h_path(next, vh, &mut h_nm);
                all_paths.push(h_nm);
            }
        }

        for path in &all_paths { self.augment_g(path); }

        /* Clean up */
        for &v in &tn {
            let db = self.find_dbase(v);
            self.contracted_into[db].clear();
            self.contracted_into[v].clear();
            self.mate_h[v] = NIL;
        }
    }

    /* ================================================================ */
    /*                      MAIN ENTRY POINT                            */
    /* ================================================================ */
    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        /* greedy init */
        for u in 0..self.n {
            if self.mate[u] != NIL { continue; }
            let neighbors: Vec<usize> = self.graph[u].clone();
            for v in neighbors {
                if self.mate[v] == NIL {
                    self.mate[u] = v as i32;
                    self.mate[v] = u as i32;
                    break;
                }
            }
        }
        while self.phase_1() { self.phase_2(); }

        let mut result = Vec::new();
        for u in 0..self.n {
            if self.mate[u] != NIL && (self.mate[u] as usize) > u {
                result.push((u, self.mate[u] as usize));
            }
        }
        result.sort_unstable();
        result
    }
}

/* ================================================================ */
/*                    VALIDATION AND MAIN                            */
/* ================================================================ */

fn validate_matching(n: usize, graph: &[Vec<usize>], matching: &[(usize, usize)]) {
    let mut deg = vec![0i32; n];
    let mut errors = 0;
    for &(u, v) in matching {
        if graph[u].binary_search(&v).is_err() {
            eprintln!("ERROR: Edge ({}, {}) not in graph!", u, v);
            errors += 1;
        }
        deg[u] += 1; deg[v] += 1;
    }
    for i in 0..n {
        if deg[i] > 1 { eprintln!("ERROR: Vertex {} in {} edges!", i, deg[i]); errors += 1; }
    }
    let matched = deg.iter().filter(|&&d| d > 0).count();
    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.len());
    println!("Matched vertices: {}", matched);
    println!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
    println!("=========================\n");
}

fn load_graph(filename: &str) -> Result<(usize, Vec<(usize, usize)>), Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let first = lines.next().ok_or("Empty file")??;
    let parts: Vec<&str> = first.split_whitespace().collect();
    let n: usize = parts[0].parse()?;
    let m: usize = parts[1].parse()?;
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
    Ok((n, edges))
}

fn main() {
    println!("Gabow's Scaling Algorithm (Optimized) - Rust Implementation");
    println!("=============================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename>", args[0]);
        std::process::exit(1);
    }

    match load_graph(&args[1]) {
        Ok((n, edges)) => {
            println!("Graph: {} vertices, {} edges", n, edges.len());
            let start = Instant::now();
            let mut gabow = GabowOptimized::new(n, &edges);
            let matching = gabow.maximum_matching();
            let duration = start.elapsed();
            validate_matching(n, &gabow.graph, &matching);
            println!("Matching size: {}", matching.len());
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

/*
 * Edmonds' Blossom Algorithm (Simple) — Unweighted Maximum Cardinality Matching
 *
 * Rust implementation. Single-source BFS (tree, not forest). Blossom IDs
 * reset to n each BFS. All indices are i32 (same type as vertex indices).
 *
 * Complexity: O(V^2 * E) worst case.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

// ── Blossom data ─────────────────────────────────────────────────────

#[derive(Clone)]
struct Blos {
    childs: Vec<i32>,          // sub-blossom IDs in cycle order
    edges: Vec<(i32, i32)>,    // edges[i] connects childs[i] to childs[(i+1)%k]
}

impl Blos {
    fn new() -> Self { Blos { childs: Vec::new(), edges: Vec::new() } }
}

// ── Solver ───────────────────────────────────────────────────────────

struct Solver {
    n: i32,
    adj: Vec<Vec<i32>>,
    mate: Vec<i32>,

    blos: Vec<Blos>,
    nblos: i32,

    inblossom: Vec<i32>,
    blossomparent: Vec<i32>,
    blossombase: Vec<i32>,

    label: Vec<i32>,              // 0=none, 1=S, 2=T, 5=breadcrumb
    labeledge: Vec<(i32, i32)>,
    queue: Vec<i32>,

    greedy_size: i32,
}

impl Solver {
    fn new(n: i32, edges: &[(i32, i32)]) -> Self {
        let nu = n as usize;
        let mut adj = vec![Vec::new(); nu];
        for &(u, v) in edges {
            if u != v && u >= 0 && u < n && v >= 0 && v < n {
                adj[u as usize].push(v);
                adj[v as usize].push(u);
            }
        }
        for a in &mut adj {
            a.sort_unstable();
            a.dedup();
        }

        let mut inblossom = vec![0i32; nu];
        let mut blossombase = vec![0i32; nu];
        let blossomparent = vec![-1i32; nu];
        for i in 0..n {
            inblossom[i as usize] = i;
            blossombase[i as usize] = i;
        }

        Solver {
            n, adj, mate: vec![-1; nu],
            blos: vec![Blos::new(); nu],
            nblos: n,
            inblossom, blossomparent, blossombase,
            label: Vec::new(), labeledge: Vec::new(), queue: Vec::new(),
            greedy_size: 0,
        }
    }

    fn is_blossom(&self, b: i32) -> bool { b >= self.n }

    fn ensure(&mut self, b: i32) {
        let needed = (b + 1) as usize;
        if self.label.len() < needed {
            self.label.resize(needed, 0);
            self.labeledge.resize(needed, (-1, -1));
            self.blossomparent.resize(needed, -1);
            self.blossombase.resize(needed, -1);
        }
    }

    fn leaves(&self, b: i32, out: &mut Vec<i32>) {
        if !self.is_blossom(b) {
            out.push(b);
            return;
        }
        for &c in &self.blos[b as usize].childs {
            self.leaves(c, out);
        }
    }

    // ── Reset for new BFS ────────────────────────────────────────────

    fn reset_blossoms(&mut self) {
        self.nblos = self.n;
        self.blos.truncate(self.n as usize);
        let nu = self.n as usize;
        self.inblossom.resize(nu, 0);
        self.blossombase.resize(nu, 0);
        self.blossomparent.resize(nu, -1);
        for i in 0..self.n {
            self.inblossom[i as usize] = i;
            self.blossombase[i as usize] = i;
            self.blossomparent[i as usize] = -1;
        }
        self.label = vec![0; nu];
        self.labeledge = vec![(-1, -1); nu];
        self.queue.clear();
    }

    // ── Tree building ────────────────────────────────────────────────

    fn assign_label(&mut self, w: i32, t: i32, v: i32) {
        let b = self.inblossom[w as usize];
        self.ensure(b);
        self.label[b as usize] = t;
        self.label[w as usize] = t;
        if v != -1 {
            self.labeledge[w as usize] = (v, w);
            self.labeledge[b as usize] = (v, w);
        } else {
            self.labeledge[w as usize] = (-1, -1);
            self.labeledge[b as usize] = (-1, -1);
        }
        if t == 1 {
            let mut lv = Vec::new();
            self.leaves(b, &mut lv);
            for u in lv { self.queue.push(u); }
        } else if t == 2 {
            let base = self.blossombase[b as usize];
            let mb = self.mate[base as usize];
            self.assign_label(mb, 1, base);
        }
    }

    // ── Blossom detection ────────────────────────────────────────────

    fn scan_blossom(&mut self, mut v: i32, mut w: i32) -> i32 {
        let mut path: Vec<i32> = Vec::new();
        let mut base = -2i32;
        let mut v_active = true;
        let mut w_active = true;

        loop {
            if !v_active && !w_active { break; }
            if v_active {
                let b = self.inblossom[v as usize];
                if self.label[b as usize] & 4 != 0 {
                    base = self.blossombase[b as usize];
                    break;
                }
                path.push(b);
                self.label[b as usize] = 5;
                let le = self.labeledge[b as usize];
                if le.0 == -1 {
                    v_active = false;
                } else {
                    v = le.0;
                    let bt = self.inblossom[v as usize];
                    v = self.labeledge[bt as usize].0;
                }
                if w_active { std::mem::swap(&mut v, &mut w); std::mem::swap(&mut v_active, &mut w_active); }
            } else {
                std::mem::swap(&mut v, &mut w);
                std::mem::swap(&mut v_active, &mut w_active);
            }
        }
        for &b in &path { self.label[b as usize] = 1; }
        base
    }

    // ── Blossom contraction ──────────────────────────────────────────

    fn add_blossom(&mut self, base: i32, mut v: i32, mut w: i32) {
        let bb = self.inblossom[base as usize];
        let mut bv = self.inblossom[v as usize];
        let mut bw = self.inblossom[w as usize];

        let bid = self.nblos;
        self.nblos += 1;
        if (bid as usize) >= self.blos.len() {
            self.blos.push(Blos::new());
        } else {
            self.blos[bid as usize].childs.clear();
            self.blos[bid as usize].edges.clear();
        }
        self.ensure(bid);
        self.blossombase[bid as usize] = base;
        self.blossomparent[bid as usize] = -1;
        self.blossomparent[bb as usize] = bid;

        // Trace from v to base
        let mut childs: Vec<i32> = Vec::new();
        let mut edges: Vec<(i32, i32)> = Vec::new();
        edges.push((v, w)); // bridge

        while bv != bb {
            self.blossomparent[bv as usize] = bid;
            childs.push(bv);
            edges.push(self.labeledge[bv as usize]);
            v = self.labeledge[bv as usize].0;
            bv = self.inblossom[v as usize];
        }
        childs.push(bb);
        childs.reverse();
        edges.reverse();

        // Trace from w to base
        while bw != bb {
            self.blossomparent[bw as usize] = bid;
            childs.push(bw);
            let le = self.labeledge[bw as usize];
            edges.push((le.1, le.0)); // reversed
            w = self.labeledge[bw as usize].0;
            bw = self.inblossom[w as usize];
        }

        self.blos[bid as usize].childs = childs;
        self.blos[bid as usize].edges = edges;

        self.label[bid as usize] = 1;
        self.labeledge[bid as usize] = self.labeledge[bb as usize];

        // Relabel: T-vertices become S
        let mut lv = Vec::new();
        self.leaves(bid, &mut lv);
        for u in lv {
            if self.label[self.inblossom[u as usize] as usize] == 2 {
                self.queue.push(u);
            }
            self.inblossom[u as usize] = bid;
        }
    }

    // ── Blossom expansion ────────────────────────────────────────────

    fn expand_blossom(&mut self, b: i32, endstage: bool) {
        struct Frame { b: i32, endstage: bool, idx: usize }
        let mut stack = vec![Frame { b, endstage, idx: 0 }];

        while !stack.is_empty() {
            let si = stack.len() - 1;
            let bl_childs_len = self.blos[stack[si].b as usize].childs.len();
            if stack[si].idx < bl_childs_len {
                let s = self.blos[stack[si].b as usize].childs[stack[si].idx];
                stack[si].idx += 1;
                self.blossomparent[s as usize] = -1;
                if self.is_blossom(s) {
                    if stack[si].endstage {
                        stack.push(Frame { b: s, endstage: true, idx: 0 });
                        continue;
                    } else {
                        let mut lv = Vec::new();
                        self.leaves(s, &mut lv);
                        for u in lv { self.inblossom[u as usize] = s; }
                    }
                } else {
                    self.inblossom[s as usize] = s;
                }
            } else {
                // All children processed
                if !stack[si].endstage && self.label[stack[si].b as usize] == 2 {
                    // Mid-stage T-blossom expansion: relabel children
                    let fb = stack[si].b;
                    self.relabel_expanded_t_blossom(fb);
                }
                let fb = stack[si].b;
                self.label[fb as usize] = 0;
                self.blos[fb as usize].childs.clear();
                self.blos[fb as usize].edges.clear();
                stack.pop();
            }
        }
    }

    fn relabel_expanded_t_blossom(&mut self, b: i32) {
        let entrychild = self.inblossom[self.labeledge[b as usize].1 as usize];
        let childs = self.blos[b as usize].childs.clone();
        let edges = self.blos[b as usize].edges.clone();
        let k = childs.len() as i32;

        let mut j: i32 = 0;
        for i in 0..childs.len() {
            if childs[i] == entrychild { j = i as i32; break; }
        }

        let jstep: i32;
        if j & 1 != 0 { j -= k; jstep = 1; } else { jstep = -1; }

        let (mut lv_, mut lw_) = self.labeledge[b as usize];

        while j != 0 {
            let (_, qq) = if jstep == 1 {
                let idx = ((j % k + k) % k) as usize;
                edges[idx]
            } else {
                let ei = (((j - 1) % k + k) % k) as usize;
                (edges[ei].1, edges[ei].0)
            };
            self.label[lw_ as usize] = 0;
            self.label[qq as usize] = 0;
            self.assign_label(lw_, 2, lv_);
            j += jstep;
            let (nlv, nlw) = if jstep == 1 {
                let idx = ((j % k + k) % k) as usize;
                edges[idx]
            } else {
                let ei = (((j - 1) % k + k) % k) as usize;
                (edges[ei].1, edges[ei].0)
            };
            lv_ = nlv;
            lw_ = nlw;
            j += jstep;
        }

        let bwi = childs[((j % k + k) % k) as usize];
        self.ensure(bwi);
        self.label[lw_ as usize] = 2;
        self.label[bwi as usize] = 2;
        self.labeledge[lw_ as usize] = (lv_, lw_);
        self.labeledge[bwi as usize] = (lv_, lw_);

        j += jstep;
        while childs[((j % k + k) % k) as usize] != entrychild {
            let bvi = childs[((j % k + k) % k) as usize];
            self.ensure(bvi);
            if self.label[bvi as usize] == 1 {
                j += jstep;
                continue;
            }
            let mut found_v = -1i32;
            if self.is_blossom(bvi) {
                let mut lvs = Vec::new();
                self.leaves(bvi, &mut lvs);
                for u in lvs {
                    if self.label[u as usize] != 0 { found_v = u; break; }
                }
            } else {
                found_v = bvi;
            }
            if found_v != -1 && self.label[found_v as usize] != 0 {
                self.label[found_v as usize] = 0;
                let mb = self.mate[self.blossombase[bvi as usize] as usize];
                self.label[mb as usize] = 0;
                let fv_le = self.labeledge[found_v as usize].0;
                self.assign_label(found_v, 2, fv_le);
            }
            j += jstep;
        }
    }

    // ── Augmentation through blossoms ────────────────────────────────

    fn augment_blossom(&mut self, b: i32, v: i32) {
        #[derive(Clone, Copy)]
        struct Frame { b: i32, v: i32, phase: u8, i: i32, j: i32, jstep: i32 }
        let mut stack = vec![Frame { b, v, phase: 0, i: 0, j: 0, jstep: 0 }];

        while !stack.is_empty() {
            let si = stack.len() - 1;
            match stack[si].phase {
                0 => {
                    let fv = stack[si].v;
                    let fb = stack[si].b;
                    let mut t = fv;
                    while self.blossomparent[t as usize] != fb { t = self.blossomparent[t as usize]; }
                    let k = self.blos[fb as usize].childs.len() as i32;
                    let mut i = 0i32;
                    for idx in 0..k { if self.blos[fb as usize].childs[idx as usize] == t { i = idx; break; } }
                    stack[si].i = i;
                    if self.is_blossom(t) {
                        stack[si].phase = 1;
                        stack.push(Frame { b: t, v: fv, phase: 0, i: 0, j: 0, jstep: 0 });
                        continue;
                    }
                    stack[si].phase = 2;
                    if i & 1 != 0 { stack[si].j = i - k; stack[si].jstep = 1; }
                    else           { stack[si].j = i;     stack[si].jstep = -1; }
                }
                1 => {
                    let fb = stack[si].b;
                    let fi = stack[si].i;
                    let k = self.blos[fb as usize].childs.len() as i32;
                    stack[si].phase = 2;
                    if fi & 1 != 0 { stack[si].j = fi - k; stack[si].jstep = 1; }
                    else            { stack[si].j = fi;     stack[si].jstep = -1; }
                }
                2 => {
                    let fb = stack[si].b;
                    let fi = stack[si].i;
                    let fj = stack[si].j;
                    let fv = stack[si].v;
                    let fjstep = stack[si].jstep;
                    let k = self.blos[fb as usize].childs.len() as i32;
                    if fj == 0 {
                        if fi > 0 {
                            let ii = fi as usize;
                            let bl = &mut self.blos[fb as usize];
                            let mut nc = bl.childs[ii..].to_vec();
                            nc.extend_from_slice(&bl.childs[..ii]);
                            let mut ne = bl.edges[ii..].to_vec();
                            ne.extend_from_slice(&bl.edges[..ii]);
                            bl.childs = nc;
                            bl.edges = ne;
                        }
                        self.blossombase[fb as usize] = fv;
                        stack.pop();
                        continue;
                    }
                    let nj = fj + fjstep;
                    stack[si].j = nj;
                    let idx1 = ((nj % k + k) % k) as usize;
                    let c1 = self.blos[fb as usize].childs[idx1];
                    let ww = if fjstep == 1 {
                        self.blos[fb as usize].edges[idx1].0
                    } else {
                        let ei = (((nj - 1) % k + k) % k) as usize;
                        self.blos[fb as usize].edges[ei].1
                    };
                    if self.is_blossom(c1) {
                        stack[si].phase = 3;
                        stack.push(Frame { b: c1, v: ww, phase: 0, i: 0, j: 0, jstep: 0 });
                        continue;
                    }
                    stack[si].phase = 3;
                }
                3 => {
                    let fb = stack[si].b;
                    let fj = stack[si].j;
                    let fjstep = stack[si].jstep;
                    let k = self.blos[fb as usize].childs.len() as i32;
                    let idx1 = ((fj % k + k) % k) as usize;
                    let xx = if fjstep == 1 {
                        self.blos[fb as usize].edges[idx1].1
                    } else {
                        let ei = (((fj - 1) % k + k) % k) as usize;
                        self.blos[fb as usize].edges[ei].0
                    };
                    let nj = fj + fjstep;
                    stack[si].j = nj;
                    let idx2 = ((nj % k + k) % k) as usize;
                    let c2 = self.blos[fb as usize].childs[idx2];
                    if self.is_blossom(c2) {
                        stack[si].phase = 4;
                        stack.push(Frame { b: c2, v: xx, phase: 0, i: 0, j: 0, jstep: 0 });
                        continue;
                    }
                    stack[si].phase = 4;
                }
                4 => {
                    let fb = stack[si].b;
                    let fj = stack[si].j;
                    let fjstep = stack[si].jstep;
                    let k = self.blos[fb as usize].childs.len() as i32;
                    let prev_j = fj - fjstep;
                    let (ww, xx) = if fjstep == 1 {
                        let idx1 = ((prev_j % k + k) % k) as usize;
                        self.blos[fb as usize].edges[idx1]
                    } else {
                        let ei = (((prev_j - 1) % k + k) % k) as usize;
                        (self.blos[fb as usize].edges[ei].1, self.blos[fb as usize].edges[ei].0)
                    };
                    self.mate[ww as usize] = xx;
                    self.mate[xx as usize] = ww;
                    stack[si].phase = 2;
                }
                _ => unreachable!(),
            }
        }
    }

    // ── Augmenting path ──────────────────────────────────────────────

    fn augment_path(&mut self, v: i32, w: i32) {
        let mut s = v;
        let mut j = w;
        loop {
            let bs = self.inblossom[s as usize];
            if self.is_blossom(bs) { self.augment_blossom(bs, s); }
            self.mate[s as usize] = j;
            let le = self.labeledge[bs as usize];
            if le.0 == -1 { break; }
            let t = le.0;
            let bt = self.inblossom[t as usize];
            let le2 = self.labeledge[bt as usize];
            s = le2.0;
            j = le2.1;
            if self.is_blossom(bt) { self.augment_blossom(bt, j); }
            self.mate[j as usize] = s;
        }
        self.mate[w as usize] = v;
    }

    // ── Greedy initialization ────────────────────────────────────────

    fn greedy_init(&mut self) -> i32 {
        let mut cnt = 0;
        for u in 0..self.n {
            if self.mate[u as usize] != -1 { continue; }
            for &v in &self.adj[u as usize].clone() {
                if self.mate[v as usize] == -1 {
                    self.mate[u as usize] = v;
                    self.mate[v as usize] = u;
                    cnt += 1;
                    break;
                }
            }
        }
        cnt
    }

    fn greedy_init_md(&mut self) -> i32 {
        let mut cnt = 0;
        let mut deg = vec![0i32; self.n as usize];
        for u in 0..self.n as usize {
            for &v in &self.adj[u] {
                deg[v as usize] += 1;
            }
        }
        let mut order: Vec<i32> = (0..self.n).collect();
        order.sort_unstable_by(|&a, &b| deg[a as usize].cmp(&deg[b as usize]).then(a.cmp(&b)));
        for u in order {
            if self.mate[u as usize] != -1 { continue; }
            let mut best = -1i32;
            let mut bd = i32::MAX;
            for &v in &self.adj[u as usize] {
                if self.mate[v as usize] == -1 && deg[v as usize] < bd {
                    best = v;
                    bd = deg[v as usize];
                }
            }
            if best >= 0 {
                self.mate[u as usize] = best;
                self.mate[best as usize] = u;
                cnt += 1;
            }
        }
        cnt
    }

    // ── Main solver ──────────────────────────────────────────────────

    fn solve(&mut self, greedy_mode: i32) -> Vec<(i32, i32)> {
        if greedy_mode == 1 { self.greedy_size = self.greedy_init(); }
        else if greedy_mode == 2 { self.greedy_size = self.greedy_init_md(); }

        let mut improved = true;
        while improved {
            improved = false;
            for root in 0..self.n {
                if self.mate[root as usize] != -1 { continue; }

                self.reset_blossoms();
                self.assign_label(root, 1, -1);

                let mut augmented = false;
                let mut qi = 0usize;
                while qi < self.queue.len() && !augmented {
                    let v = self.queue[qi];
                    qi += 1;
                    if self.label[self.inblossom[v as usize] as usize] != 1 { continue; }

                    let neighbors = self.adj[v as usize].clone();
                    for &w in &neighbors {
                        let bv = self.inblossom[v as usize];
                        let bw = self.inblossom[w as usize];
                        if bv == bw { continue; }
                        self.ensure(bw);

                        let lbw = self.label[bw as usize];
                        if lbw == 0 {
                            if self.mate[w as usize] == -1 {
                                self.augment_path(v, w);
                                augmented = true;
                                break;
                            }
                            self.assign_label(w, 2, v);
                        } else if lbw == 1 {
                            let base = self.scan_blossom(v, w);
                            if base >= 0 {
                                self.add_blossom(base, v, w);
                            }
                        }
                    }
                }

                // Expand all remaining blossoms
                for b in self.n..self.nblos {
                    if !self.blos[b as usize].childs.is_empty()
                        && self.blossomparent[b as usize] == -1
                    {
                        self.expand_blossom(b, true);
                    }
                }

                if augmented { improved = true; break; }
            }
        }

        let mut result = Vec::new();
        for u in 0..self.n {
            let m = self.mate[u as usize];
            if m != -1 && m > u { result.push((u, m)); }
        }
        result.sort_unstable();
        result
    }
}

// ── Validation and main ──────────────────────────────────────────────

fn validate_matching(n: i32, graph: &[Vec<i32>], matching: &[(i32, i32)]) {
    let mut deg = vec![0i32; n as usize];
    let mut errors = 0;
    for &(u, v) in matching {
        if graph[u as usize].binary_search(&v).is_err() {
            eprintln!("ERROR: Edge ({},{}) not in graph!", u, v);
            errors += 1;
        }
        deg[u as usize] += 1;
        deg[v as usize] += 1;
    }
    for i in 0..n as usize {
        if deg[i] > 1 {
            eprintln!("ERROR: Vertex {} in {} edges!", i, deg[i]);
            errors += 1;
        }
    }
    let matched = deg.iter().filter(|&&d| d > 0).count();
    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.len());
    println!("Matched vertices: {}", matched);
    println!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
    println!("=========================\n");
}

fn load_graph(filename: &str) -> Result<(i32, Vec<(i32, i32)>), Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let first = lines.next().ok_or("Empty file")??;
    let parts: Vec<&str> = first.split_whitespace().collect();
    let n: i32 = parts[0].parse()?;
    let _m: i32 = parts[1].parse()?;

    let mut edges = Vec::new();
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let u: i32 = parts[0].parse()?;
            let v: i32 = parts[1].parse()?;
            edges.push((u, v));
        }
    }
    Ok((n, edges))
}

fn main() {
    println!("Edmonds' Blossom Algorithm (Simple) - Rust Implementation");
    println!("==========================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename> [--greedy|--greedy-md]", args[0]);
        std::process::exit(1);
    }

    let mut gm = 0;
    for a in &args[2..] {
        match a.as_str() {
            "--greedy" => gm = 1,
            "--greedy-md" => gm = 2,
            _ => {}
        }
    }

    match load_graph(&args[1]) {
        Ok((n, edges)) => {
            println!("Graph: {} vertices, {} edges", n, edges.len());

            let start = Instant::now();
            let mut sol = Solver::new(n, &edges);
            let matching = sol.solve(gm);
            let duration = start.elapsed();

            validate_matching(n, &sol.adj, &matching);

            println!("Matching size: {}", matching.len());
            if gm > 0 {
                println!("Greedy init size: {}", sol.greedy_size);
                if !matching.is_empty() {
                    println!("Greedy/Final: {:.2}%", 100.0 * sol.greedy_size as f64 / matching.len() as f64);
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

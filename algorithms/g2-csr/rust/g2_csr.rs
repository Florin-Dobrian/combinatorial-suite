// Gabow MCM with proper duals and edge IDs. CSR variant.
// O(E * sqrt(V)) Maximum Cardinality Matching.
//
// Edge model: each undirected edge gets a unique ID with fixed
// (src, tgt) endpoints. Adjacency is stored as CSR (offsets + flat
// edge ID buffer), not as Vec<Vec>. Iteration sites use index-based
// loops over adj_off[v]..adj_off[v+1], avoiding the per-site .clone()
// that the Vec<Vec> version required to satisfy the borrow checker.
// The H-side structure (contracted_into) remains Vec<Vec> because it
// grows dynamically during each phase.
//
// All integers. No floating point in algorithm. No dependencies.

use std::io::{self, BufRead};
use std::time::Instant;

const NIL: i32 = -1;
const EVEN: i32 = 0;
const ODD: i32 = 1;
const UNLABELED: i32 = 2;

struct GabowMCM {
    n: usize,
    m: usize,
    greedy_size: usize,

    // Edge storage
    esrc: Vec<i32>,
    etgt: Vec<i32>,
    // adj as CSR:
    //   adj_off[v] .. adj_off[v+1] is the range of edge IDs incident to v
    //   adj_edges[j] is the edge ID at position j
    adj_edges: Vec<usize>,
    adj_off: Vec<usize>,
    mate: Vec<i32>,

    // Phase 1
    label: Vec<i32>,
    parent: Vec<i32>,
    source_bridge: Vec<i32>,
    target_bridge: Vec<i32>,
    bd: Vec<i32>,
    b_delta: Vec<i32>,

    base_par: Vec<usize>,
    dbase_par: Vec<usize>,
    dbase_rank: Vec<u32>,

    max_pq: usize,
    pq: Vec<Vec<usize>>,  // L[d] = edge IDs becoming tight at Delta=d

    path1: Vec<f64>,
    path2: Vec<f64>,
    strue: f64,

    tree_nodes: Vec<usize>,
    delta: i32,

    // H state
    rep: Vec<usize>,
    mate_h: Vec<i32>,
    is_h: Vec<bool>,

    // Phase 2
    label_h: Vec<i32>,
    parent_h: Vec<i32>,    // edge ID or NIL
    bridge_h: Vec<i32>,    // edge ID or NIL
    dir_h: Vec<i32>,       // 1 or -1
    even_time_h: Vec<i32>,
    t_h: i32,
    contracted_into: Vec<Vec<usize>>,
    size_of_m: usize,

    max_delta_used: usize,
    prev_tree_nodes: Vec<usize>,
    free_vertices: Vec<usize>,
    free_list_built: bool,
}

impl GabowMCM {
    fn new(n: usize, edges: &[(i32, i32)]) -> Self {
        let max_pq = n / 2 + 2;

        // Dedup edges
        let mut sorted: Vec<(i32, i32)> = Vec::with_capacity(edges.len());
        for &(u, v) in edges {
            if u >= 0 && (u as usize) < n && v >= 0 && (v as usize) < n && u != v {
                let (a, b) = if u < v { (u, v) } else { (v, u) };
                sorted.push((a, b));
            }
        }
        sorted.sort_unstable();
        sorted.dedup();

        let m = sorted.len();
        let mut esrc = vec![0i32; m];
        let mut etgt = vec![0i32; m];
        for (i, &(u, v)) in sorted.iter().enumerate() {
            esrc[i] = u;
            etgt[i] = v;
        }

        // Build CSR adjacency: adj_off is prefix-sum of degrees,
        // adj_edges holds edge IDs grouped by source vertex.
        let mut adj_off: Vec<usize> = vec![0; n + 1];
        for i in 0..m {
            adj_off[esrc[i] as usize + 1] += 1;
            adj_off[etgt[i] as usize + 1] += 1;
        }
        for v in 0..n {
            adj_off[v + 1] += adj_off[v];
        }
        let mut adj_edges: Vec<usize> = vec![0; adj_off[n]];
        // tmp_pos tracks next write position per vertex during the fill
        let mut tmp_pos: Vec<usize> = adj_off[0..n].to_vec();
        for i in 0..m {
            let u = esrc[i] as usize;
            let v = etgt[i] as usize;
            adj_edges[tmp_pos[u]] = i;
            tmp_pos[u] += 1;
            adj_edges[tmp_pos[v]] = i;
            tmp_pos[v] += 1;
        }

        let base_par: Vec<usize> = (0..n).collect();
        let dbase_par: Vec<usize> = (0..n).collect();

        GabowMCM {
            n, m,
            greedy_size: 0,
            esrc, etgt, adj_edges, adj_off,
            mate: vec![NIL; n],
            label: vec![UNLABELED; n],
            parent: vec![NIL; n],
            source_bridge: vec![NIL; n],
            target_bridge: vec![NIL; n],
            bd: vec![1; n],
            b_delta: vec![0; n],
            base_par,
            dbase_par,
            dbase_rank: vec![0; n],
            max_pq,
            pq: vec![Vec::new(); max_pq],
            path1: vec![0.0; n],
            path2: vec![0.0; n],
            strue: 0.0,
            tree_nodes: Vec::new(),
            delta: 0,
            rep: vec![0; n],
            mate_h: vec![NIL; n],
            is_h: vec![false; m],
            label_h: vec![UNLABELED; n],
            parent_h: vec![NIL; n],
            bridge_h: vec![NIL; n],
            dir_h: vec![0; n],
            even_time_h: vec![0; n],
            t_h: 0,
            contracted_into: vec![Vec::new(); n],
            size_of_m: 0,
            max_delta_used: 0,
            prev_tree_nodes: Vec::new(),
            free_vertices: Vec::new(),
            free_list_built: false,
        }
    }

    #[inline]
    fn opposite(&self, v: i32, eid: usize) -> i32 {
        if self.esrc[eid] == v { self.etgt[eid] } else { self.esrc[eid] }
    }

    #[inline]
    fn w_edge(&self, eid: usize) -> i32 {
        if self.mate[self.esrc[eid] as usize] == self.etgt[eid] { 2 } else { 0 }
    }

    fn find_base(&mut self, mut v: usize) -> usize {
        while self.base_par[v] != v {
            self.base_par[v] = self.base_par[self.base_par[v]];
            v = self.base_par[v];
        }
        v
    }

    fn find_dbase(&mut self, mut v: usize) -> usize {
        while self.dbase_par[v] != v {
            self.dbase_par[v] = self.dbase_par[self.dbase_par[v]];
            v = self.dbase_par[v];
        }
        v
    }

    fn union_dbase(&mut self, a: usize, b: usize) {
        let a = self.find_dbase(a);
        let b = self.find_dbase(b);
        if a == b { return; }
        if self.dbase_rank[a] < self.dbase_rank[b] {
            self.dbase_par[a] = b;
        } else {
            self.dbase_par[b] = a;
            if self.dbase_rank[a] == self.dbase_rank[b] {
                self.dbase_rank[a] += 1;
            }
        }
    }

    fn make_rep_dbase(&mut self, v: usize) {
        let r = self.find_dbase(v);
        if r != v {
            self.dbase_par[r] = v;
            self.dbase_par[v] = v;
        }
    }

    fn d(&mut self, v: usize) -> i32 {
        let bv = self.find_base(v);
        if self.label[bv] == UNLABELED { return 1; }
        if self.label[bv] == EVEN {
            return self.bd[v] - (self.delta - self.b_delta[v]);
        }
        self.bd[v] + (self.delta - self.b_delta[v]) // ODD
    }

    fn scan_edge(&mut self, eid: usize, z: usize) {
        let u = self.opposite(z as i32, eid) as usize;
        if self.mate[u] == z as i32 { return; }
        let bu = self.find_base(u);
        if self.label[bu] == ODD { return; }
        let p = self.d(z) + self.d(u);
        let tight_at = if self.label[bu] == UNLABELED {
            self.delta + p
        } else {
            self.delta + p / 2
        };
        if tight_at >= 0 && (tight_at as usize) < self.max_pq {
            let ta = tight_at as usize;
            self.pq[ta].push(eid);
            if ta > self.max_delta_used { self.max_delta_used = ta; }
        }
    }

    fn shrink_path(&mut self, b: usize, x: usize, y: usize,
                   dunions: &mut Vec<(usize, usize)>) {
        let mut v = self.find_base(x);
        while v != b {
            self.base_par[v] = b;
            dunions.push((v, b));
            let mv = self.mate[v] as usize;
            self.base_par[mv] = b;
            dunions.push((mv, b));
            self.base_par[b] = b;
            self.source_bridge[mv] = x as i32;
            self.target_bridge[mv] = y as i32;
            self.bd[mv] = self.bd[mv] + (self.delta - self.b_delta[mv]);
            self.b_delta[mv] = self.delta;
            let s = self.adj_off[mv];
            let e = self.adj_off[mv + 1];
            for j in s..e {
                let eid = self.adj_edges[j];
                self.scan_edge(eid, mv);
            }
            let pmv = self.parent[mv] as usize;
            v = self.find_base(pmv);
        }
        dunions.push((b, b));
    }

    fn build_free_list(&mut self) {
        self.free_vertices.clear();
        for v in 0..self.n {
            if self.mate[v] == NIL {
                self.free_vertices.push(v);
            }
        }
        self.free_list_built = true;
    }

    fn update_free_list(&mut self) {
        let mate = &self.mate;
        self.free_vertices.retain(|&v| mate[v] == NIL);
    }

    fn phase_1(&mut self) -> bool {
        self.delta = 0;
        self.tree_nodes.clear();
        let clear_limit = (self.max_delta_used + 1).min(self.max_pq);
        for i in 0..clear_limit {
            self.pq[i].clear();
        }
        self.max_delta_used = 0;
        let mut dunions: Vec<(usize, usize)> = Vec::new();

        // Reset only previous tree nodes
        let prev = std::mem::take(&mut self.prev_tree_nodes);
        for &v in &prev {
            self.base_par[v] = v;
            self.dbase_par[v] = v;
            self.dbase_rank[v] = 0;
            self.label[v] = UNLABELED;
            self.parent[v] = NIL;
            self.source_bridge[v] = NIL;
            self.target_bridge[v] = NIL;
            self.bd[v] = 1;
            self.b_delta[v] = 0;
            let s = self.adj_off[v];
            let e = self.adj_off[v + 1];
            for j in s..e {
                let eid = self.adj_edges[j];
                self.is_h[eid] = false;
            }
        }
        drop(prev);

        // Build or update free vertex list
        if !self.free_list_built {
            self.build_free_list();
        } else {
            self.update_free_list();
        }

        // Label free vertices EVEN, then scan
        let free: Vec<usize> = self.free_vertices.clone();
        for &v in &free {
            self.label[v] = EVEN;
            self.tree_nodes.push(v);
        }
        for &v in &free {
            let s = self.adj_off[v];
            let e = self.adj_off[v + 1];
            for j in s..e {
                let eid = self.adj_edges[j];
                self.scan_edge(eid, v);
            }
        }

        let mut found_sap = false;

        while (self.delta as usize) <= self.max_delta_used {
            // Skip empty levels
            while (self.delta as usize) <= self.max_delta_used
                && self.pq[self.delta as usize].is_empty()
            {
                self.delta += 1;
            }
            if (self.delta as usize) > self.max_delta_used { break; }

            let di = self.delta as usize;
            let mut qi = 0;
            while qi < self.pq[di].len() {
                let eid = self.pq[di][qi];
                qi += 1;
                let mut x = self.esrc[eid];
                let mut y = self.etgt[eid];

                // Stale-entry guard: the priority queue schedules edges
                // based on predicted d-trajectories. A later label change
                // on either endpoint can invalidate the prediction. Every
                // label change re-scans and enqueues a fresh correct entry,
                // so discarding stale entries here loses no correct SAP.
                let dx = self.d(x as usize);
                let dy = self.d(y as usize);
                if dx + dy != self.w_edge(eid) {
                    continue;
                }

                let bx = self.find_base(x as usize);
                if self.label[bx] != EVEN {
                    std::mem::swap(&mut x, &mut y);
                }
                let bx2 = self.find_base(x as usize);
                let by = self.find_base(y as usize);
                if y == self.mate[x as usize] || bx2 == by || self.label[by] == ODD {
                    continue;
                }

                if self.label[by] == UNLABELED {
                    let z = self.mate[y as usize];
                    self.bd[y as usize] = 1;
                    self.b_delta[y as usize] = self.delta;
                    self.bd[z as usize] = 1;
                    self.b_delta[z as usize] = self.delta;
                    self.parent[z as usize] = y;
                    self.parent[y as usize] = x;
                    self.label[y as usize] = ODD;
                    self.label[z as usize] = EVEN;
                    self.tree_nodes.push(y as usize);
                    self.tree_nodes.push(z as usize);
                    let s = self.adj_off[z as usize];
                    let e = self.adj_off[z as usize + 1];
                    for j in s..e {
                        let e2 = self.adj_edges[j];
                        self.scan_edge(e2, z as usize);
                    }
                } else if self.label[by] == EVEN {
                    self.strue += 1.0;
                    let mut hx = self.find_base(x as usize);
                    let mut hy = self.find_base(y as usize);
                    self.path1[hx] = self.strue;
                    self.path2[hy] = self.strue;
                    let mut lca: i32 = NIL;
                    loop {
                        if self.path1[hy] == self.strue { lca = hy as i32; break; }
                        if self.path2[hx] == self.strue { lca = hx as i32; break; }
                        let hxr = self.mate[hx] == NIL
                            || self.parent[self.mate[hx] as usize] == NIL;
                        let hyr = self.mate[hy] == NIL
                            || self.parent[self.mate[hy] as usize] == NIL;
                        if hxr && hyr { break; }
                        if !hxr {
                            let mhx = self.mate[hx] as usize;
                            let pmhx = self.parent[mhx] as usize;
                            hx = self.find_base(pmhx);
                            self.path1[hx] = self.strue;
                        }
                        if !hyr {
                            let mhy = self.mate[hy] as usize;
                            let pmhy = self.parent[mhy] as usize;
                            hy = self.find_base(pmhy);
                            self.path2[hy] = self.strue;
                        }
                    }
                    if lca != NIL {
                        let lca_u = lca as usize;
                        let xu = x as usize;
                        let yu = y as usize;
                        self.shrink_path(lca_u, xu, yu, &mut dunions);
                        self.shrink_path(lca_u, yu, xu, &mut dunions);
                    } else {
                        found_sap = true;
                    }
                }
            }
            self.pq[di].clear();

            if found_sap {
                // Build H
                let tn: Vec<usize> = self.tree_nodes.clone();
                for &v in &tn {
                    let db = self.find_dbase(v);
                    self.contracted_into[db].push(v);
                    self.mate_h[v] = NIL;
                }
                // Mark tight edges
                for &u in &tn {
                    let uh = self.find_dbase(u);
                    let s = self.adj_off[u];
                    let e = self.adj_off[u + 1];
                    for j in s..e {
                        let eid = self.adj_edges[j];
                        let v = self.opposite(u as i32, eid) as usize;
                        let vh = self.find_dbase(v);
                        if uh != vh {
                            let du = self.d(u);
                            let dv = self.d(v);
                            let we = self.w_edge(eid);
                            if du + dv == we {
                                self.is_h[eid] = true;
                                if we == 2 {
                                    self.mate_h[uh] = vh as i32;
                                    self.mate_h[vh] = uh as i32;
                                }
                            }
                        }
                    }
                }
                self.prev_tree_nodes = self.tree_nodes.clone();
                return true;
            }

            for &(a, b) in &dunions {
                if a == b {
                    self.make_rep_dbase(a);
                } else {
                    self.union_dbase(a, b);
                }
            }
            dunions.clear();
            self.delta += 1;
        }
        self.prev_tree_nodes = self.tree_nodes.clone();
        false
    }

    // Phase 2: find_apHG - recursive DFS on H
    fn find_ap_hg(&mut self, vh: usize) -> i32 {
        let ci: Vec<usize> = self.contracted_into[vh].clone();
        for v in ci {
            let a_beg = self.adj_off[v];
            let a_end = self.adj_off[v + 1];
            for j in a_beg..a_end {
                let eid = self.adj_edges[j];
                if !self.is_h[eid] { continue; }
                let w = self.opposite(v as i32, eid) as usize;
                let uh = self.rep[w];
                if self.mate_h[vh] == uh as i32 { continue; }

                if self.label_h[uh] == UNLABELED {
                    let muh = self.mate_h[uh];
                    if muh == NIL {
                        self.label_h[uh] = ODD;
                        self.parent_h[uh] = eid as i32;
                        return uh as i32;
                    }
                    self.label_h[uh] = ODD;
                    self.label_h[muh as usize] = EVEN;
                    self.parent_h[uh] = eid as i32;
                    self.even_time_h[muh as usize] = self.t_h;
                    self.t_h += 1;
                    let s = self.find_ap_hg(muh as usize);
                    if s != NIL { return s; }
                } else {
                    let bh = self.find_dbase(vh);
                    let mut zh = self.find_dbase(uh);
                    if self.even_time_h[bh] < self.even_time_h[zh] {
                        let mut tmp: Vec<usize> = Vec::new();
                        let mut endpoints: Vec<usize> = Vec::new();
                        while zh != bh {
                            endpoints.push(zh);
                            zh = self.mate_h[zh] as usize;
                            endpoints.push(zh);
                            tmp.insert(0, zh);
                            let pe = self.parent_h[zh] as usize;
                            let rep_src = self.rep[self.esrc[pe] as usize];
                            let next = if rep_src == zh {
                                self.etgt[pe] as usize
                            } else {
                                self.esrc[pe] as usize
                            };
                            zh = self.find_dbase(self.rep[next]);
                        }
                        for &nd in &endpoints {
                            self.union_dbase(nd, bh);
                        }
                        self.make_rep_dbase(bh);
                        for &odd_node in &tmp {
                            self.bridge_h[odd_node] = eid as i32;
                            self.dir_h[odd_node] =
                                if self.etgt[eid] == v as i32 { 1 } else { -1 };
                        }
                        for &odd_node in &tmp {
                            let s = self.find_ap_hg(odd_node);
                            if s != NIL { return s; }
                        }
                    }
                }
            }
        }
        NIL
    }

    fn find_path_in_hg(&self, path: &mut Vec<usize>, vh: usize, uh: usize) {
        if vh == uh { return; }
        if self.label_h[vh] == EVEN {
            let mvh = self.mate_h[vh] as usize;
            let pe = self.parent_h[mvh] as usize;
            path.push(pe);
            let rep_src = self.rep[self.esrc[pe] as usize];
            let next = if rep_src == mvh {
                self.etgt[pe] as usize
            } else {
                self.esrc[pe] as usize
            };
            self.find_path_in_hg(path, self.rep[next], uh);
        } else {
            // ODD: use bridge
            let be = self.bridge_h[vh] as usize;
            let (mate_side, uh_side) = if self.dir_h[vh] == 1 {
                (self.rep[self.esrc[be] as usize], self.rep[self.etgt[be] as usize])
            } else {
                (self.rep[self.etgt[be] as usize], self.rep[self.esrc[be] as usize])
            };
            let mt = if self.mate_h[vh] != NIL {
                self.rep[self.mate_h[vh] as usize]
            } else {
                vh
            };
            self.find_path_in_hg(path, mate_side, mt);
            path.push(be);
            self.find_path_in_hg(path, uh_side, uh);
        }
    }

    fn find_path_in_g(&self, pairs: &mut Vec<(usize, usize)>, v: usize, u: usize) {
        if v == u { return; }
        if self.label[v] == EVEN {
            let mv = self.mate[v] as usize;
            let pmv = self.parent[mv] as usize;
            pairs.push((mv, pmv));
            self.find_path_in_g(pairs, pmv, u);
        } else {
            let sb = self.source_bridge[v] as usize;
            let tb = self.target_bridge[v] as usize;
            let mv = self.mate[v] as usize;
            self.find_path_in_g(pairs, sb, mv);
            pairs.push((sb, tb));
            self.find_path_in_g(pairs, tb, u);
        }
    }

    fn augment_g(&mut self, h_edge_ids: &[usize]) {
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        for &eid in h_edge_ids {
            let u = self.esrc[eid] as usize;
            let v = self.etgt[eid] as usize;
            pairs.push((u, v));
            let rep_u = self.rep[u];
            let rep_v = self.rep[v];
            self.find_path_in_g(&mut pairs, u, rep_u);
            self.find_path_in_g(&mut pairs, v, rep_v);
        }
        for &(a, b) in &pairs {
            self.mate[a] = b as i32;
            self.mate[b] = a as i32;
        }
        self.size_of_m += 1;
    }

    fn phase_2(&mut self) {
        self.t_h = 0;
        // DO NOT reset dbase here — keep Phase 1 state
        let tn: Vec<usize> = self.tree_nodes.clone();
        for &v in &tn {
            self.rep[v] = self.find_dbase(v);
            self.label_h[v] = UNLABELED;
            self.parent_h[v] = NIL;
            self.bridge_h[v] = NIL;
            self.dir_h[v] = 0;
            self.even_time_h[v] = 0;
        }

        let mut all_paths: Vec<Vec<usize>> = Vec::new();

        for &vh in &tn {
            if vh != self.rep[vh] { continue; }
            if self.label_h[vh] == UNLABELED && self.mate_h[vh] == NIL {
                self.label_h[vh] = EVEN;
                self.even_time_h[vh] = self.t_h;
                self.t_h += 1;
                let found = self.find_ap_hg(vh);
                if found != NIL {
                    let found_u = found as usize;
                    let pe = self.parent_h[found_u] as usize;
                    let mut path = vec![pe];
                    let rep_src = self.rep[self.esrc[pe] as usize];
                    let next = if rep_src == found_u {
                        self.etgt[pe] as usize
                    } else {
                        self.esrc[pe] as usize
                    };
                    self.find_path_in_hg(&mut path, self.rep[next], vh);
                    all_paths.push(path);
                }
            }
        }

        for path in &all_paths {
            self.augment_g(path);
        }

        for &v in &tn {
            self.contracted_into[v].clear();
            self.mate_h[v] = NIL;
        }
    }

    fn greedy_init(&mut self) -> usize {
        let mut cnt = 0;
        for u in 0..self.n {
            if self.mate[u] != NIL { continue; }
            let s = self.adj_off[u];
            let e = self.adj_off[u + 1];
            for j in s..e {
                let eid = self.adj_edges[j];
                let v = self.opposite(u as i32, eid) as usize;
                if self.mate[v] == NIL {
                    self.mate[u] = v as i32;
                    self.mate[v] = u as i32;
                    cnt += 1;
                    break;
                }
            }
        }
        cnt
    }

    fn greedy_init_md(&mut self) -> usize {
        let mut cnt = 0;
        let mut deg = vec![0u32; self.n];
        for i in 0..self.m {
            deg[self.esrc[i] as usize] += 1;
            deg[self.etgt[i] as usize] += 1;
        }
        let mut order: Vec<usize> = (0..self.n).collect();
        order.sort_unstable_by(|&a, &b| {
            deg[a].cmp(&deg[b]).then(a.cmp(&b))
        });
        for u in order {
            if self.mate[u] != NIL { continue; }
            let mut best: i32 = NIL;
            let mut best_deg = u32::MAX;
            let s = self.adj_off[u];
            let e = self.adj_off[u + 1];
            for j in s..e {
                let eid = self.adj_edges[j];
                let v = self.opposite(u as i32, eid) as usize;
                if self.mate[v] == NIL && deg[v] < best_deg {
                    best = v as i32;
                    best_deg = deg[v];
                }
            }
            if best != NIL {
                self.mate[u] = best;
                self.mate[best as usize] = u as i32;
                cnt += 1;
            }
        }
        cnt
    }

    fn solve(&mut self, greedy_mode: i32) -> Vec<(usize, usize)> {
        if greedy_mode == 1 {
            self.greedy_size = self.greedy_init();
            self.size_of_m = self.greedy_size;
        } else if greedy_mode == 2 {
            self.greedy_size = self.greedy_init_md();
            self.size_of_m = self.greedy_size;
        }

        let mut phase_count = 0;
        loop {
            let has_sap = self.phase_1();
            if !has_sap { break; }
            self.phase_2();
            phase_count += 1;
        }
        println!("Phases: {}", phase_count);

        let mut result: Vec<(usize, usize)> = Vec::new();
        for u in 0..self.n {
            if self.mate[u] != NIL && self.mate[u] as usize > u {
                result.push((u, self.mate[u] as usize));
            }
        }
        result.sort_unstable();
        result
    }
}

fn validate_matching(n: usize, adj_edges: &[usize], adj_off: &[usize],
                     esrc: &[i32], etgt: &[i32],
                     matching: &[(usize, usize)]) {
    let mut deg = vec![0u32; n];
    let mut errors = 0;
    for &(u, v) in matching {
        let mut found = false;
        let s = adj_off[u];
        let e = adj_off[u + 1];
        for j in s..e {
            let eid = adj_edges[j];
            if (esrc[eid] == u as i32 && etgt[eid] == v as i32)
                || (esrc[eid] == v as i32 && etgt[eid] == u as i32) {
                found = true;
                break;
            }
        }
        if !found {
            eprintln!("ERROR: Edge ({},{}) not in graph!", u, v);
            errors += 1;
        }
        deg[u] += 1;
        deg[v] += 1;
    }
    for i in 0..n {
        if deg[i] > 1 {
            eprintln!("ERROR: Vertex {} in {} edges!", i, deg[i]);
            errors += 1;
        }
    }
    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.len());
    println!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
    println!("=========================\n");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    println!("Gabow MCM (duals + edge IDs) - Rust");
    println!("====================================\n");
    if args.len() < 2 {
        println!("Usage: {} <filename> [--greedy|--greedy-md]", args[0]);
        return;
    }
    let mut greedy_mode = 0i32;
    for arg in &args[2..] {
        if arg == "--greedy" { greedy_mode = 1; }
        else if arg == "--greedy-md" { greedy_mode = 2; }
    }

    let file = std::fs::File::open(&args[1]).expect("Cannot open file");
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();

    let header = lines.next().unwrap().unwrap();
    let hdr: Vec<i32> = header.split_whitespace()
        .map(|s| s.parse().unwrap()).collect();
    let n = hdr[0] as usize;
    let m_in = hdr[1] as usize;

    let mut edges: Vec<(i32, i32)> = Vec::with_capacity(m_in);
    for line in lines {
        let line = line.unwrap();
        let parts: Vec<i32> = line.split_whitespace()
            .map(|s| s.parse().unwrap()).collect();
        if parts.len() >= 2 {
            edges.push((parts[0], parts[1]));
        }
    }

    println!("Graph: {} vertices, {} edges (input)", n, edges.len());
    let t0 = Instant::now();
    let mut gabow = GabowMCM::new(n, &edges);
    println!("Graph: {} vertices, {} edges (deduped)", gabow.n, gabow.m);
    let matching = gabow.solve(greedy_mode);
    let elapsed = t0.elapsed();
    validate_matching(n, &gabow.adj_edges, &gabow.adj_off, &gabow.esrc, &gabow.etgt, &matching);
    println!("Matching size: {}", matching.len());
    if greedy_mode > 0 {
        println!("Greedy init size: {}", gabow.greedy_size);
        if !matching.is_empty() {
            println!("Greedy/Final: {:.2}%", 100.0 * gabow.greedy_size as f64 / matching.len() as f64);
        }
    }
    println!("Time: {} ms", elapsed.as_millis());
}

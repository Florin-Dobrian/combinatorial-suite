/*
 * Gabow MCM (duals + edge IDs) - O(E√V) Maximum Cardinality Matching
 *
 * CSR adjacency: contiguous flat arrays.
 * Bucket-PQ Edmonds search (phase 1) builds H minor; explicit DFS on
 * H (phase 2) extracts a maximal vertex-disjoint set of saps.
 *
 * Refactored to separate input (GeneralGraph), output (GeneralMatching),
 * and algorithm state (G2State).
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;
const EVEN: i32 = 0;
const ODD: i32 = 1;
const UNLABELED: i32 = 2;

/* ---------- Input: GeneralGraph ---------- */

#[allow(dead_code)]
struct GeneralGraph {
    num_vtxs: usize,
    num_edgs: usize,            // number of undirected edges
    idx: Vec<usize>,            // row-pointer, length num_vtxs+1
    adj: Vec<i32>,              // length 2*num_edgs, neighbor vertex IDs
}

/* ---------- Output: GeneralMatching ---------- */

#[allow(dead_code)]
struct GeneralMatching {
    mate: Vec<i32>,
    num_edgs: usize,
}

/* ---------- State: G2State ---------- */

struct G2State {
    /* Edge storage: s_vtx[eid], t_vtx[eid] are fixed endpoints */
    s_vtx: Vec<i32>,
    t_vtx: Vec<i32>,
    /* edg_adj: vertex -> edge IDs incident to it. Uses graph.idx as
     * the row pointer (graph.idx[v]..graph.idx[v+1] gives the same
     * row range, since each undirected edge contributes once to
     * each endpoint's adjacency in both views). */
    edg_adj: Vec<i32>,

    /* ---- Phase 1: alternating tree ---- */
    label: Vec<i32>,
    parent: Vec<i32>,             // parent vertex in alternating tree
    source_bridge: Vec<i32>,
    target_bridge: Vec<i32>,

    /* ---- Phase 1: dual checkpoints ---- */
    bd: Vec<i32>,
    b_delta: Vec<i32>,
    delta: i32,

    /* ---- Phase 1: blossom DSU (current) ---- */
    base_par: Vec<i32>,

    /* ---- Phase 1/2: blossom DSU (snapshotted across H-construction) ---- */
    dbase_par: Vec<i32>,
    dbase_rank: Vec<i32>,

    /* ---- Phase 1: bucket priority queue ---- */
    max_pq: i32,
    pq: Vec<Vec<i32>>,            // pq[d] = edge IDs becoming tight at delta=d
    max_delta_used: i32,

    /* ---- Phase 1: augmenting-path-vs-blossom walking trick ---- */
    path1: Vec<f64>,
    path2: Vec<f64>,
    strue: f64,

    /* ---- Phase 1: visited tracking ---- */
    tree_nodes: Vec<i32>,
    prev_tree_nodes: Vec<i32>,

    /* ---- Free vertex tracking ---- */
    free_vertices: Vec<i32>,
    free_list_built: bool,

    /* ---- H minor ---- */
    rep: Vec<i32>,
    mate_h: Vec<i32>,
    is_h: Vec<bool>,              // per edge ID: is this edge in H?
    contracted_into: Vec<Vec<i32>>,

    /* ---- Phase 2: search on H ---- */
    label_h: Vec<i32>,
    parent_h: Vec<i32>,           // parent_h[uh] = edge ID
    bridge_h: Vec<i32>,           // bridge_h[vh] = edge ID
    dir_h: Vec<i32>,              // dir_h[vh] = 1 or -1
    even_time_h: Vec<i32>,
    t_h: i32,
}

/* ---------- GeneralGraph construction ---------- */

fn build_general_graph(num_vtxs: usize, edges: &[(i32, i32)]) -> GeneralGraph {
    let mut tmp: Vec<Vec<i32>> = vec![Vec::new(); num_vtxs];
    for &(u, v) in edges {
        if u >= 0 && (u as usize) < num_vtxs && v >= 0 && (v as usize) < num_vtxs && u != v {
            tmp[u as usize].push(v);
            tmp[v as usize].push(u);
        }
    }
    for adj in &mut tmp { adj.sort_unstable(); adj.dedup(); }

    let mut idx = vec![0usize; num_vtxs + 1];
    for v in 0..num_vtxs { idx[v + 1] = idx[v] + tmp[v].len(); }
    let mut adj: Vec<i32> = Vec::with_capacity(idx[num_vtxs]);
    for v in 0..num_vtxs { adj.extend_from_slice(&tmp[v]); }

    /* Each undirected edge is stored twice (once per endpoint),
     * so the total adj length is 2 * num_edgs. */
    let num_edgs = idx[num_vtxs] / 2;

    GeneralGraph { num_vtxs, num_edgs, idx, adj }
}

/* ---------- GeneralMatching construction ---------- */

fn empty_general_matching(graph: &GeneralGraph) -> GeneralMatching {
    GeneralMatching {
        mate: vec![NIL; graph.num_vtxs],
        num_edgs: 0,
    }
}

/* ---------- G2State construction ---------- */

fn empty_g2_state(graph: &GeneralGraph) -> G2State {
    let n = graph.num_vtxs;
    let m = graph.num_edgs;

    /* Build s_vtx/t_vtx from the graph's adjacency. We pick each
     * undirected edge exactly once by iterating u < v. */
    let mut s_vtx: Vec<i32> = Vec::with_capacity(m);
    let mut t_vtx: Vec<i32> = Vec::with_capacity(m);
    for u in 0..n {
        let u_begin = graph.idx[u];
        let u_end = graph.idx[u + 1];
        for k in u_begin..u_end {
            let v = graph.adj[k];
            if (u as i32) < v {
                s_vtx.push(u as i32);
                t_vtx.push(v);
            }
        }
    }

    /* Build edg_adj (vertex -> edge IDs), reusing graph.idx as the
     * row pointer (graph.idx[v]..graph.idx[v+1] gives the same row
     * range, since each undirected edge contributes once to each
     * endpoint). */
    let mut edg_adj = vec![0i32; graph.idx[n]];
    let mut tmp_pos: Vec<usize> = graph.idx[0..n].to_vec();
    for i in 0..m {
        let u = s_vtx[i] as usize;
        let v = t_vtx[i] as usize;
        edg_adj[tmp_pos[u]] = i as i32; tmp_pos[u] += 1;
        edg_adj[tmp_pos[v]] = i as i32; tmp_pos[v] += 1;
    }

    let max_pq: i32 = (n as i32) / 2 + 2;

    G2State {
        s_vtx,
        t_vtx,
        edg_adj,

        label: vec![UNLABELED; n],
        parent: vec![NIL; n],
        source_bridge: vec![NIL; n],
        target_bridge: vec![NIL; n],

        bd: vec![1; n],
        b_delta: vec![0; n],
        delta: 0,

        base_par: (0..n as i32).collect(),

        dbase_par: (0..n as i32).collect(),
        dbase_rank: vec![0; n],

        max_pq,
        pq: vec![Vec::new(); max_pq as usize],
        max_delta_used: 0,

        path1: vec![0.0; n],
        path2: vec![0.0; n],
        strue: 0.0,

        tree_nodes: Vec::new(),
        prev_tree_nodes: Vec::new(),

        free_vertices: Vec::new(),
        free_list_built: false,

        is_h: vec![false; m],

        rep: vec![0; n],
        mate_h: vec![NIL; n],
        contracted_into: vec![Vec::new(); n],

        label_h: vec![UNLABELED; n],
        parent_h: vec![NIL; n],
        bridge_h: vec![NIL; n],
        dir_h: vec![0; n],
        even_time_h: vec![0; n],
        t_h: 0,
    }
}

/* ---------- Helpers ---------- */

#[inline]
fn opposite(state: &G2State, v: i32, eid: i32) -> i32 {
    let e = eid as usize;
    if state.s_vtx[e] == v { state.t_vtx[e] } else { state.s_vtx[e] }
}

#[inline]
fn wght(matching: &GeneralMatching, state: &G2State, eid: i32) -> i32 {
    let e = eid as usize;
    if matching.mate[state.s_vtx[e] as usize] == state.t_vtx[e] { 2 } else { 0 }
}

fn find_base(state: &mut G2State, v: i32) -> i32 {
    let mut v = v;
    while state.base_par[v as usize] != v {
        let p = state.base_par[v as usize];
        let pp = state.base_par[p as usize];
        state.base_par[v as usize] = pp;
        v = pp;
    }
    v
}

fn find_dbase(state: &mut G2State, v: i32) -> i32 {
    let mut v = v;
    while state.dbase_par[v as usize] != v {
        let p = state.dbase_par[v as usize];
        let pp = state.dbase_par[p as usize];
        state.dbase_par[v as usize] = pp;
        v = pp;
    }
    v
}

fn union_dbase(state: &mut G2State, a: i32, b: i32) {
    let a = find_dbase(state, a);
    let b = find_dbase(state, b);
    if a == b { return; }
    if state.dbase_rank[a as usize] < state.dbase_rank[b as usize] {
        state.dbase_par[a as usize] = b;
    } else {
        state.dbase_par[b as usize] = a;
        if state.dbase_rank[a as usize] == state.dbase_rank[b as usize] {
            state.dbase_rank[a as usize] += 1;
        }
    }
}

fn make_rep_dbase(state: &mut G2State, v: i32) {
    let r = find_dbase(state, v);
    if r != v {
        state.dbase_par[r as usize] = v;
        state.dbase_par[v as usize] = v;
    }
}

/* dual(v): dual variable, computed on demand from label and bd/b_delta checkpoints */
fn dual(state: &mut G2State, v: i32) -> i32 {
    let bv = find_base(state, v);
    let l = state.label[bv as usize];
    if l == UNLABELED { return 1; }
    let vu = v as usize;
    if l == EVEN {
        return state.bd[vu] - (state.delta - state.b_delta[vu]);
    }
    state.bd[vu] + (state.delta - state.b_delta[vu])
}

fn scan_edge(matching: &GeneralMatching, state: &mut G2State, eid: i32, z: i32) {
    let u = opposite(state, z, eid);
    if matching.mate[u as usize] == z { return; }
    let bu = find_base(state, u);
    if state.label[bu as usize] == ODD { return; }
    let p = dual(state, z) + dual(state, u);
    let lbu = state.label[bu as usize];
    let tight_at = if lbu == UNLABELED { state.delta + p } else { state.delta + p / 2 };
    if tight_at >= 0 && tight_at < state.max_pq {
        let ta = tight_at as usize;
        state.pq[ta].push(eid);
        if tight_at > state.max_delta_used { state.max_delta_used = tight_at; }
    }
}

/* shrink_path: contract a half-cycle into the blossom rooted at b. */
fn shrink_path(graph: &GeneralGraph, matching: &GeneralMatching, state: &mut G2State,
               b: i32, x: i32, y: i32, dunions: &mut Vec<(i32, i32)>) {
    let mut v = find_base(state, x);
    while v != b {
        state.base_par[v as usize] = b;
        dunions.push((v, b));
        v = matching.mate[v as usize];
        state.base_par[v as usize] = b;
        dunions.push((v, b));
        state.base_par[b as usize] = b;
        let vu = v as usize;
        state.source_bridge[vu] = x;
        state.target_bridge[vu] = y;
        state.bd[vu] = state.bd[vu] + (state.delta - state.b_delta[vu]);
        state.b_delta[vu] = state.delta;
        let s = graph.idx[vu];
        let e = graph.idx[vu + 1];
        for j in s..e {
            let eid = state.edg_adj[j];
            scan_edge(matching, state, eid, v);
        }
        v = find_base(state, state.parent[vu]);
    }
    dunions.push((b, b));
}

/* ---- Free-vertex list maintenance ---- */
fn build_free_list(graph: &GeneralGraph, matching: &GeneralMatching, state: &mut G2State) {
    state.free_vertices.clear();
    for v in 0..graph.num_vtxs {
        if matching.mate[v] == NIL {
            state.free_vertices.push(v as i32);
        }
    }
    state.free_list_built = true;
}

fn update_free_list(matching: &GeneralMatching, state: &mut G2State) {
    state.free_vertices.retain(|&v| matching.mate[v as usize] == NIL);
}

/* ---------- Phase 1 ---------- */

fn phase1(graph: &GeneralGraph, matching: &GeneralMatching, state: &mut G2State) -> bool {
    state.delta = 0;

    /* Move previous run's tree_nodes into prev_tree_nodes (the reset
     * list), and reset tree_nodes for the new run. The swap is
     * needed because phase2 reads tree_nodes between phase1 calls,
     * so we can't clear it at the end of phase1. */
    std::mem::swap(&mut state.prev_tree_nodes, &mut state.tree_nodes);
    state.tree_nodes.clear();

    /* Only clear pq entries up to max_delta_used from previous phase */
    let clear_limit = (state.max_delta_used + 1).min(state.max_pq);
    for i in 0..clear_limit as usize {
        state.pq[i].clear();
    }
    state.max_delta_used = 0;

    let mut dunions: Vec<(i32, i32)> = Vec::new();

    /* Reset only previous tree nodes */
    let prev = std::mem::take(&mut state.prev_tree_nodes);
    for &v in &prev {
        let vu = v as usize;
        state.base_par[vu] = v;
        state.dbase_par[vu] = v;
        state.dbase_rank[vu] = 0;
        state.label[vu] = UNLABELED;
        state.parent[vu] = NIL;
        state.source_bridge[vu] = NIL;
        state.target_bridge[vu] = NIL;
        state.bd[vu] = 1;
        state.b_delta[vu] = 0;
        let s = graph.idx[vu];
        let e = graph.idx[vu + 1];
        for j in s..e {
            let eid = state.edg_adj[j];
            state.is_h[eid as usize] = false;
        }
    }
    drop(prev);

    /* Build or update free vertex list */
    if !state.free_list_built {
        build_free_list(graph, matching, state);
    } else {
        update_free_list(matching, state);
    }

    /* Label free vertices EVEN, then scan */
    let free = std::mem::take(&mut state.free_vertices);
    for &v in &free {
        state.label[v as usize] = EVEN;
        state.tree_nodes.push(v);
    }
    for &v in &free {
        let s = graph.idx[v as usize];
        let e = graph.idx[v as usize + 1];
        for j in s..e {
            let eid = state.edg_adj[j];
            scan_edge(matching, state, eid, v);
        }
    }
    state.free_vertices = free;

    let mut found_sap = false;

    while state.delta <= state.max_delta_used {
        /* Skip empty levels */
        while state.delta <= state.max_delta_used
              && state.pq[state.delta as usize].is_empty()
        {
            state.delta += 1;
        }
        if state.delta > state.max_delta_used { break; }

        let di = state.delta as usize;
        let mut qi = 0;
        while qi < state.pq[di].len() {
            let eid = state.pq[di][qi];
            qi += 1;
            let mut x = state.s_vtx[eid as usize];
            let mut y = state.t_vtx[eid as usize];

            /* Stale-entry guard: discarding stale entries here loses
             * no correct sap because every label change re-enqueues
             * a fresh correct entry. */
            let dx = dual(state, x);
            let dy = dual(state, y);
            if dx + dy != wght(matching, state, eid) { continue; }

            let bx = find_base(state, x);
            if state.label[bx as usize] != EVEN {
                std::mem::swap(&mut x, &mut y);
            }
            let bx2 = find_base(state, x);
            let by = find_base(state, y);
            if y == matching.mate[x as usize] || bx2 == by || state.label[by as usize] == ODD {
                continue;
            }

            if state.label[by as usize] == UNLABELED {
                let z = matching.mate[y as usize];
                let yu = y as usize;
                let zu = z as usize;
                state.bd[yu] = 1;
                state.b_delta[yu] = state.delta;
                state.bd[zu] = 1;
                state.b_delta[zu] = state.delta;
                state.parent[zu] = y;
                state.parent[yu] = x;
                state.label[yu] = ODD;
                state.label[zu] = EVEN;
                state.tree_nodes.push(y);
                state.tree_nodes.push(z);
                let s = graph.idx[zu];
                let e = graph.idx[zu + 1];
                for j in s..e {
                    let e2 = state.edg_adj[j];
                    scan_edge(matching, state, e2, z);
                }
            } else if state.label[by as usize] == EVEN {
                state.strue += 1.0;
                let mut hx = find_base(state, x);
                let mut hy = find_base(state, y);
                state.path1[hx as usize] = state.strue;
                state.path2[hy as usize] = state.strue;
                let mut lca: i32 = NIL;
                loop {
                    if state.path1[hy as usize] == state.strue { lca = hy; break; }
                    if state.path2[hx as usize] == state.strue { lca = hx; break; }
                    let mhx = matching.mate[hx as usize];
                    let mhy = matching.mate[hy as usize];
                    let hxr = mhx == NIL || state.parent[mhx as usize] == NIL;
                    let hyr = mhy == NIL || state.parent[mhy as usize] == NIL;
                    if hxr && hyr { break; }
                    if !hxr {
                        let pmhx = state.parent[mhx as usize];
                        hx = find_base(state, pmhx);
                        state.path1[hx as usize] = state.strue;
                    }
                    if !hyr {
                        let pmhy = state.parent[mhy as usize];
                        hy = find_base(state, pmhy);
                        state.path2[hy as usize] = state.strue;
                    }
                }
                if lca != NIL {
                    shrink_path(graph, matching, state, lca, x, y, &mut dunions);
                    shrink_path(graph, matching, state, lca, y, x, &mut dunions);
                } else {
                    found_sap = true;
                }
            }
        }
        state.pq[di].clear();

        if found_sap {
            /* Build H */
            let tn = std::mem::take(&mut state.tree_nodes);
            for &v in &tn {
                let db = find_dbase(state, v);
                state.contracted_into[db as usize].push(v);
                state.mate_h[v as usize] = NIL;
            }
            /* Mark tight edges */
            for &u in &tn {
                let uh = find_dbase(state, u);
                let s = graph.idx[u as usize];
                let e = graph.idx[u as usize + 1];
                for j in s..e {
                    let eid = state.edg_adj[j];
                    let v = opposite(state, u, eid);
                    let vh = find_dbase(state, v);
                    if uh != vh {
                        let du = dual(state, u);
                        let dv = dual(state, v);
                        let we = wght(matching, state, eid);
                        if du + dv == we {
                            state.is_h[eid as usize] = true;
                            if we == 2 {
                                state.mate_h[uh as usize] = vh;
                                state.mate_h[vh as usize] = uh;
                            }
                        }
                    }
                }
            }
            state.tree_nodes = tn;   /* restore for phase2 */
            return true;
        }

        for &(a, b) in &dunions {
            if a == b {
                make_rep_dbase(state, a);
            } else {
                union_dbase(state, a, b);
            }
        }
        dunions.clear();
        state.delta += 1;
    }
    false
}

/* ---------- Phase 2 ---------- */

/* find_ap_hg: recursive DFS on H. */
fn find_ap_hg(graph: &GeneralGraph, matching: &GeneralMatching, state: &mut G2State,
              vh: i32) -> i32 {
    let ci = state.contracted_into[vh as usize].clone();
    for v in ci {
        let a_beg = graph.idx[v as usize];
        let a_end = graph.idx[v as usize + 1];
        for j in a_beg..a_end {
            let eid = state.edg_adj[j];
            if !state.is_h[eid as usize] { continue; }
            let w = opposite(state, v, eid);
            let uh = state.rep[w as usize];
            if state.mate_h[vh as usize] == uh { continue; }

            if state.label_h[uh as usize] == UNLABELED {
                let muh = state.mate_h[uh as usize];
                if muh == NIL {
                    state.label_h[uh as usize] = ODD;
                    state.parent_h[uh as usize] = eid;
                    return uh;
                }
                state.label_h[uh as usize] = ODD;
                state.label_h[muh as usize] = EVEN;
                state.parent_h[uh as usize] = eid;
                state.even_time_h[muh as usize] = state.t_h;
                state.t_h += 1;
                let s = find_ap_hg(graph, matching, state, muh);
                if s != NIL { return s; }
            } else {
                let bh = find_dbase(state, vh);
                let mut zh = find_dbase(state, uh);
                if state.even_time_h[bh as usize] < state.even_time_h[zh as usize] {
                    let mut tmp: Vec<i32> = Vec::new();
                    let mut endpoints: Vec<i32> = Vec::new();
                    while zh != bh {
                        endpoints.push(zh);
                        zh = state.mate_h[zh as usize];
                        endpoints.push(zh);
                        tmp.insert(0, zh);
                        let pe = state.parent_h[zh as usize];
                        let rep_src = state.rep[state.s_vtx[pe as usize] as usize];
                        let next = if rep_src == zh {
                            state.t_vtx[pe as usize]
                        } else {
                            state.s_vtx[pe as usize]
                        };
                        zh = find_dbase(state, state.rep[next as usize]);
                    }
                    for &nd in &endpoints { union_dbase(state, nd, bh); }
                    make_rep_dbase(state, bh);
                    for &odd_node in &tmp {
                        state.bridge_h[odd_node as usize] = eid;
                        state.dir_h[odd_node as usize] =
                            if state.t_vtx[eid as usize] == v { 1 } else { -1 };
                    }
                    for &odd_node in &tmp {
                        let s = find_ap_hg(graph, matching, state, odd_node);
                        if s != NIL { return s; }
                    }
                }
            }
        }
    }
    NIL
}

/* find_path_in_hg: trace augmenting path in H from vh to uh. */
fn find_path_in_hg(state: &G2State, path: &mut Vec<i32>, vh: i32, uh: i32) {
    if vh == uh { return; }
    if state.label_h[vh as usize] == EVEN {
        let mvh = state.mate_h[vh as usize];
        let pe = state.parent_h[mvh as usize];
        path.push(pe);
        let rep_src = state.rep[state.s_vtx[pe as usize] as usize];
        let next = if rep_src == mvh {
            state.t_vtx[pe as usize]
        } else {
            state.s_vtx[pe as usize]
        };
        find_path_in_hg(state, path, state.rep[next as usize], uh);
    } else {
        let be = state.bridge_h[vh as usize];
        let (mate_side, uh_side) = if state.dir_h[vh as usize] == 1 {
            (state.rep[state.s_vtx[be as usize] as usize],
             state.rep[state.t_vtx[be as usize] as usize])
        } else {
            (state.rep[state.t_vtx[be as usize] as usize],
             state.rep[state.s_vtx[be as usize] as usize])
        };
        let mt = if state.mate_h[vh as usize] != NIL {
            state.rep[state.mate_h[vh as usize] as usize]
        } else {
            vh
        };
        find_path_in_hg(state, path, mate_side, mt);
        path.push(be);
        find_path_in_hg(state, path, uh_side, uh);
    }
}

/* find_path_in_g: unfold within Phase 1 blossom. */
fn find_path_in_g(matching: &GeneralMatching, state: &G2State,
                  pairs: &mut Vec<(i32, i32)>, v: i32, u: i32) {
    if v == u { return; }
    if state.label[v as usize] == EVEN {
        let mv = matching.mate[v as usize];
        let pmv = state.parent[mv as usize];
        pairs.push((mv, pmv));
        find_path_in_g(matching, state, pairs, pmv, u);
    } else {
        let sb = state.source_bridge[v as usize];
        let tb = state.target_bridge[v as usize];
        let mv = matching.mate[v as usize];
        find_path_in_g(matching, state, pairs, sb, mv);
        pairs.push((sb, tb));
        find_path_in_g(matching, state, pairs, tb, u);
    }
}

/* augment_g: unfold H-path edges to G, augment matching. */
fn augment_g(matching: &mut GeneralMatching, state: &G2State, h_edge_ids: &[i32]) {
    let mut pairs: Vec<(i32, i32)> = Vec::new();
    for &eid in h_edge_ids {
        let u = state.s_vtx[eid as usize];
        let v = state.t_vtx[eid as usize];
        pairs.push((u, v));
        find_path_in_g(matching, state, &mut pairs, u, state.rep[u as usize]);
        find_path_in_g(matching, state, &mut pairs, v, state.rep[v as usize]);
    }
    for &(a, b) in &pairs {
        matching.mate[a as usize] = b;
        matching.mate[b as usize] = a;
    }
    matching.num_edgs += 1;
}

fn phase2(graph: &GeneralGraph, matching: &mut GeneralMatching, state: &mut G2State) {
    state.t_h = 0;
    let tn = std::mem::take(&mut state.tree_nodes);
    for &v in &tn {
        let db = find_dbase(state, v);
        state.rep[v as usize] = db;
        state.label_h[v as usize] = UNLABELED;
        state.parent_h[v as usize] = NIL;
        state.bridge_h[v as usize] = NIL;
        state.dir_h[v as usize] = 0;
        state.even_time_h[v as usize] = 0;
    }

    let mut all_paths: Vec<Vec<i32>> = Vec::new();

    for &vh in &tn {
        if vh != state.rep[vh as usize] { continue; }
        if state.label_h[vh as usize] == UNLABELED && state.mate_h[vh as usize] == NIL {
            state.label_h[vh as usize] = EVEN;
            state.even_time_h[vh as usize] = state.t_h;
            state.t_h += 1;
            let found = find_ap_hg(graph, matching, state, vh);
            if found != NIL {
                let pe = state.parent_h[found as usize];
                let mut path = vec![pe];
                let rep_src = state.rep[state.s_vtx[pe as usize] as usize];
                let next = if rep_src == found {
                    state.t_vtx[pe as usize]
                } else {
                    state.s_vtx[pe as usize]
                };
                find_path_in_hg(state, &mut path, state.rep[next as usize], vh);
                all_paths.push(path);
            }
        }
    }

    for path in &all_paths {
        augment_g(matching, state, path);
    }

    for &v in &tn {
        state.contracted_into[v as usize].clear();
        state.mate_h[v as usize] = NIL;
    }
    state.tree_nodes = tn;
}

/* ---------- Greedy initial matching: simple ---------- */

fn greedy_init(graph: &GeneralGraph, matching: &mut GeneralMatching) -> usize {
    let mut num_edgs: usize = 0;
    for u in 0..graph.num_vtxs {
        if matching.mate[u] != NIL { continue; }
        let u_begin = graph.idx[u];
        let u_end = graph.idx[u + 1];
        for k in u_begin..u_end {
            let v = graph.adj[k];
            if matching.mate[v as usize] == NIL {
                matching.mate[u] = v;
                matching.mate[v as usize] = u as i32;
                num_edgs += 1;
                break;
            }
        }
    }
    matching.num_edgs += num_edgs;
    num_edgs
}

/* ---------- Greedy initial matching: min-degree ---------- */

fn greedy_init_md(graph: &GeneralGraph, matching: &mut GeneralMatching) -> usize {
    let mut num_edgs: usize = 0;
    let mut deg = vec![0usize; graph.num_vtxs];
    for u in 0..graph.num_vtxs {
        deg[u] = graph.idx[u + 1] - graph.idx[u];
    }
    let mut order: Vec<usize> = (0..graph.num_vtxs).collect();
    /* Sort vertices in increasing order of degree, breaking ties by vertex label. */
    order.sort_unstable_by(|&u1, &u2| {
        deg[u1].cmp(&deg[u2]).then(u1.cmp(&u2))
    });
    for u in order {
        if matching.mate[u] != NIL { continue; }
        let mut best: i32 = NIL;
        let mut best_deg = usize::MAX;
        let u_begin = graph.idx[u];
        let u_end = graph.idx[u + 1];
        for k in u_begin..u_end {
            let v = graph.adj[k];
            if matching.mate[v as usize] == NIL && deg[v as usize] < best_deg {
                best = v;
                best_deg = deg[v as usize];
            }
        }
        if best != NIL {
            matching.mate[u] = best;
            matching.mate[best as usize] = u as i32;
            num_edgs += 1;
        }
    }
    matching.num_edgs += num_edgs;
    num_edgs
}

/* ---------- Top-level driver ---------- */

fn g2_mcm(graph: &GeneralGraph, matching: &mut GeneralMatching) -> i32 {
    let mut state = empty_g2_state(graph);
    let mut num_phases: i32 = 0;
    loop {
        let has_sap = phase1(graph, matching, &mut state);
        if !has_sap { break; }
        phase2(graph, matching, &mut state);
        num_phases += 1;
    }
    num_phases
}

/* ---------- Validation ---------- */

fn validate_general_matching(graph: &GeneralGraph, matching: &GeneralMatching) {
    let mut errors = 0;
    let mut num_matched = 0;

    for u in 0..graph.num_vtxs {
        let v = matching.mate[u];
        if v == NIL { continue; }
        num_matched += 1;
        if v < 0 || (v as usize) >= graph.num_vtxs {
            eprintln!("ERROR: mate[{}] = {} out of range", u, v);
            errors += 1;
        } else if matching.mate[v as usize] != u as i32 {
            eprintln!("ERROR: mate[{}]={} but mate[{}]={}", u, v, v, matching.mate[v as usize]);
            errors += 1;
        } else if (u as i32) < v {
            /* Check edge presence in graph (each undirected edge once, when u < v). */
            let u_begin = graph.idx[u];
            let u_end = graph.idx[u + 1];
            if !graph.adj[u_begin..u_end].binary_search(&v).is_ok() {
                eprintln!("ERROR: edge ({},{}) not in graph", u, v);
                errors += 1;
            }
        }
    }

    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.num_edgs);
    println!("Vertices matched: {}", num_matched);
    println!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
    println!("=========================\n");
}

/* ---------- Graph loader ---------- */

fn load_graph(filename: &str) -> Result<(usize, Vec<(i32, i32)>), Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let first = lines.next().ok_or("Empty file")??;
    let parts: Vec<&str> = first.split_whitespace().collect();
    if parts.len() != 2 {
        return Err("First line must have 2 numbers".into());
    }
    let num_vtxs: usize = parts[0].parse()?;
    let m: usize = parts[1].parse()?;

    let mut edges = Vec::with_capacity(m);
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let u: i32 = parts[0].parse()?;
            let v: i32 = parts[1].parse()?;
            edges.push((u, v));
        }
    }
    Ok((num_vtxs, edges))
}

/* ---------- Main ---------- */

fn main() {
    println!("Gabow MCM (duals + edge IDs) - Rust Implementation (CSR)");
    println!("=========================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename> [--greedy|--greedy-md]", args[0]);
        std::process::exit(1);
    }
    let greedy_mode: i32 = if args.iter().any(|a| a == "--greedy-md") { 2 }
                          else if args.iter().any(|a| a == "--greedy") { 1 }
                          else { 0 };

    match load_graph(&args[1]) {
        Ok((num_vtxs, edges)) => {
            println!("Graph: {} vertices, {} edges", num_vtxs, edges.len());

            let graph = build_general_graph(num_vtxs, &edges);
            let mut matching = empty_general_matching(&graph);

            let start = Instant::now();

            let greedy_size: usize = match greedy_mode {
                1 => greedy_init(&graph, &mut matching),
                2 => greedy_init_md(&graph, &mut matching),
                _ => 0,
            };

            let num_phases = g2_mcm(&graph, &mut matching);

            let duration = start.elapsed();

            validate_general_matching(&graph, &matching);

            println!("Phases: {}", num_phases);
            println!("Matching size: {}", matching.num_edgs);
            if greedy_mode > 0 {
                let fs = matching.num_edgs;
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

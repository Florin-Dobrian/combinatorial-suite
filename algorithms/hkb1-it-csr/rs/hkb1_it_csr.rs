/*
 * Hopcroft-Karp Bidirectional, Iterative (HKB1) - O(E sqrt(V)) Maximum Bipartite Matching
 *
 * Variant 1 ("BG2"): bidirectional bucket-PQ keyed by BFS level. No FIFO.
 * Free vertices on BOTH partitions are seeded as EVEN tree roots at level 0,
 * with S-roots and T-roots growing distinct trees. Edges are scheduled into
 * buckets at the level where they become "ready":
 *
 *   EVEN-UNLABELED edge (will GROW):  bucket[level(EVEN) + 1].
 *   EVEN-EVEN cross-tree (a BRIDGE):  bucket[level(u) + level(v)].
 *
 * A bridge is recorded only when tight at its bucket (sLvl + tLvl == d), and
 * the augment crosses a bridge only at the committed level, so each phase
 * fixes one SAP length and augments a maximal set of shortest APs (MSAP).
 *
 * CSR adjacency, three-object architecture (BipartiteGraph / BipartiteMatching
 * / Hkb1State), top-level returns num_phases.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;

const LBL_UNLABELED: i8 = 0;
const LBL_EVEN: i8 = 1;
const LBL_ODD: i8 = 2;

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

/* ---------- State: Hkb1State ---------- */

struct Hkb1State {
    s_label: Vec<i8>,
    t_label: Vec<i8>,
    s_level: Vec<i32>,
    t_level: Vec<i32>,
    s_visited: Vec<bool>,
    t_visited: Vec<bool>,

    buckets: Vec<Vec<(i32, i32)>>, // buckets[d]: (s, t) edges ready at level d
    bridges: Vec<(i32, i32)>,      // bridges recorded in the committed bucket
    committed_d: i32,              // bucket level the bridges fired at
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

/* ---------- Greedy initial matching: simple ---------- */

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

/* ---------- Greedy initial matching: min-degree ---------- */

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

/* ---------- HKB1 search: bidirectional bucket-PQ by level ---------- */

fn schedule_s(s: i32, g: &BipartiteGraph, m: &BipartiteMatching, st: &mut Hkb1State) {
    /* Schedule outgoing unmatched edges of a newly-EVEN S-vertex `s`. */
    let level_s = st.s_level[s as usize];
    let mated_t = m.s_mate[s as usize];   /* NIL for free S */
    let begin = g.s_idx[s as usize];
    let end = g.s_idx[s as usize + 1];
    for k in begin..end {
        let t = g.s_adj[k];
        if t == mated_t { continue; }     /* skip matched edge */
        let t_lbl = st.t_label[t as usize];
        if t_lbl == LBL_ODD { continue; } /* already in tree on T-side */
        let target = if t_lbl == LBL_EVEN {
            level_s + st.t_level[t as usize]   /* EVEN-EVEN cross-tree bridge */
        } else {
            level_s + 1                        /* EVEN-UNLABELED GROW */
        } as usize;
        while st.buckets.len() <= target { st.buckets.push(Vec::new()); }
        st.buckets[target].push((s, t));
    }
}

fn schedule_t(t: i32, g: &BipartiteGraph, m: &BipartiteMatching, st: &mut Hkb1State) {
    /* Schedule outgoing unmatched edges of a newly-EVEN T-vertex `t`. */
    let level_t = st.t_level[t as usize];
    let mated_s = m.t_mate[t as usize];   /* NIL for free T */
    let begin = g.t_idx[t as usize];
    let end = g.t_idx[t as usize + 1];
    for k in begin..end {
        let s = g.t_adj[k];
        if s == mated_s { continue; }
        let s_lbl = st.s_label[s as usize];
        if s_lbl == LBL_ODD { continue; }
        let target = if s_lbl == LBL_EVEN {
            st.s_level[s as usize] + level_t
        } else {
            level_t + 1
        } as usize;
        while st.buckets.len() <= target { st.buckets.push(Vec::new()); }
        st.buckets[target].push((s, t));
    }
}

fn search(g: &BipartiteGraph, m: &BipartiteMatching, st: &mut Hkb1State) -> bool {
    /* Reset per-iteration state. */
    st.s_label.fill(LBL_UNLABELED);
    st.t_label.fill(LBL_UNLABELED);
    st.s_level.fill(NIL);
    st.t_level.fill(NIL);
    st.buckets.clear();
    st.bridges.clear();

    /* SEED: free S and free T as EVEN level 0. Label first, schedule second. */
    for s in 0..g.s_num_vtxs {
        if m.s_mate[s] == NIL { st.s_label[s] = LBL_EVEN; st.s_level[s] = 0; }
    }
    for t in 0..g.t_num_vtxs {
        if m.t_mate[t] == NIL { st.t_label[t] = LBL_EVEN; st.t_level[t] = 0; }
    }
    for s in 0..g.s_num_vtxs {
        if m.s_mate[s] == NIL { schedule_s(s as i32, g, m, st); }
    }
    for t in 0..g.t_num_vtxs {
        if m.t_mate[t] == NIL { schedule_t(t as i32, g, m, st); }
    }

    /* Drain buckets in level-increasing order; stop at first level with bridges. */
    let mut d: usize = 0;
    while d < st.buckets.len() {
        let mut found_this_level = false;
        let mut qi = 0;
        while qi < st.buckets[d].len() {
            let (s, t) = st.buckets[d][qi];
            qi += 1;
            let s_lbl = st.s_label[s as usize];
            let t_lbl = st.t_label[t as usize];

            /* Doubly-booked or partner-already-in-tree: skip. */
            if s_lbl == LBL_ODD || t_lbl == LBL_ODD { continue; }

            if s_lbl == LBL_EVEN && t_lbl == LBL_EVEN {
                /* Bridge candidate; record only when tight at this level. */
                if st.s_level[s as usize] + st.t_level[t as usize] == d as i32 {
                    st.bridges.push((s, t));
                    found_this_level = true;
                }
                continue;
            }

            /* GROW: one EVEN, one UNLABELED (which has a mate). */
            if s_lbl == LBL_EVEN {
                let mate_t = m.t_mate[t as usize];
                if mate_t == NIL { continue; }   /* defensive */
                st.t_label[t as usize] = LBL_ODD;
                st.t_level[t as usize] = d as i32;
                st.s_label[mate_t as usize] = LBL_EVEN;
                st.s_level[mate_t as usize] = d as i32;
                schedule_s(mate_t, g, m, st);
            } else {
                let mate_s = m.s_mate[s as usize];
                if mate_s == NIL { continue; }
                st.s_label[s as usize] = LBL_ODD;
                st.s_level[s as usize] = d as i32;
                st.t_label[mate_s as usize] = LBL_EVEN;
                st.t_level[mate_s as usize] = d as i32;
                schedule_t(mate_s, g, m, st);
            }
        }
        if found_this_level { st.committed_d = d as i32; break; }
        d += 1;
    }

    !st.bridges.is_empty()
}

/* ---------- Augment: HK-style forward DFS from each free S ---------- */

fn dfs_ap(s: i32, g: &BipartiteGraph, m: &mut BipartiteMatching, st: &mut Hkb1State) -> bool {
    if st.s_visited[s as usize] { return false; }
    st.s_visited[s as usize] = true;

    let s_lvl = st.s_level[s as usize];
    let s_lbl = st.s_label[s as usize];
    let begin = g.s_idx[s as usize];
    let end = g.s_idx[s as usize + 1];

    for k in begin..end {
        let t = g.s_adj[k];
        if t == m.s_mate[s as usize] { continue; }   /* skip matched edge */
        if st.t_visited[t as usize] { continue; }

        let t_lbl = st.t_label[t as usize];
        let t_lvl = st.t_level[t as usize];
        let mut admit = false;

        if s_lbl == LBL_EVEN {
            /* S-tree side. */
            if t_lbl == LBL_ODD && t_lvl == s_lvl + 1 {
                admit = true;                                   /* forward in S-tree */
            } else if t_lbl == LBL_EVEN && s_lvl + t_lvl == st.committed_d {
                admit = true;                                   /* bridge at committed level */
            }
        } else if s_lbl == LBL_ODD {
            /* T-tree side (post-bridge). */
            if t_lbl == LBL_EVEN && t_lvl == s_lvl - 1 {
                admit = true;
            }
        }
        if !admit { continue; }

        st.t_visited[t as usize] = true;
        let next_s = m.t_mate[t as usize];

        if next_s == NIL {
            /* t is free T -- AP found. */
            m.s_mate[s as usize] = t;
            m.t_mate[t as usize] = s;
            return true;
        }
        if dfs_ap(next_s, g, m, st) {
            m.s_mate[s as usize] = t;
            m.t_mate[t as usize] = s;
            return true;
        }
    }
    false
}

fn augment(g: &BipartiteGraph, m: &mut BipartiteMatching, st: &mut Hkb1State) -> i32 {
    st.s_visited.fill(false);
    st.t_visited.fill(false);

    let mut new_edgs = 0;
    for s in 0..g.s_num_vtxs {
        if m.s_mate[s] != NIL { continue; }            /* not free */
        if st.s_visited[s] { continue; }
        if st.s_label[s] != LBL_EVEN { continue; }     /* must be S-tree EVEN at sLvl 0 */
        if dfs_ap(s as i32, g, m, st) {
            new_edgs += 1;
        }
    }
    new_edgs
}

/* ---------- Top-level HKB1 ---------- */

fn hkb1_iterative_mcm(g: &BipartiteGraph, m: &mut BipartiteMatching) -> i32 {
    let mut st = Hkb1State {
        s_label: vec![LBL_UNLABELED; g.s_num_vtxs],
        t_label: vec![LBL_UNLABELED; g.t_num_vtxs],
        s_level: vec![NIL; g.s_num_vtxs],
        t_level: vec![NIL; g.t_num_vtxs],
        s_visited: vec![false; g.s_num_vtxs],
        t_visited: vec![false; g.t_num_vtxs],
        buckets: Vec::new(),
        bridges: Vec::new(),
        committed_d: 0,
    };

    let mut num_phases: i32 = 0;
    let mut new_edgs: i32 = 0;
    while search(g, m, &mut st) {
        num_phases += 1;
        new_edgs += augment(g, m, &mut st);
    }
    m.num_edgs += new_edgs as usize;
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
    println!("Hopcroft-Karp Bidirectional Iterative (HKB1) - Rust Implementation (CSR)");
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

            let num_phases = hkb1_iterative_mcm(&g, &mut m);

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

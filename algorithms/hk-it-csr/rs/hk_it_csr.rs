/*
 * Hopcroft-Karp Iterative Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * CSR adjacency: contiguous flat arrays.
 * Old HK's lean BFS (single s_level[] array, sentinel trick) + iterative
 * stack-based DFS with edge index array (no recursion, no rescan).
 *
 * Refactored to separate input (BipartiteGraph), output (BipartiteMatching),
 * and algorithm state (HKIState).
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;

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

/* ---------- State: HKIState ---------- */

struct HKIState {
    s_level: Vec<i32>,         // length s_num_vtxs+1; s_level[s_num_vtxs] is the NIL sentinel
    s_idx: Vec<usize>,         // length s_num_vtxs; relative offset within s's adjacency, persistent within a phase
    s_prcb_stk: Vec<i32>,      // length s_num_vtxs; DFS stack (s-vertices)
    t_prcb_stk: Vec<i32>,      // length s_num_vtxs; t-vertex chosen at each depth
    stk_top: usize,
}

/* ---------- BipartiteGraph construction ---------- */

fn build_bipartite_graph(s_num_vtxs: usize, t_num_vtxs: usize,
                         edges: &[(usize, usize)]) -> BipartiteGraph {
    let mut s_tmp: Vec<Vec<i32>> = vec![Vec::new(); s_num_vtxs];
    let mut t_tmp: Vec<Vec<i32>> = vec![Vec::new(); t_num_vtxs];
    for &(s, t) in edges {
        if s < s_num_vtxs && t < t_num_vtxs {
            s_tmp[s].push(t as i32);
            t_tmp[t].push(s as i32);
        }
    }
    for adj in &mut s_tmp { adj.sort_unstable(); adj.dedup(); }
    for adj in &mut t_tmp { adj.sort_unstable(); adj.dedup(); }

    let mut s_idx = vec![0usize; s_num_vtxs + 1];
    for s in 0..s_num_vtxs { s_idx[s + 1] = s_idx[s] + s_tmp[s].len(); }
    let mut s_adj: Vec<i32> = Vec::with_capacity(s_idx[s_num_vtxs]);
    for s in 0..s_num_vtxs {
        s_adj.extend_from_slice(&s_tmp[s]);
    }

    let mut t_idx = vec![0usize; t_num_vtxs + 1];
    for t in 0..t_num_vtxs { t_idx[t + 1] = t_idx[t] + t_tmp[t].len(); }
    let mut t_adj: Vec<i32> = Vec::with_capacity(t_idx[t_num_vtxs]);
    for t in 0..t_num_vtxs {
        t_adj.extend_from_slice(&t_tmp[t]);
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

fn empty_bipartite_matching(graph: &BipartiteGraph) -> BipartiteMatching {
    BipartiteMatching {
        s_num_vtxs: graph.s_num_vtxs,
        t_num_vtxs: graph.t_num_vtxs,
        num_edgs: 0,
        s_mate: vec![NIL; graph.s_num_vtxs],
        t_mate: vec![NIL; graph.t_num_vtxs],
    }
}

/* ---------- Greedy initial matching: simple ---------- */

fn greedy_init(graph: &BipartiteGraph, matching: &mut BipartiteMatching) -> usize {
    let mut num_edgs: usize = 0;
    for s in 0..graph.s_num_vtxs {
        if matching.s_mate[s] != NIL { continue; }
        let s_begin = graph.s_idx[s];
        let s_end = graph.s_idx[s + 1];
        for k in s_begin..s_end {
            let t = graph.s_adj[k];
            if matching.t_mate[t as usize] == NIL {
                matching.s_mate[s] = t;
                matching.t_mate[t as usize] = s as i32;
                num_edgs += 1;
                break;
            }
        }
    }
    matching.num_edgs += num_edgs;
    num_edgs
}

/* ---------- Greedy initial matching: min-degree ---------- */

fn greedy_init_md(graph: &BipartiteGraph, matching: &mut BipartiteMatching) -> usize {
    let mut num_edgs: usize = 0;
    let mut deg = vec![0usize; graph.t_num_vtxs];
    for s in 0..graph.s_num_vtxs {
        let s_begin = graph.s_idx[s];
        let s_end = graph.s_idx[s + 1];
        for k in s_begin..s_end { deg[graph.s_adj[k] as usize] += 1; }
    }
    let mut s_order: Vec<usize> = (0..graph.s_num_vtxs).collect();
    /* Sort s-vertices in increasing order of degree, breaking ties by vertex label. */
    s_order.sort_unstable_by(|&s1, &s2| {
        let s1_deg = graph.s_idx[s1 + 1] - graph.s_idx[s1];
        let s2_deg = graph.s_idx[s2 + 1] - graph.s_idx[s2];
        s1_deg.cmp(&s2_deg).then(s1.cmp(&s2))
    });
    for s in s_order {
        if matching.s_mate[s] != NIL { continue; }
        let mut best: i32 = NIL;
        let mut best_deg = usize::MAX;
        let s_begin = graph.s_idx[s];
        let s_end = graph.s_idx[s + 1];
        for k in s_begin..s_end {
            let t = graph.s_adj[k];
            if matching.t_mate[t as usize] == NIL && deg[t as usize] < best_deg {
                best = t;
                best_deg = deg[t as usize];
            }
        }
        if best != NIL {
            matching.s_mate[s] = best;
            matching.t_mate[best as usize] = s as i32;
            num_edgs += 1;
        }
    }
    matching.num_edgs += num_edgs;
    num_edgs
}

/* ---------- HK BFS ---------- */

fn bfs(graph: &BipartiteGraph, matching: &BipartiteMatching, state: &mut HKIState) -> bool {
    let mut s_prcb_que: Vec<usize> = Vec::with_capacity(graph.s_num_vtxs);

    for s in 0..graph.s_num_vtxs {
        if matching.s_mate[s] == NIL {
            state.s_level[s] = 0;
            s_prcb_que.push(s);
        } else {
            state.s_level[s] = i32::MAX;
        }
    }
    state.s_level[graph.s_num_vtxs] = i32::MAX;

    let mut que_head = 0;
    while que_head < s_prcb_que.len() {
        let s = s_prcb_que[que_head];
        que_head += 1;
        if state.s_level[s] < state.s_level[graph.s_num_vtxs] {
            let s_begin = graph.s_idx[s];
            let s_end = graph.s_idx[s + 1];
            for k in s_begin..s_end {
                let t = graph.s_adj[k];
                let ss = if matching.t_mate[t as usize] == NIL {
                    graph.s_num_vtxs
                } else {
                    matching.t_mate[t as usize] as usize
                };
                if state.s_level[ss] == i32::MAX {
                    state.s_level[ss] = state.s_level[s] + 1;
                    if matching.t_mate[t as usize] != NIL {
                        s_prcb_que.push(matching.t_mate[t as usize] as usize);
                    }
                }
            }
        }
    }
    state.s_level[graph.s_num_vtxs] != i32::MAX
}

/*
 * DFS: iterative with edge index.
 *
 * state.s_idx[s] is an offset WITHIN s's adjacency range [graph.s_idx[s], graph.s_idx[s+1]).
 * So the "current candidate edge" is graph.s_adj[graph.s_idx[s] + state.s_idx[s]].
 */
fn dfs(s_first: usize, graph: &BipartiteGraph, matching: &mut BipartiteMatching, state: &mut HKIState) -> bool {
    state.stk_top = 0;
    state.s_prcb_stk[state.stk_top] = s_first as i32;
    state.stk_top += 1;

    while state.stk_top > 0 {
        let s = state.s_prcb_stk[state.stk_top - 1] as usize;
        let s_begin = graph.s_idx[s];
        let s_end = graph.s_idx[s + 1];
        let s_num_edgs = s_end - s_begin;

        let mut pushed = false;
        while state.s_idx[s] < s_num_edgs {
            let t = graph.s_adj[s_begin + state.s_idx[s]];
            let ss = if matching.t_mate[t as usize] == NIL {
                graph.s_num_vtxs
            } else {
                matching.t_mate[t as usize] as usize
            };
            if state.s_level[ss] != state.s_level[s] + 1 {
                state.s_idx[s] += 1;
                continue;
            }

            state.t_prcb_stk[state.stk_top - 1] = t;
            state.s_idx[s] += 1;

            if matching.t_mate[t as usize] == NIL {
                /* Found augmenting path — augment all the way back */
                let mut k = state.stk_top as isize - 1;
                while k >= 0 {
                    let ks = state.s_prcb_stk[k as usize] as usize;
                    let kt = state.t_prcb_stk[k as usize] as usize;
                    matching.t_mate[kt] = ks as i32;
                    matching.s_mate[ks] = kt as i32;
                    k -= 1;
                }
                return true;
            }

            state.s_prcb_stk[state.stk_top] = matching.t_mate[t as usize];
            state.stk_top += 1;
            pushed = true;
            break;
        }

        if !pushed {
            state.s_level[s] = i32::MAX;
            state.stk_top -= 1;
        }
    }
    false
}

/* ---------- Top-level Hopcroft-Karp Iterative ---------- */

fn hk_iterative(graph: &BipartiteGraph, matching: &mut BipartiteMatching) -> i32 {
    let n = graph.s_num_vtxs.max(1);
    let mut state = HKIState {
        s_level: vec![0i32; graph.s_num_vtxs + 1],
        s_idx: vec![0usize; graph.s_num_vtxs],
        s_prcb_stk: vec![0i32; n],
        t_prcb_stk: vec![0i32; n],
        stk_top: 0,
    };

    let mut num_phases: i32 = 0;
    let mut new_edgs: usize = 0;
    while bfs(graph, matching, &mut state) {
        num_phases += 1;
        for s in 0..graph.s_num_vtxs { state.s_idx[s] = 0; }
        for s in 0..graph.s_num_vtxs {
            if matching.s_mate[s] == NIL && dfs(s, graph, matching, &mut state) {
                new_edgs += 1;
            }
        }
    }
    matching.num_edgs += new_edgs;
    num_phases
}

/* ---------- Validation ---------- */

fn validate_bipartite_matching(graph: &BipartiteGraph, matching: &BipartiteMatching) {
    let mut errors = 0;
    let mut matched_s = 0;
    let mut matched_t = 0;

    for s in 0..graph.s_num_vtxs {
        if matching.s_mate[s] != NIL {
            matched_s += 1;
            let t = matching.s_mate[s];
            if t < 0 || (t as usize) >= graph.t_num_vtxs {
                eprintln!("ERROR: s_mate[{}] = {} out of range", s, t);
                errors += 1;
            } else if matching.t_mate[t as usize] != s as i32 {
                eprintln!("ERROR: s_mate[{}]={} but t_mate[{}]={}", s, t, t, matching.t_mate[t as usize]);
                errors += 1;
            } else {
                let s_begin = graph.s_idx[s];
                let s_end = graph.s_idx[s + 1];
                if !graph.s_adj[s_begin..s_end].binary_search(&t).is_ok() {
                    eprintln!("ERROR: edge ({},{}) not in graph", s, t);
                    errors += 1;
                }
            }
        }
    }
    for t in 0..graph.t_num_vtxs {
        if matching.t_mate[t] != NIL { matched_t += 1; }
    }

    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.num_edgs);
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
            let s: usize = parts[0].parse()?;
            let t: usize = parts[1].parse()?;
            edges.push((s, t));
        }
    }
    Ok((s_num_vtxs, t_num_vtxs, edges))
}

/* ---------- Main ---------- */

fn main() {
    println!("Hopcroft-Karp Iterative Algorithm - Rust Implementation (CSR)");
    println!("===============================================================\n");

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

            let graph = build_bipartite_graph(s_num_vtxs, t_num_vtxs, &edges);
            let mut matching = empty_bipartite_matching(&graph);

            let start = Instant::now();

            let greedy_size: usize = match greedy_mode {
                1 => greedy_init(&graph, &mut matching),
                2 => greedy_init_md(&graph, &mut matching),
                _ => 0,
            };

            let num_phases = hk_iterative(&graph, &mut matching);

            let duration = start.elapsed();

            validate_bipartite_matching(&graph, &matching);

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

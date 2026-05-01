/*
 * Hopcroft-Karp Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * VV (vector-of-vectors) adjacency.
 *
 * Refactored to separate input (BipartiteGraph), output (BipartiteMatching),
 * and algorithm state (HKState). The decoupling lets us drop the .clone()
 * calls on the inner adjacency that the original (single-struct) Rust port
 * needed to satisfy the borrow checker.
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
    s_adj: Vec<Vec<i32>>,    // s_adj[u] = t-neighbors of s-vertex u
    t_adj: Vec<Vec<i32>>,    // t_adj[v] = s-neighbors of t-vertex v
}

/* ---------- Output: BipartiteMatching ---------- */

#[allow(dead_code)]
struct BipartiteMatching {
    s_num_vtxs: usize,
    t_num_vtxs: usize,
    num_edgs: usize,
    s_mate: Vec<i32>,        // s_mate[u] = t-partner of s-vertex u, or NIL
    t_mate: Vec<i32>,        // t_mate[v] = s-partner of t-vertex v, or NIL
}

/* ---------- State: HKState ---------- */

struct HKState {
    dist: Vec<i32>,          // length s_num_vtxs+1; dist[s_num_vtxs] is the NIL sentinel
}

/* ---------- BipartiteGraph construction ---------- */

fn build_bipartite_graph(s_num_vtxs: usize, t_num_vtxs: usize,
                         edges: &[(usize, usize)]) -> BipartiteGraph {
    let mut s_adj: Vec<Vec<i32>> = vec![Vec::new(); s_num_vtxs];
    let mut t_adj: Vec<Vec<i32>> = vec![Vec::new(); t_num_vtxs];
    for &(u, v) in edges {
        if u < s_num_vtxs && v < t_num_vtxs {
            s_adj[u].push(v as i32);
            t_adj[v].push(u as i32);
        }
    }
    for adj in &mut s_adj { adj.sort_unstable(); adj.dedup(); }
    for adj in &mut t_adj { adj.sort_unstable(); adj.dedup(); }

    let num_edgs: usize = s_adj.iter().map(|a| a.len()).sum();

    BipartiteGraph {
        s_num_vtxs,
        t_num_vtxs,
        num_edgs,
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

/* ---------- Greedy initial matching: simple ---------- */

fn greedy_init(g: &BipartiteGraph, m: &mut BipartiteMatching) -> usize {
    let mut cnt: usize = 0;
    for u in 0..g.s_num_vtxs {
        if m.s_mate[u] != NIL { continue; }
        for &v in &g.s_adj[u] {
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
        for &v in &g.s_adj[u] {
            deg[v as usize] += 1;
        }
    }

    let mut order: Vec<usize> = (0..g.s_num_vtxs).collect();
    order.sort_unstable_by(|&a, &b| {
        g.s_adj[a].len().cmp(&g.s_adj[b].len()).then(a.cmp(&b))
    });

    for u in order {
        if m.s_mate[u] != NIL { continue; }
        let mut best: i32 = NIL;
        let mut best_deg = usize::MAX;
        for &v in &g.s_adj[u] {
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

/* ---------- HK BFS ---------- */

fn bfs(g: &BipartiteGraph, m: &BipartiteMatching, s: &mut HKState) -> bool {
    let mut queue: Vec<usize> = Vec::with_capacity(g.s_num_vtxs);

    for u in 0..g.s_num_vtxs {
        if m.s_mate[u] == NIL {
            s.dist[u] = 0;
            queue.push(u);
        } else {
            s.dist[u] = i32::MAX;
        }
    }
    s.dist[g.s_num_vtxs] = i32::MAX;

    let mut qi = 0;
    while qi < queue.len() {
        let u = queue[qi];
        qi += 1;
        if s.dist[u] < s.dist[g.s_num_vtxs] {
            for &v in &g.s_adj[u] {
                let paired = if m.t_mate[v as usize] == NIL {
                    g.s_num_vtxs
                } else {
                    m.t_mate[v as usize] as usize
                };
                if s.dist[paired] == i32::MAX {
                    s.dist[paired] = s.dist[u] + 1;
                    if m.t_mate[v as usize] != NIL {
                        queue.push(m.t_mate[v as usize] as usize);
                    }
                }
            }
        }
    }
    s.dist[g.s_num_vtxs] != i32::MAX
}

/* ---------- HK DFS ---------- */

fn dfs(u_opt: i32, g: &BipartiteGraph, m: &mut BipartiteMatching, s: &mut HKState) -> bool {
    if u_opt == NIL { return true; }
    let u = u_opt as usize;

    /* No clone: g.s_adj is borrowed via &g, m.t_mate/s_mate via &mut m -- disjoint borrows. */
    for i in 0..g.s_adj[u].len() {
        let v = g.s_adj[u][i];
        let paired = if m.t_mate[v as usize] == NIL {
            g.s_num_vtxs as i32
        } else {
            m.t_mate[v as usize]
        };
        if s.dist[paired as usize] == s.dist[u] + 1 {
            if dfs(m.t_mate[v as usize], g, m, s) {
                m.t_mate[v as usize] = u as i32;
                m.s_mate[u] = v;
                return true;
            }
        }
    }
    s.dist[u] = i32::MAX;
    false
}

/* ---------- Top-level Hopcroft-Karp ---------- */

fn hopcroft_karp(g: &BipartiteGraph, m: &mut BipartiteMatching) {
    let mut s = HKState {
        dist: vec![0i32; g.s_num_vtxs + 1],
    };

    let mut phases: i32 = 0;
    let mut new_edgs: usize = 0;
    while bfs(g, m, &mut s) {
        phases += 1;
        for u in 0..g.s_num_vtxs {
            if m.s_mate[u] == NIL && dfs(u as i32, g, m, &mut s) {
                new_edgs += 1;
            }
        }
    }
    m.num_edgs += new_edgs;
    println!("Phases: {}", phases);
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
            } else if !g.s_adj[u].binary_search(&v).is_ok() {
                eprintln!("ERROR: edge ({},{}) not in graph", u, v);
                errors += 1;
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
    println!("Hopcroft-Karp Algorithm - Rust Implementation (VV)");
    println!("====================================================\n");

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

            let bipartite_graph = build_bipartite_graph(s_num_vtxs, t_num_vtxs, &edges);
            let mut bipartite_matching = empty_bipartite_matching(&bipartite_graph);

            let start = Instant::now();

            let greedy_size: usize = match greedy_mode {
                1 => greedy_init(&bipartite_graph, &mut bipartite_matching),
                2 => greedy_init_md(&bipartite_graph, &mut bipartite_matching),
                _ => 0,
            };

            hopcroft_karp(&bipartite_graph, &mut bipartite_matching);

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

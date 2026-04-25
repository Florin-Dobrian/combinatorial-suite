/*
 * Hopcroft-Karp Hybrid Algorithm - O(E√V) Maximum Bipartite Matching
 *
 * CSR variant: adjacency stored as flat (adj_off, adj_edges) arrays
 * instead of Vec<Vec<usize>>. Algorithm unchanged.
 *
 * Old HK's lean BFS (single dist[] array, sentinel trick, no status enums)
 * + iterative stack-based DFS with edge index array (no recursion, no rescan).
 *
 * Rust port. All integers, no hash containers, fully deterministic.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;

struct HopcroftKarpHybrid {
    left_count: usize,
    right_count: usize,
    greedy_size: usize,
    /* CSR adjacency: adj_edges[adj_off[u] .. adj_off[u+1]) */
    adj_off: Vec<usize>,
    adj_edges: Vec<usize>,
    pair_left: Vec<i32>,
    pair_right: Vec<i32>,
    dist: Vec<i32>,

    /* Edge index: per-vertex offset within adjacency range (persistent within a phase) */
    edge_idx: Vec<usize>,

    /* DFS stack */
    stk_u: Vec<i32>,
    stk_v: Vec<i32>,
    stk_top: usize,
}

impl HopcroftKarpHybrid {
    fn new(left_count: usize, right_count: usize, edges: &[(usize, usize)]) -> Self {
        let mut tmp: Vec<Vec<usize>> = vec![Vec::new(); left_count];
        for &(u, v) in edges {
            if u < left_count && v < right_count { tmp[u].push(v); }
        }
        for adj in &mut tmp { adj.sort_unstable(); adj.dedup(); }

        let mut adj_off = vec![0usize; left_count + 1];
        for i in 0..left_count { adj_off[i + 1] = adj_off[i] + tmp[i].len(); }
        let mut adj_edges = Vec::with_capacity(adj_off[left_count]);
        for i in 0..left_count { adj_edges.extend_from_slice(&tmp[i]); }

        HopcroftKarpHybrid {
            left_count,
            right_count,
            greedy_size: 0,
            adj_off,
            adj_edges,
            pair_left: vec![NIL; left_count],
            pair_right: vec![NIL; right_count],
            dist: vec![0; left_count + 1],
            edge_idx: vec![0usize; left_count],
            stk_u: vec![0; left_count.max(1)],
            stk_v: vec![0; left_count.max(1)],
            stk_top: 0,
        }
    }

    fn bfs(&mut self) -> bool {
        let mut queue: Vec<usize> = Vec::with_capacity(self.left_count);

        for u in 0..self.left_count {
            if self.pair_left[u] == NIL { self.dist[u] = 0; queue.push(u); }
            else { self.dist[u] = i32::MAX; }
        }
        self.dist[self.left_count] = i32::MAX;

        let mut qi = 0usize;
        while qi < queue.len() {
            let u = queue[qi]; qi += 1;
            if self.dist[u] < self.dist[self.left_count] {
                let s = self.adj_off[u];
                let e = self.adj_off[u + 1];
                for j in s..e {
                    let v = self.adj_edges[j];
                    let pn = if self.pair_right[v] == NIL { self.left_count } else { self.pair_right[v] as usize };
                    if self.dist[pn] == i32::MAX {
                        self.dist[pn] = self.dist[u] + 1;
                        if self.pair_right[v] != NIL { queue.push(self.pair_right[v] as usize); }
                    }
                }
            }
        }
        self.dist[self.left_count] != i32::MAX
    }

    /*
     * DFS: iterative with edge index.
     * edge_idx[u] is an offset within u's adjacency range [adj_off[u], adj_off[u+1]).
     */
    fn dfs(&mut self, root: usize) -> bool {
        self.stk_top = 0;
        self.stk_u[0] = root as i32;
        self.stk_top = 1;

        while self.stk_top > 0 {
            let u = self.stk_u[self.stk_top - 1] as usize;
            let s = self.adj_off[u];
            let e = self.adj_off[u + 1];
            let sz = e - s;

            let mut pushed = false;
            while self.edge_idx[u] < sz {
                let v = self.adj_edges[s + self.edge_idx[u]];
                let pn = if self.pair_right[v] == NIL { self.left_count } else { self.pair_right[v] as usize };
                if self.dist[pn] != self.dist[u] + 1 {
                    self.edge_idx[u] += 1;
                    continue;
                }

                self.stk_v[self.stk_top - 1] = v as i32;
                self.edge_idx[u] += 1;

                if self.pair_right[v] == NIL {
                    let mut d = self.stk_top as isize - 1;
                    while d >= 0 {
                        let du = self.stk_u[d as usize] as usize;
                        let dv = self.stk_v[d as usize] as usize;
                        self.pair_right[dv] = du as i32;
                        self.pair_left[du] = dv as i32;
                        d -= 1;
                    }
                    return true;
                }

                self.stk_u[self.stk_top] = self.pair_right[v];
                self.stk_top += 1;
                pushed = true;
                break;
            }

            if !pushed {
                self.dist[u] = i32::MAX;
                self.stk_top -= 1;
            }
        }
        false
    }

    fn greedy_init(&mut self) -> usize {
        let mut cnt = 0usize;
        for u in 0..self.left_count {
            if self.pair_left[u] != NIL { continue; }
            let s = self.adj_off[u];
            let e = self.adj_off[u + 1];
            for j in s..e {
                let v = self.adj_edges[j];
                if self.pair_right[v] == NIL {
                    self.pair_left[u] = v as i32;
                    self.pair_right[v] = u as i32;
                    cnt += 1;
                    break;
                }
            }
        }
        cnt
    }

    fn greedy_init_md(&mut self) -> usize {
        let mut cnt = 0usize;
        let mut deg = vec![0usize; self.right_count];
        for u in 0..self.left_count {
            let s = self.adj_off[u];
            let e = self.adj_off[u + 1];
            for j in s..e { deg[self.adj_edges[j]] += 1; }
        }
        let mut order: Vec<usize> = (0..self.left_count).collect();
        order.sort_unstable_by(|&a, &b| {
            let da = self.adj_off[a + 1] - self.adj_off[a];
            let db = self.adj_off[b + 1] - self.adj_off[b];
            da.cmp(&db).then(a.cmp(&b))
        });
        for u in order {
            if self.pair_left[u] != NIL { continue; }
            let mut best: i32 = -1;
            let mut best_deg = usize::MAX;
            let s = self.adj_off[u];
            let e = self.adj_off[u + 1];
            for j in s..e {
                let v = self.adj_edges[j];
                if self.pair_right[v] == NIL && deg[v] < best_deg {
                    best = v as i32;
                    best_deg = deg[v];
                }
            }
            if best >= 0 {
                self.pair_left[u] = best;
                self.pair_right[best as usize] = u as i32;
                cnt += 1;
            }
        }
        cnt
    }

    fn maximum_matching(&mut self, greedy_mode: i32) -> Vec<(usize, usize)> {
        let greedy_count = match greedy_mode {
            1 => self.greedy_init(),
            2 => self.greedy_init_md(),
            _ => 0,
        };
        self.greedy_size = greedy_count;

        let mut phases = 0;
        while self.bfs() {
            phases += 1;
            for u in 0..self.left_count { self.edge_idx[u] = 0; }
            for u in 0..self.left_count {
                if self.pair_left[u] == NIL { self.dfs(u); }
            }
        }

        println!("Phases: {}", phases);

        let mut matching: Vec<(usize, usize)> = Vec::new();
        for u in 0..self.left_count {
            if self.pair_left[u] != NIL { matching.push((u, self.pair_left[u] as usize)); }
        }
        matching.sort_unstable();
        matching
    }
}

fn validate_matching(
    left_count: usize, right_count: usize,
    adj_off: &[usize], adj_edges: &[usize],
    matching: &[(usize, usize)],
) {
    let mut ldeg = vec![0i32; left_count];
    let mut rdeg = vec![0i32; right_count];
    let mut errors = 0;

    for &(u, v) in matching {
        let s = adj_off[u];
        let e = adj_off[u + 1];
        if !adj_edges[s..e].binary_search(&v).is_ok() {
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
    let matched_l = ldeg.iter().filter(|&&d| d > 0).count();
    let matched_r = rdeg.iter().filter(|&&d| d > 0).count();

    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.len());
    println!("Left matched: {}, Right matched: {}", matched_l, matched_r);
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
    println!("Hopcroft-Karp Hybrid CSR Algorithm - Rust Implementation");
    println!("===========================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename> [--greedy|--greedy-md]", args[0]);
        std::process::exit(1);
    }

    let greedy_mode: i32 = if args.iter().any(|a| a == "--greedy-md") { 2 }
                          else if args.iter().any(|a| a == "--greedy") { 1 }
                          else { 0 };

    match load_graph(&args[1]) {
        Ok((left_count, right_count, edges)) => {
            println!("Graph: {} left, {} right, {} edges", left_count, right_count, edges.len());

            let start = Instant::now();
            let mut hk = HopcroftKarpHybrid::new(left_count, right_count, &edges);
            let matching = hk.maximum_matching(greedy_mode);
            let duration = start.elapsed();

            validate_matching(left_count, right_count, &hk.adj_off, &hk.adj_edges, &matching);

            println!("Matching size: {}", matching.len());
            if greedy_mode > 0 {
                let gs = hk.greedy_size;
                let fs = matching.len();
                println!("Greedy init size: {}", gs);
                if fs > 0 { println!("Greedy/Final: {:.2}%", 100.0 * gs as f64 / fs as f64); }
                else { println!("Greedy/Final: NA"); }
            }
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

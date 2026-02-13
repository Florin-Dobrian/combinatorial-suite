/*
 * Hopcroft-Karp Algorithm - O(EâˆšV) Maximum Bipartite Matching
 *
 * Rust implementation â€” fully deterministic, no hash containers.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;

struct HopcroftKarp {
    left_count: usize,
    greedy_size: usize,
    right_count: usize,
    graph: Vec<Vec<usize>>,
    pair_left: Vec<i32>,
    pair_right: Vec<i32>,
    dist: Vec<i32>,
}

impl HopcroftKarp {
    fn new(left_count: usize, right_count: usize, edges: &[(usize, usize)]) -> Self {
        let mut graph = vec![Vec::new(); left_count];
        for &(u, v) in edges {
            if u < left_count && v < right_count {
                graph[u].push(v);
            }
        }
        for adj in &mut graph { adj.sort_unstable(); adj.dedup(); }

        HopcroftKarp {
            left_count,
            greedy_size: 0,
            right_count,
            graph,
            pair_left: vec![NIL; left_count],
            pair_right: vec![NIL; right_count],
            dist: vec![0; left_count + 1],
        }
    }

    fn bfs(&mut self) -> bool {
        let mut queue = Vec::new();
        let mut qi = 0;

        for u in 0..self.left_count {
            if self.pair_left[u] == NIL {
                self.dist[u] = 0;
                queue.push(u);
            } else {
                self.dist[u] = i32::MAX;
            }
        }
        self.dist[self.left_count] = i32::MAX;

        while qi < queue.len() {
            let u = queue[qi];
            qi += 1;
            if self.dist[u] < self.dist[self.left_count] {
                for &v in &self.graph[u] {
                    let paired = if self.pair_right[v] == NIL {
                        self.left_count
                    } else {
                        self.pair_right[v] as usize
                    };
                    if self.dist[paired] == i32::MAX {
                        self.dist[paired] = self.dist[u] + 1;
                        if self.pair_right[v] != NIL {
                            queue.push(self.pair_right[v] as usize);
                        }
                    }
                }
            }
        }
        self.dist[self.left_count] != i32::MAX
    }

    fn dfs(&mut self, u_opt: i32) -> bool {
        if u_opt == NIL { return true; }
        let u = u_opt as usize;

        let neighbors = self.graph[u].clone();
        for &v in &neighbors {
            let paired = if self.pair_right[v] == NIL {
                self.left_count
            } else {
                self.pair_right[v] as usize
            };
            if self.dist[paired] == self.dist[u] + 1 {
                if self.dfs(self.pair_right[v]) {
                    self.pair_right[v] = u as i32;
                    self.pair_left[u] = v as i32;
                    return true;
                }
            }
        }
        self.dist[u] = i32::MAX;
        false
    }


    fn greedy_init(&mut self) -> usize {
        let mut cnt: usize = 0;
        for u in 0..self.left_count {
            if self.pair_left[u] != NIL { continue; }
            let neighbors: Vec<usize> = self.graph[u].clone();
            for &v in &neighbors {
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

    /* Min-degree greedy: match each exposed left vertex with lowest-degree unmatched right neighbor */
    fn greedy_init_md(&mut self) -> usize {
        let mut cnt: usize = 0;
        let mut deg = vec![0usize; self.right_count];
        for u in 0..self.left_count {
            for &v in &self.graph[u] {
                deg[v] += 1;
            }
        }
        let mut order: Vec<usize> = (0..self.left_count).collect();
        order.sort_unstable_by(|&a, &b| self.graph[a].len().cmp(&self.graph[b].len()).then(a.cmp(&b)));
        for u in order {
            if self.pair_left[u] != NIL { continue; }
            let mut best: i32 = -1;
            let mut best_deg = usize::MAX;
            let neighbors: Vec<usize> = self.graph[u].clone();
            for &v in &neighbors {
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
        self.greedy_size = match greedy_mode {
            1 => self.greedy_init(),
            2 => self.greedy_init_md(),
            _ => 0,
        };
        while self.bfs() {
            for u in 0..self.left_count {
                if self.pair_left[u] == NIL {
                    self.dfs(u as i32);
                }
            }
        }

        let mut matching = Vec::new();
        for u in 0..self.left_count {
            if self.pair_left[u] != NIL {
                matching.push((u, self.pair_left[u] as usize));
            }
        }
        matching.sort_unstable();
        matching
    }
}

fn validate_matching(
    left_count: usize, right_count: usize,
    graph: &[Vec<usize>], matching: &[(usize, usize)],
) {
    let mut left_deg = vec![0i32; left_count];
    let mut right_deg = vec![0i32; right_count];
    let mut errors = 0;

    for &(u, v) in matching {
        if !graph[u].binary_search(&v).is_ok() {
            eprintln!("ERROR: Edge ({}, {}) not in graph!", u, v);
            errors += 1;
        }
        left_deg[u] += 1;
        right_deg[v] += 1;
    }
    for i in 0..left_count {
        if left_deg[i] > 1 { eprintln!("ERROR: Left vertex {} in {} edges!", i, left_deg[i]); errors += 1; }
    }
    for i in 0..right_count {
        if right_deg[i] > 1 { eprintln!("ERROR: Right vertex {} in {} edges!", i, right_deg[i]); errors += 1; }
    }
    let matched_left = left_deg.iter().filter(|&&d| d > 0).count();
    let matched_right = right_deg.iter().filter(|&&d| d > 0).count();

    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.len());
    println!("Matched vertices: {} left, {} right", matched_left, matched_right);
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
    println!("Hopcroft-Karp Algorithm - Rust Implementation");
    println!("================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename> [--greedy|--greedy-md]", args[0]);
        std::process::exit(1);
    }

    let greedy_mode: i32 = if args.iter().any(|a| a == "--greedy-md") { 2 } else if args.iter().any(|a| a == "--greedy") { 1 } else { 0 };
    match load_graph(&args[1]) {
        Ok((left_count, right_count, edges)) => {
            println!("Graph: {} left, {} right, {} edges", left_count, right_count, edges.len());

            let start = Instant::now();
            let mut hk = HopcroftKarp::new(left_count, right_count, &edges);
            let matching = hk.maximum_matching(greedy_mode);
            let duration = start.elapsed();

            validate_matching(left_count, right_count, &hk.graph, &matching);

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

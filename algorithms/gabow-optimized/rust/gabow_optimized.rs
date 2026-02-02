/*
 * Gabow's Scaling Algorithm (Optimized) - O(E√V) Maximum Matching
 *
 * Rust implementation — fully deterministic, no hash containers.
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
    base: Vec<usize>,
    parent: Vec<i32>,
    source_bridge: Vec<i32>,
    target_bridge: Vec<i32>,
    edge_queue: Vec<Vec<(usize, usize)>>,
    delta: usize,
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
        for adj in &mut graph { adj.sort_unstable(); }

        GabowOptimized {
            n,
            graph,
            mate: vec![NIL; n],
            label: vec![UNLABELED; n],
            base: (0..n).collect(),
            parent: vec![NIL; n],
            source_bridge: vec![NIL; n],
            target_bridge: vec![NIL; n],
            edge_queue: vec![Vec::new(); n + 1],
            delta: 0,
        }
    }

    fn find_base(&mut self, v: usize) -> usize {
        if self.base[v] != v {
            self.base[v] = self.find_base(self.base[v]);
        }
        self.base[v]
    }

    fn find_lca(&mut self, u: usize, v: usize) -> Option<usize> {
        let mut marked = vec![false; self.n];

        let mut x = self.find_base(u);
        while self.mate[x] != NIL {
            marked[x] = true;
            if self.parent[self.mate[x] as usize] == NIL { break; }
            x = self.find_base(self.parent[self.mate[x] as usize] as usize);
        }
        marked[x] = true;

        let mut y = self.find_base(v);
        while self.mate[y] != NIL {
            if marked[y] { return Some(y); }
            if self.parent[self.mate[y] as usize] == NIL { break; }
            y = self.find_base(self.parent[self.mate[y] as usize] as usize);
        }

        if marked[y] { Some(y) } else { None }
    }

    fn shrink_path(&mut self, lca: usize, x: usize, y: usize) {
        let mut v = self.find_base(x);
        while v != lca {
            self.base[v] = lca;
            let mv = self.mate[v];
            if mv == NIL { break; }
            self.base[mv as usize] = lca;
            self.source_bridge[mv as usize] = x as i32;
            self.target_bridge[mv as usize] = y as i32;
            if self.parent[mv as usize] == NIL { break; }
            v = self.find_base(self.parent[mv as usize] as usize);
        }
    }

    fn scan_edge(&mut self, u: usize, v: usize) {
        if self.delta < self.edge_queue.len() {
            self.edge_queue[self.delta].push((u, v));
        }
    }

    fn phase_1(&mut self) -> bool {
        self.delta = 0;
        for q in &mut self.edge_queue { q.clear(); }

        for i in 0..self.n {
            self.base[i] = i;
            self.label[i] = if self.mate[i] == NIL { EVEN } else { UNLABELED };
            self.parent[i] = NIL;
            self.source_bridge[i] = NIL;
            self.target_bridge[i] = NIL;
        }

        for v in 0..self.n {
            if self.mate[v] == NIL {
                let neighbors = self.graph[v].clone();
                for u in neighbors { self.scan_edge(v, u); }
            }
        }

        while self.delta <= self.n {
            while let Some((mut x, mut y)) = self.edge_queue[self.delta].pop() {
                let mut bx = self.find_base(x);
                let mut by = self.find_base(y);

                if self.label[bx] != EVEN {
                    std::mem::swap(&mut x, &mut y);
                    std::mem::swap(&mut bx, &mut by);
                }

                if bx == by || self.label[bx] != EVEN { continue; }
                if y as i32 == self.mate[x] || self.label[by] == ODD { continue; }

                if self.label[by] == UNLABELED {
                    let z = self.mate[y];
                    if z != NIL {
                        let z = z as usize;
                        self.label[y] = ODD;
                        self.label[z] = EVEN;
                        self.parent[y] = x as i32;
                        self.parent[z] = y as i32;

                        let neighbors = self.graph[z].clone();
                        for w in neighbors { self.scan_edge(z, w); }
                    }
                } else if self.label[by] == EVEN {
                    if let Some(lca) = self.find_lca(x, y) {
                        self.shrink_path(lca, x, y);
                        self.shrink_path(lca, y, x);
                    } else {
                        return true;
                    }
                }
            }
            self.delta += 1;
        }
        false
    }

    fn phase_2(&mut self) {
        for start in 0..self.n {
            if self.mate[start] != NIL || self.label[start] != EVEN { continue; }

            let mut queue = vec![start];
            let mut pred = vec![NIL; self.n];
            let mut vis = vec![false; self.n];
            let mut qi = 0;

            vis[self.find_base(start)] = true;
            let mut endpoint: Option<usize> = None;

            while qi < queue.len() && endpoint.is_none() {
                let u = queue[qi];
                qi += 1;

                let neighbors = self.graph[u].clone();
                for v in neighbors {
                    let bu = self.find_base(u);
                    let bv = self.find_base(v);
                    if bu == bv || vis[bv] { continue; }

                    if self.mate[v] == NIL && v != start {
                        pred[v] = u as i32;
                        endpoint = Some(v);
                        break;
                    }

                    if self.label[bv] != ODD {
                        pred[v] = u as i32;
                        vis[bv] = true;
                        let mv = self.mate[v];
                        if mv != NIL && !vis[self.find_base(mv as usize)] {
                            pred[mv as usize] = v as i32;
                            vis[self.find_base(mv as usize)] = true;
                            queue.push(mv as usize);
                        }
                    }
                }
            }

            if let Some(mut curr) = endpoint {
                let mut path = Vec::new();
                loop {
                    path.push(curr);
                    if pred[curr] == NIL { break; }
                    curr = pred[curr] as usize;
                }
                path.reverse();

                let mut i = 0;
                while i + 1 < path.len() {
                    self.mate[path[i]] = path[i + 1] as i32;
                    self.mate[path[i + 1]] = path[i] as i32;
                    i += 2;
                }
            }
        }
    }

    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        while self.phase_1() {
            self.phase_2();
        }

        let mut matching = Vec::new();
        for u in 0..self.n {
            if self.mate[u] != NIL && (self.mate[u] as usize) > u {
                matching.push((u, self.mate[u] as usize));
            }
        }
        matching.sort_unstable();
        matching
    }
}

fn validate_matching(n: usize, graph: &[Vec<usize>], matching: &[(usize, usize)]) {
    let mut deg = vec![0i32; n];
    let mut errors = 0;

    for &(u, v) in matching {
        if !graph[u].binary_search(&v).is_ok() {
            eprintln!("ERROR: Edge ({}, {}) not in graph!", u, v);
            errors += 1;
        }
        deg[u] += 1;
        deg[v] += 1;
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

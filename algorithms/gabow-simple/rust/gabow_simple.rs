/*
 * Gabow's Algorithm (Simple) - O(VE) Maximum Matching
 *
 * Rust implementation — fully deterministic, no hash containers.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;

struct GabowSimple {
    n: usize,
    graph: Vec<Vec<usize>>,
    mate: Vec<i32>,
    base: Vec<usize>,
    parent: Vec<i32>,
    blossom: Vec<bool>,
    visited: Vec<bool>,
}

impl GabowSimple {
    fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut graph = vec![Vec::new(); n];
        for &(u, v) in edges {
            if u < n && v < n && u != v {
                graph[u].push(v);
                graph[v].push(u);
            }
        }
        for adj in &mut graph { adj.sort_unstable(); }

        GabowSimple {
            n,
            graph,
            mate: vec![NIL; n],
            base: vec![0; n],
            parent: vec![NIL; n],
            blossom: vec![false; n],
            visited: vec![false; n],
        }
    }

    fn find_base(&mut self, v: usize) -> usize {
        if self.base[v] != v {
            self.base[v] = self.find_base(self.base[v]);
        }
        self.base[v]
    }

    fn find_lca(&mut self, mut u: usize, mut v: usize) -> Option<usize> {
        let mut path = vec![false; self.n];

        for _ in 0..self.n {
            u = self.find_base(u);
            path[u] = true;
            if self.mate[u] == NIL { break; }
            let mu = self.mate[u] as usize;
            if self.parent[mu] == NIL { break; }
            u = self.parent[mu] as usize;
        }

        for _ in 0..self.n {
            v = self.find_base(v);
            if path[v] { return Some(v); }
            if self.mate[v] == NIL { break; }
            let mv = self.mate[v] as usize;
            if self.parent[mv] == NIL { break; }
            v = self.parent[mv] as usize;
        }

        None
    }

    fn mark_blossom(&mut self, mut u: usize, lca: usize, queue: &mut Vec<usize>, qi: &mut usize) {
        let _ = qi; /* queue grows, qi stays valid */
        for _ in 0..self.n {
            if self.find_base(u) == lca { break; }
            let bv = self.find_base(u);
            let mu = self.mate[u] as usize;
            let bw = self.find_base(mu);

            self.blossom[bv] = true;
            self.blossom[bw] = true;

            if !self.visited[bw] {
                self.visited[bw] = true;
                queue.push(bw);
            }

            if self.parent[mu] == NIL { break; }
            u = self.parent[mu] as usize;
        }
    }

    fn contract_blossom(&mut self, u: usize, v: usize, queue: &mut Vec<usize>, qi: &mut usize) {
        if let Some(lca) = self.find_lca(u, v) {
            self.blossom = vec![false; self.n];
            self.mark_blossom(u, lca, queue, qi);
            self.mark_blossom(v, lca, queue, qi);

            for i in 0..self.n {
                let bi = self.find_base(i);
                if self.blossom[bi] {
                    self.base[i] = lca;
                    if !self.visited[i] {
                        self.visited[i] = true;
                        queue.push(i);
                    }
                }
            }
        }
    }

    fn find_augmenting_path(&mut self, start: usize) -> bool {
        for i in 0..self.n { self.base[i] = i; self.parent[i] = NIL; }
        self.visited = vec![false; self.n];

        let mut queue = vec![start];
        self.visited[start] = true;
        let mut qi = 0;

        while qi < queue.len() {
            let u = queue[qi];
            qi += 1;

            let neighbors = self.graph[u].clone();
            for &v in &neighbors {
                let bu = self.find_base(u);
                let bv = self.find_base(v);
                if bu == bv { continue; }

                if self.mate[v] == NIL {
                    self.parent[v] = u as i32;
                    return true;
                }

                if !self.visited[bv] {
                    self.parent[v] = u as i32;
                    self.visited[bv] = true;
                    let w = self.mate[v] as usize;
                    let bw = self.find_base(w);
                    self.visited[bw] = true;
                    queue.push(w);
                } else {
                    /* check if same tree → blossom */
                    let mut ru = bu;
                    for _ in 0..self.n {
                        if self.mate[ru] == NIL { break; }
                        let mru = self.mate[ru] as usize;
                        if self.parent[mru] == NIL { break; }
                        ru = self.find_base(self.parent[mru] as usize);
                    }
                    let mut rv = bv;
                    for _ in 0..self.n {
                        if self.mate[rv] == NIL { break; }
                        let mrv = self.mate[rv] as usize;
                        if self.parent[mrv] == NIL { break; }
                        rv = self.find_base(self.parent[mrv] as usize);
                    }
                    if ru == rv {
                        self.contract_blossom(u, v, &mut queue, &mut qi);
                    }
                }
            }
        }
        false
    }

    fn augment_path(&mut self, mut v: usize) {
        while self.parent[v] != NIL {
            let pv = self.parent[v] as usize;
            let ppv = self.mate[pv];
            self.mate[v] = pv as i32;
            self.mate[pv] = v as i32;
            if ppv == NIL { break; }
            v = ppv as usize;
        }
    }

    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        let mut found = true;
        while found {
            found = false;
            for v in 0..self.n {
                if self.mate[v] == NIL {
                    if self.find_augmenting_path(v) {
                        for u in 0..self.n {
                            if self.mate[u] == NIL && self.parent[u] != NIL {
                                self.augment_path(u);
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
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
    println!("Gabow's Algorithm (Simple) - Rust Implementation");
    println!("==================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename>", args[0]);
        std::process::exit(1);
    }

    match load_graph(&args[1]) {
        Ok((n, edges)) => {
            println!("Graph: {} vertices, {} edges", n, edges.len());

            let start = Instant::now();
            let mut gabow = GabowSimple::new(n, &edges);
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

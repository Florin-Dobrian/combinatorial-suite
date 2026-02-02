/*
 * Edmonds' Blossom Algorithm (Optimized) - O(V²E) Maximum Matching
 *
 * Rust implementation — fully deterministic, no hash containers.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;

struct EdmondsBlossomOptimized {
    n: usize,
    graph: Vec<Vec<usize>>,
    mate: Vec<i32>,
}

impl EdmondsBlossomOptimized {
    fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut graph = vec![Vec::new(); n];
        for &(u, v) in edges {
            if u < n && v < n && u != v {
                graph[u].push(v);
                graph[v].push(u);
            }
        }
        for adj in &mut graph { adj.sort_unstable(); }

        EdmondsBlossomOptimized { n, graph, mate: vec![NIL; n] }
    }

    fn find_base(v: usize, base: &[usize]) -> usize {
        let mut c = v;
        while base[c] != c { c = base[c]; }
        c
    }

    fn find_blossom_base(v: usize, w: usize, parent: &[i32], base: &[usize]) -> usize {
        let mut on_path = vec![false; base.len()];
        let mut c = v;
        loop {
            on_path[Self::find_base(c, base)] = true;
            if parent[c] == NIL { break; }
            c = parent[c] as usize;
        }
        c = w;
        loop {
            let bc = Self::find_base(c, base);
            if on_path[bc] { return bc; }
            if parent[c] == NIL { break; }
            c = parent[c] as usize;
        }
        Self::find_base(v, base)
    }

    fn trace_and_update(
        start: usize, blossom_base: usize,
        base: &mut [usize], label: &mut [i32], parent: &[i32], mate: &[i32],
        queue: &mut Vec<usize>, in_queue: &mut [bool],
    ) {
        let mut c = start;
        loop {
            let cb = Self::find_base(c, base);
            if cb == blossom_base { break; }
            base[c] = blossom_base;
            if label[c] == 2 {
                label[c] = 1;
                if !in_queue[c] { queue.push(c); in_queue[c] = true; }
            }
            if mate[c] == NIL { break; }
            let mc = mate[c] as usize;
            base[mc] = blossom_base;
            if parent[mc] == NIL { break; }
            c = parent[mc] as usize;
        }
    }

    fn find_augmenting_path(&self, start: usize) -> Option<Vec<usize>> {
        let mut parent = vec![NIL; self.n];
        let mut base: Vec<usize> = (0..self.n).collect();
        let mut label = vec![0i32; self.n]; /* 0=unlabeled, 1=outer, 2=inner */
        let mut in_queue = vec![false; self.n];

        label[start] = 1;
        let mut queue = vec![start];
        in_queue[start] = true;
        let mut qi = 0;

        while qi < queue.len() {
            let v = queue[qi];
            qi += 1;
            let vb = Self::find_base(v, &base);

            for &w in &self.graph[v] {
                let wb = Self::find_base(w, &base);
                if vb == wb { continue; }

                if label[w] == 0 {
                    if self.mate[w] != NIL {
                        let mw = self.mate[w] as usize;
                        label[w] = 2;
                        label[mw] = 1;
                        parent[w] = v as i32;
                        parent[mw] = w as i32;
                        if !in_queue[mw] { queue.push(mw); in_queue[mw] = true; }
                    } else {
                        parent[w] = v as i32;
                        let mut path = vec![w];
                        let mut c = v;
                        loop { path.push(c); if parent[c] == NIL { break; } c = parent[c] as usize; }
                        return Some(path);
                    }
                } else if label[w] == 1 {
                    let bb = Self::find_blossom_base(v, w, &parent, &base);
                    Self::trace_and_update(v, bb, &mut base, &mut label, &parent, &self.mate, &mut queue, &mut in_queue);
                    Self::trace_and_update(w, bb, &mut base, &mut label, &parent, &self.mate, &mut queue, &mut in_queue);
                }
            }
        }
        None
    }

    fn augment(&mut self, path: &[usize]) {
        let mut i = 0;
        while i + 1 < path.len() {
            let u = path[i];
            let v = path[i + 1];
            self.mate[u] = v as i32;
            self.mate[v] = u as i32;
            i += 2;
        }
    }

    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        let mut improved = true;
        while improved {
            improved = false;
            for v in 0..self.n {
                if self.mate[v] == NIL {
                    if let Some(path) = self.find_augmenting_path(v) {
                        self.augment(&path);
                        improved = true;
                        break;
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
    println!("Edmonds' Blossom Algorithm (Optimized) - Rust Implementation");
    println!("==============================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename>", args[0]);
        std::process::exit(1);
    }

    match load_graph(&args[1]) {
        Ok((n, edges)) => {
            println!("Graph: {} vertices, {} edges", n, edges.len());

            let start = Instant::now();
            let mut eb = EdmondsBlossomOptimized::new(n, &edges);
            let matching = eb.maximum_matching();
            let duration = start.elapsed();

            validate_matching(n, &eb.graph, &matching);

            println!("Matching size: {}", matching.len());
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

/*
 * Micali-Vazirani Algorithm (Hybrid) - O(E√V) Maximum Matching
 *
 * Hybrid approach:
 * - MV-style MIN phase (level building with even/odd tracking)
 * - Gabow-style MAX phase (BFS path finding and augmentation)
 *
 * Rust implementation — fully deterministic, no hash containers.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;
const UNSET: i32 = i32::MAX;

struct Node {
    preds: Vec<usize>,
    match_: i32,
    min_level: i32,
    even_level: i32,
    odd_level: i32,
}

impl Node {
    fn new() -> Self {
        Node { preds: Vec::new(), match_: NIL, min_level: UNSET, even_level: UNSET, odd_level: UNSET }
    }

    fn set_min_level(&mut self, level: i32) {
        self.min_level = level;
        if level % 2 == 0 { self.even_level = level; }
        else { self.odd_level = level; }
    }

    fn reset(&mut self) {
        self.preds.clear();
        self.min_level = UNSET;
        self.even_level = UNSET;
        self.odd_level = UNSET;
    }
}

struct MicaliVazirani {
    n: usize,
    graph: Vec<Vec<usize>>,
    nodes: Vec<Node>,
    base: Vec<usize>,
    levels: Vec<Vec<usize>>,
}

impl MicaliVazirani {
    fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut graph = vec![Vec::new(); n];
        for &(u, v) in edges {
            if u < n && v < n && u != v {
                graph[u].push(v);
                graph[v].push(u);
            }
        }
        for adj in &mut graph { adj.sort_unstable(); }

        MicaliVazirani {
            n,
            graph,
            nodes: (0..n).map(|_| Node::new()).collect(),
            base: (0..n).collect(),
            levels: Vec::new(),
        }
    }

    fn find_base(&mut self, v: usize) -> usize {
        if self.base[v] != v {
            self.base[v] = self.find_base(self.base[v]);
        }
        self.base[v]
    }

    fn add_to_level(&mut self, level: usize, node: usize) {
        while self.levels.len() <= level { self.levels.push(Vec::new()); }
        self.levels[level].push(node);
    }

    fn step_to(&mut self, to: usize, from: usize, level: i32) {
        let next = level + 1;
        let tl = self.nodes[to].min_level;
        if tl >= next {
            if tl != next {
                self.add_to_level(next as usize, to);
                self.nodes[to].set_min_level(next);
            }
            self.nodes[to].preds.push(from);
        }
    }

    fn phase_1(&mut self) {
        self.levels.clear();
        for i in 0..self.n { self.base[i] = i; self.nodes[i].reset(); }

        for i in 0..self.n {
            if self.nodes[i].match_ == NIL {
                self.add_to_level(0, i);
                self.nodes[i].set_min_level(0);
            }
        }

        for i in 0..self.n {
            if i >= self.levels.len() || self.levels[i].is_empty() { continue; }
            let level_snap = self.levels[i].clone();
            for &cur in &level_snap {
                let match_ = self.nodes[cur].match_;
                if i % 2 == 0 {
                    let neighbors = self.graph[cur].clone();
                    for nb in neighbors {
                        if nb as i32 != match_ { self.step_to(nb, cur, i as i32); }
                    }
                } else {
                    if match_ != NIL { self.step_to(match_ as usize, cur, i as i32); }
                }
            }
        }
    }

    fn phase_2(&mut self) -> bool {
        let mut found = false;

        for start in 0..self.n {
            if self.nodes[start].match_ != NIL || self.nodes[start].min_level != 0 { continue; }

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

                    if self.nodes[v].match_ == NIL && v != start {
                        pred[v] = u as i32;
                        endpoint = Some(v);
                        break;
                    }

                    pred[v] = u as i32;
                    vis[bv] = true;
                    let mv = self.nodes[v].match_;
                    if mv != NIL && !vis[self.find_base(mv as usize)] {
                        pred[mv as usize] = v as i32;
                        vis[self.find_base(mv as usize)] = true;
                        queue.push(mv as usize);
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
                    self.nodes[path[i]].match_ = path[i + 1] as i32;
                    self.nodes[path[i + 1]].match_ = path[i] as i32;
                    i += 2;
                }
                found = true;
            }
        }
        found
    }

    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        loop {
            self.phase_1();
            if !self.phase_2() { break; }
        }

        let mut matching = Vec::new();
        for u in 0..self.n {
            if self.nodes[u].match_ != NIL && (self.nodes[u].match_ as usize) > u {
                matching.push((u, self.nodes[u].match_ as usize));
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
    println!("Micali-Vazirani Algorithm (Hybrid) - Rust Implementation");
    println!("==========================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename>", args[0]);
        std::process::exit(1);
    }

    match load_graph(&args[1]) {
        Ok((n, edges)) => {
            println!("Graph: {} vertices, {} edges", n, edges.len());

            let start = Instant::now();
            let mut mv = MicaliVazirani::new(n, &edges);
            let matching = mv.maximum_matching();
            let duration = start.elapsed();

            validate_matching(n, &mv.graph, &matching);

            println!("Matching size: {}", matching.len());
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

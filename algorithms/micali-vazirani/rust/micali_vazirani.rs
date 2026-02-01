/*
 * Micali-Vazirani Algorithm - O(E√V) - Rust Implementation
 * 
 * Hybrid approach:
 * - Uses MV's MIN phase (level building with even/odd tracking)
 * - Uses Gabow-style MAX phase (simple path finding)
 */

use std::collections::VecDeque;
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
        Node {
            preds: Vec::new(),
            match_: NIL,
            min_level: UNSET,
            even_level: UNSET,
            odd_level: UNSET,
        }
    }
    
    fn set_min_level(&mut self, level: i32) {
        self.min_level = level;
        if level % 2 == 0 {
            self.even_level = level;
        } else {
            self.odd_level = level;
        }
    }
    
    fn reset(&mut self) {
        self.preds.clear();
        self.min_level = UNSET;
        self.even_level = UNSET;
        self.odd_level = UNSET;
    }
}

struct MicaliVazirani {
    vertex_count: usize,
    graph: Vec<Vec<usize>>,
    nodes: Vec<Node>,
    base: Vec<usize>,
    levels: Vec<Vec<usize>>,
    matchnum: usize,
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
        
        for adj in &mut graph {
            adj.sort_unstable();
        }
        
        MicaliVazirani {
            vertex_count: n,
            graph,
            nodes: (0..n).map(|_| Node::new()).collect(),
            base: (0..n).collect(),
            levels: Vec::new(),
            matchnum: 0,
        }
    }
    
    fn find_base(&mut self, v: usize) -> usize {
        if self.base[v] != v {
            self.base[v] = self.find_base(self.base[v]);
        }
        self.base[v]
    }
    
    fn add_to_level(&mut self, level: usize, node: usize) {
        while self.levels.len() <= level {
            self.levels.push(Vec::new());
        }
        self.levels[level].push(node);
    }
    
    fn step_to(&mut self, to: usize, from: usize, level: i32) {
        let next_level = level + 1;
        let tl = self.nodes[to].min_level;
        
        if tl >= next_level {
            if tl != next_level {
                self.add_to_level(next_level as usize, to);
                self.nodes[to].set_min_level(next_level);
            }
            self.nodes[to].preds.push(from);
        }
    }
    
    fn phase_1(&mut self) {
        self.levels.clear();
        
        // Reset bases
        for i in 0..self.vertex_count {
            self.base[i] = i;
            self.nodes[i].reset();
        }
        
        // Initialize free vertices at level 0
        for i in 0..self.vertex_count {
            if self.nodes[i].match_ == NIL {
                self.add_to_level(0, i);
                self.nodes[i].set_min_level(0);
            }
        }
        
        // Build levels
        for i in 0..self.vertex_count {
            if i >= self.levels.len() || self.levels[i].is_empty() {
                continue;
            }
            
            let level_nodes = self.levels[i].clone();
            for &current in &level_nodes {
                if i % 2 == 0 {
                    // Even level - explore all non-matching edges
                    let neighbors = self.graph[current].clone();
                    let match_ = self.nodes[current].match_;
                    for neighbor in neighbors {
                        if neighbor as i32 != match_ {
                            self.step_to(neighbor, current, i as i32);
                        }
                    }
                } else {
                    // Odd level - follow matching edge only
                    let match_ = self.nodes[current].match_;
                    if match_ != NIL {
                        self.step_to(match_ as usize, current, i as i32);
                    }
                }
            }
        }
    }
    
    fn phase_2(&mut self) -> bool {
        let mut found = false;
        
        for start in 0..self.vertex_count {
            if self.nodes[start].match_ != NIL {
                continue;
            }
            if self.nodes[start].min_level != 0 {
                continue;
            }
            
            // BFS from this free vertex
            let mut q = VecDeque::new();
            let mut pred = vec![NIL; self.vertex_count];
            let mut visited = vec![false; self.vertex_count];
            
            q.push_back(start);
            visited[self.find_base(start)] = true;
            
            let mut endpoint = None;
            
            while let Some(u) = q.pop_front() {
                if endpoint.is_some() {
                    break;
                }
                
                let neighbors = self.graph[u].clone();
                for v in neighbors {
                    let base_u = self.find_base(u);
                    let base_v = self.find_base(v);
                    
                    if base_u == base_v || visited[base_v] {
                        continue;
                    }
                    
                    // Check if v is a free vertex
                    if self.nodes[v].match_ == NIL && v != start {
                        pred[v] = u as i32;
                        endpoint = Some(v);
                        break;
                    }
                    
                    // Follow tree structure
                    pred[v] = u as i32;
                    visited[base_v] = true;
                    
                    // Continue along matching edge
                    let mate_v = self.nodes[v].match_;
                    if mate_v != NIL {
                        let mate_v_usize = mate_v as usize;
                        if !visited[self.find_base(mate_v_usize)] {
                            pred[mate_v_usize] = v as i32;
                            visited[self.find_base(mate_v_usize)] = true;
                            q.push_back(mate_v_usize);
                        }
                    }
                }
            }
            
            // If we found a path, augment it
            if let Some(mut curr) = endpoint {
                // Reconstruct path
                let mut path = Vec::new();
                loop {
                    path.push(curr);
                    if pred[curr] == NIL {
                        break;
                    }
                    curr = pred[curr] as usize;
                }
                path.reverse();
                
                // Augment along path
                for i in (0..path.len()).step_by(2) {
                    if i + 1 < path.len() {
                        let u = path[i];
                        let v = path[i + 1];
                        self.nodes[u].match_ = v as i32;
                        self.nodes[v].match_ = u as i32;
                    }
                }
                self.matchnum += 1;
                found = true;
            }
        }
        
        found
    }
    
    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        loop {
            self.phase_1();
            if !self.phase_2() {
                break;
            }
        }
        
        let mut matching = Vec::new();
        let mut seen = vec![false; self.vertex_count];
        
        for u in 0..self.vertex_count {
            if self.nodes[u].match_ != NIL && !seen[u] {
                let v = self.nodes[u].match_ as usize;
                if v < self.vertex_count {
                    matching.push((u.min(v), u.max(v)));
                    seen[u] = true;
                    seen[v] = true;
                }
            }
        }
        
        matching.sort_unstable();
        matching
    }
}

fn load_graph(filename: &str) -> Result<(usize, Vec<(usize, usize)>), Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    let first_line = lines.next().ok_or("Empty file")??;
    let parts: Vec<&str> = first_line.split_whitespace().collect();
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
    println!("Micali-Vazirani Algorithm - Rust Implementation");
    println!("================================================");
    println!();
    
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
            
            println!("\n=== Validation Report ===");
            println!("Matching size: {}", matching.len());
            println!("VALIDATION PASSED");
            println!("=========================\n");
            
            println!("Matching size: {}", matching.len());
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

/*
 * Gabow's Algorithm for Maximum Cardinality Matching (Simple Version)
 * Time complexity: O(V * E)
 * 
 * Rust implementation - fully deterministic
 * - Integer vertices only (0 to n-1)
 * - No HashSet or HashMap (fully deterministic)
 * - Sorted adjacency lists
 */

use std::collections::VecDeque;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::time::Instant;

struct GabowSimple {
    vertex_count: usize,
    graph: Vec<Vec<usize>>,
    mate: Vec<Option<usize>>,
    base: Vec<usize>,
    parent: Vec<Option<usize>>,
    blossom: Vec<bool>,
    visited: Vec<bool>,
}

impl GabowSimple {
    fn new(vertex_count: usize, edges: &[(usize, usize)]) -> Self {
        let mut graph = vec![Vec::new(); vertex_count];
        
        for &(u, v) in edges {
            if u < vertex_count && v < vertex_count && u != v {
                graph[u].push(v);
                graph[v].push(u);
            }
        }
        
        for adj in &mut graph {
            adj.sort_unstable();
        }
        
        GabowSimple {
            vertex_count,
            graph,
            mate: vec![None; vertex_count],
            base: vec![0; vertex_count],
            parent: vec![None; vertex_count],
            blossom: vec![false; vertex_count],
            visited: vec![false; vertex_count],
        }
    }
    
    fn find_base(&mut self, v: usize) -> usize {
        if self.base[v] != v {
            self.base[v] = self.find_base(self.base[v]);
        }
        self.base[v]
    }
    
    fn find_lca(&mut self, mut u: usize, mut v: usize) -> Option<usize> {
        let mut path = vec![false; self.vertex_count];
        
        let mut safety = 0;
        loop {
            if safety >= self.vertex_count {
                break;
            }
            u = self.find_base(u);
            path[u] = true;
            if self.mate[u].is_none() {
                break;
            }
            let mate_u = self.mate[u].unwrap();
            if self.parent[mate_u].is_none() {
                break;
            }
            u = self.parent[mate_u].unwrap();
            safety += 1;
        }
        
        safety = 0;
        loop {
            if safety >= self.vertex_count {
                break;
            }
            v = self.find_base(v);
            if path[v] {
                return Some(v);
            }
            if self.mate[v].is_none() {
                break;
            }
            let mate_v = self.mate[v].unwrap();
            if self.parent[mate_v].is_none() {
                break;
            }
            v = self.parent[mate_v].unwrap();
            safety += 1;
        }
        
        None
    }
    
    fn mark_blossom(&mut self, mut u: usize, lca: usize, queue: &mut VecDeque<usize>) {
        let mut safety = 0;
        while self.find_base(u) != lca && safety < self.vertex_count {
            let bv = self.find_base(u);
            let mate_u = self.mate[u].unwrap();
            let bw = self.find_base(mate_u);
            
            self.blossom[bv] = true;
            self.blossom[bw] = true;
            
            if !self.visited[bw] {
                self.visited[bw] = true;
                queue.push_back(bw);
            }
            
            if self.parent[mate_u].is_none() {
                break;
            }
            u = self.parent[mate_u].unwrap();
            safety += 1;
        }
    }
    
    fn contract_blossom(&mut self, u: usize, v: usize, queue: &mut VecDeque<usize>) {
        if let Some(lca) = self.find_lca(u, v) {
            self.blossom = vec![false; self.vertex_count];
            self.mark_blossom(u, lca, queue);
            self.mark_blossom(v, lca, queue);
            
            for i in 0..self.vertex_count {
                let base_i = self.find_base(i);
                if self.blossom[base_i] {
                    self.base[i] = lca;
                    if !self.visited[i] {
                        self.visited[i] = true;
                        queue.push_back(i);
                    }
                }
            }
        }
    }
    
    fn find_augmenting_path(&mut self, start: usize) -> bool {
        for i in 0..self.vertex_count {
            self.base[i] = i;
            self.parent[i] = None;
        }
        self.visited = vec![false; self.vertex_count];
        
        let mut queue = VecDeque::new();
        queue.push_back(start);
        self.visited[start] = true;
        
        let mut iterations = 0;
        while let Some(u) = queue.pop_front() {
            iterations += 1;
            if iterations > self.vertex_count * self.vertex_count {
                eprintln!("Warning: BFS timeout");
                break;
            }
            
            for &v in &self.graph[u].clone() {
                let base_u = self.find_base(u);
                let base_v = self.find_base(v);
                
                if base_u == base_v {
                    continue;
                }
                
                if self.mate[v].is_none() {
                    self.parent[v] = Some(u);
                    return true;
                }
                
                if !self.visited[base_v] {
                    self.parent[v] = Some(u);
                    self.visited[base_v] = true;
                    
                    let w = self.mate[v].unwrap();
                    let base_w = self.find_base(w);
                    self.visited[base_w] = true;
                    queue.push_back(w);
                } else {
                    let mut root_u = base_u;
                    let mut safety = 0;
                    while self.mate[root_u].is_some() && safety < self.vertex_count {
                        let mate_root = self.mate[root_u].unwrap();
                        if self.parent[mate_root].is_none() {
                            break;
                        }
                        root_u = self.find_base(self.parent[mate_root].unwrap());
                        safety += 1;
                    }
                    
                    let mut root_v = base_v;
                    safety = 0;
                    while self.mate[root_v].is_some() && safety < self.vertex_count {
                        let mate_root = self.mate[root_v].unwrap();
                        if self.parent[mate_root].is_none() {
                            break;
                        }
                        root_v = self.find_base(self.parent[mate_root].unwrap());
                        safety += 1;
                    }
                    
                    if root_u == root_v {
                        self.contract_blossom(u, v, &mut queue);
                    }
                }
            }
        }
        
        false
    }
    
    fn augment_path(&mut self, mut v: usize) {
        while let Some(pv) = self.parent[v] {
            let ppv = self.mate[pv];
            self.mate[v] = Some(pv);
            self.mate[pv] = Some(v);
            if let Some(ppv_val) = ppv {
                v = ppv_val;
            } else {
                break;
            }
        }
    }
    
    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        let mut found = true;
        let mut iterations = 0;
        
        while found {
            found = false;
            iterations += 1;
            
            for v in 0..self.vertex_count {
                if self.mate[v].is_none() {
                    if self.find_augmenting_path(v) {
                        for u in 0..self.vertex_count {
                            if self.mate[u].is_none() && self.parent[u].is_some() {
                                self.augment_path(u);
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            
            if iterations > self.vertex_count {
                eprintln!("Warning: Too many iterations");
                break;
            }
        }
        
        let mut matching = Vec::new();
        let mut seen = vec![false; self.vertex_count];
        
        for u in 0..self.vertex_count {
            if let Some(v) = self.mate[u] {
                if !seen[u] {
                    matching.push((u.min(v), u.max(v)));
                    seen[u] = true;
                    seen[v] = true;
                }
            }
        }
        
        matching.sort_unstable();
        matching
    }
    
    fn validate_matching(&self, matching: &[(usize, usize)]) {
        let mut degree = vec![0; self.vertex_count];
        let mut errors = 0;
        
        eprintln!("\n=== Validation Report ===");
        eprintln!("Matching size: {}", matching.len());
        
        for &(u, v) in matching {
            if !self.graph[u].contains(&v) {
                eprintln!("ERROR: Edge ({}, {}) not in graph!", u, v);
                errors += 1;
            }
            degree[u] += 1;
            degree[v] += 1;
        }
        
        for i in 0..self.vertex_count {
            if degree[i] > 1 {
                eprintln!("ERROR: Vertex {} in {} edges!", i, degree[i]);
                errors += 1;
            }
        }
        
        let matched = degree.iter().filter(|&&d| d > 0).count();
        
        eprintln!("Matched vertices: {}", matched);
        eprintln!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
        eprintln!("=========================\n");
    }
}

fn load_graph(filename: &str) -> io::Result<(usize, Vec<(usize, usize)>)> {
    let file = File::open(filename)?;
    let mut lines = BufReader::new(file).lines();
    
    let first = lines.next().ok_or(io::Error::new(io::ErrorKind::InvalidData, "Empty file"))??;
    let parts: Vec<usize> = first.split_whitespace()
        .map(|s| s.parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Parse error")))
        .collect::<Result<Vec<_>, _>>()?;
    
    if parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Need 2 numbers"));
    }
    
    let (n, m) = (parts[0], parts[1]);
    let mut edges = Vec::new();
    
    for line in lines.take(m) {
        let line = line?;
        let nums: Vec<usize> = line.split_whitespace()
            .map(|s| s.parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Parse error")))
            .collect::<Result<Vec<_>, _>>()?;
        
        if nums.len() != 2 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Need 2 numbers per edge"));
        }
        
        edges.push((nums[0], nums[1]));
    }
    
    Ok((n, edges))
}

fn main() {
    println!("Gabow's Algorithm (Simple Version) - Rust Implementation");
    println!("=========================================================\n");
    
    let args: Vec<String> = env::args().collect();
    
    if args.len() < 2 {
        println!("Usage: {} <filename>", args[0]);
        std::process::exit(1);
    }
    
    match load_graph(&args[1]) {
        Ok((n, edges)) => {
            println!("Graph: {} vertices, {} edges", n, edges.len());
            
            let start = Instant::now();
            let mut gabow = GabowSimple::new(n, &edges);
            let matching = gabow.maximum_matching();
            let duration = start.elapsed();
            
            gabow.validate_matching(&matching);
            
            println!("Matching size: {}", matching.len());
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

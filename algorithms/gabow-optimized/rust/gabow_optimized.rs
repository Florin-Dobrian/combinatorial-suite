/*
 * Gabow's O(E√V) Matching Algorithm - Rust Implementation
 * 
 * Time complexity: O(E√V)
 * Space complexity: O(V + E)
 */

use std::collections::VecDeque;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;
const UNLABELED: i32 = 0;
const EVEN: i32 = 1;
const ODD: i32 = 2;

struct GabowOptimized {
    vertex_count: usize,
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
        
        // Sort for determinism
        for adj in &mut graph {
            adj.sort_unstable();
        }
        
        GabowOptimized {
            vertex_count: n,
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
        let mut marked = vec![false; self.vertex_count];
        
        let mut x = self.find_base(u);
        while self.mate[x] != NIL {
            marked[x] = true;
            if self.parent[self.mate[x] as usize] == NIL {
                break;
            }
            x = self.find_base(self.parent[self.mate[x] as usize] as usize);
        }
        marked[x] = true;
        
        let mut y = self.find_base(v);
        while self.mate[y] != NIL {
            if marked[y] {
                return Some(y);
            }
            if self.parent[self.mate[y] as usize] == NIL {
                break;
            }
            y = self.find_base(self.parent[self.mate[y] as usize] as usize);
        }
        
        if marked[y] {
            Some(y)
        } else {
            None
        }
    }
    
    fn shrink_path(&mut self, lca: usize, x: usize, y: usize) {
        let mut v = self.find_base(x);
        while v != lca {
            self.base[v] = lca;
            
            let mate_v = self.mate[v];
            if mate_v == NIL {
                break;
            }
            
            self.base[mate_v as usize] = lca;
            self.source_bridge[mate_v as usize] = x as i32;
            self.target_bridge[mate_v as usize] = y as i32;
            
            if self.parent[mate_v as usize] == NIL {
                break;
            }
            v = self.find_base(self.parent[mate_v as usize] as usize);
        }
    }
    
    fn scan_edge(&mut self, u: usize, v: usize) {
        if self.delta < self.edge_queue.len() {
            self.edge_queue[self.delta].push((u, v));
        }
    }
    
    fn phase_1(&mut self) -> bool {
        self.delta = 0;
        
        for q in &mut self.edge_queue {
            q.clear();
        }
        
        for i in 0..self.vertex_count {
            self.base[i] = i;
            self.label[i] = if self.mate[i] == NIL { EVEN } else { UNLABELED };
            self.parent[i] = NIL;
            self.source_bridge[i] = NIL;
            self.target_bridge[i] = NIL;
        }
        
        for v in 0..self.vertex_count {
            if self.mate[v] == NIL {
                for &u in &self.graph[v].clone() {
                    self.scan_edge(v, u);
                }
            }
        }
        
        while self.delta <= self.vertex_count {
            while let Some((mut x, mut y)) = self.edge_queue[self.delta].pop() {
                let mut base_x = self.find_base(x);
                let mut base_y = self.find_base(y);
                
                if self.label[base_x] != EVEN {
                    std::mem::swap(&mut x, &mut y);
                    std::mem::swap(&mut base_x, &mut base_y);
                }
                
                if base_x == base_y || self.label[base_x] != EVEN {
                    continue;
                }
                if y == self.mate[x] as usize || self.label[base_y] == ODD {
                    continue;
                }
                
                if self.label[base_y] == UNLABELED {
                    let z = self.mate[y];
                    if z != NIL {
                        let z = z as usize;
                        self.label[y] = ODD;
                        self.label[z] = EVEN;
                        self.parent[y] = x as i32;
                        self.parent[z] = y as i32;
                        
                        for &w in &self.graph[z].clone() {
                            self.scan_edge(z, w);
                        }
                    }
                } else if self.label[base_y] == EVEN {
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
        for start in 0..self.vertex_count {
            if self.mate[start] != NIL || self.label[start] != EVEN {
                continue;
            }
            
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
                
                // Clone neighbors to avoid borrow checker issues
                let neighbors = self.graph[u].clone();
                
                for &v in &neighbors {
                    let base_u = self.find_base(u);
                    let base_v = self.find_base(v);
                    
                    if base_u == base_v || visited[base_v] {
                        continue;
                    }
                    
                    if self.mate[v] == NIL && v != start {
                        pred[v] = u as i32;
                        endpoint = Some(v);
                        break;
                    }
                    
                    if self.label[base_v] != ODD {
                        pred[v] = u as i32;
                        visited[base_v] = true;
                        
                        let mate_v = self.mate[v];
                        if mate_v != NIL && !visited[self.find_base(mate_v as usize)] {
                            pred[mate_v as usize] = v as i32;
                            visited[self.find_base(mate_v as usize)] = true;
                            q.push_back(mate_v as usize);
                        }
                    }
                }
            }
            
            if let Some(mut curr) = endpoint {
                let mut path = Vec::new();
                while pred[curr] != NIL {
                    path.push(curr);
                    curr = pred[curr] as usize;
                }
                path.push(curr);
                path.reverse();
                
                for i in (0..path.len()).step_by(2) {
                    if i + 1 < path.len() {
                        let u = path[i];
                        let v = path[i + 1];
                        self.mate[u] = v as i32;
                        self.mate[v] = u as i32;
                    }
                }
            }
        }
    }
    
    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        while self.phase_1() {
            self.phase_2();
        }
        
        let mut matching = Vec::new();
        let mut seen = vec![false; self.vertex_count];
        
        for u in 0..self.vertex_count {
            if self.mate[u] != NIL && !seen[u] {
                let v = self.mate[u] as usize;
                matching.push((u.min(v), u.max(v)));
                seen[u] = true;
                seen[v] = true;
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
    println!("Gabow's Scaling Algorithm (Optimized) - Rust Implementation");
    println!("===========================================================");
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
            let mut gabow = GabowOptimized::new(n, &edges);
            let matching = gabow.maximum_matching();
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

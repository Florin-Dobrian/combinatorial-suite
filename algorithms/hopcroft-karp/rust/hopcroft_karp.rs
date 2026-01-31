/*
 * Hopcroft-Karp Algorithm for Maximum Bipartite Matching (Rust Implementation)
 * Time complexity: O(E * sqrt(V))
 * 
 * Uses integers for vertices and deterministic data structures.
 */

use std::collections::VecDeque;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::time::Instant;

struct HopcroftKarp {
    left_count: usize,
    right_count: usize,
    graph: Vec<Vec<usize>>,  // graph[left_node] = list of right nodes
    pair_left: Vec<Option<usize>>,   // pair_left[u] = matched right node
    pair_right: Vec<Option<usize>>,  // pair_right[v] = matched left node
    dist: Vec<i32>,
}

impl HopcroftKarp {
    fn new(left_count: usize, right_count: usize, edges: &[(usize, usize)]) -> Self {
        let mut graph = vec![Vec::new(); left_count];
        
        // Build adjacency list
        for &(u, v) in edges {
            if u < left_count && v < right_count {
                graph[u].push(v);
            }
        }
        
        // Sort adjacency lists for deterministic iteration
        for adj in &mut graph {
            adj.sort_unstable();
        }
        
        HopcroftKarp {
            left_count,
            right_count,
            graph,
            pair_left: vec![None; left_count],
            pair_right: vec![None; right_count],
            dist: vec![0; left_count + 1],  // +1 for NIL at index left_count
        }
    }
    
    fn bfs(&mut self) -> bool {
        let mut queue = VecDeque::new();
        
        // Initialize distances and queue with unmatched left nodes
        for u in 0..self.left_count {
            if self.pair_left[u].is_none() {
                self.dist[u] = 0;
                queue.push_back(u);
            } else {
                self.dist[u] = i32::MAX;
            }
        }
        
        self.dist[self.left_count] = i32::MAX;  // NIL node at index left_count
        
        // BFS
        while let Some(u) = queue.pop_front() {
            if self.dist[u] < self.dist[self.left_count] {
                for &v in &self.graph[u] {
                    let paired_node = self.pair_right[v].unwrap_or(self.left_count);
                    
                    if self.dist[paired_node] == i32::MAX {
                        self.dist[paired_node] = self.dist[u] + 1;
                        if self.pair_right[v].is_some() {
                            queue.push_back(self.pair_right[v].unwrap());
                        }
                    }
                }
            }
        }
        
        self.dist[self.left_count] != i32::MAX
    }
    
    fn dfs(&mut self, u: Option<usize>) -> bool {
        if u.is_none() {
            return true;
        }
        
        let u = u.unwrap();
        let neighbors = self.graph[u].clone();  // Clone to avoid borrow checker issues
        
        for &v in &neighbors {
            let paired_node = self.pair_right[v].unwrap_or(self.left_count);
            
            if self.dist[paired_node] == self.dist[u] + 1 {
                if self.dfs(self.pair_right[v]) {
                    self.pair_right[v] = Some(u);
                    self.pair_left[u] = Some(v);
                    return true;
                }
            }
        }
        
        self.dist[u] = i32::MAX;
        false
    }
    
    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        // While there exist augmenting paths
        while self.bfs() {
            for u in 0..self.left_count {
                if self.pair_left[u].is_none() {
                    self.dfs(Some(u));
                }
            }
        }
        
        // Build the matching vector
        let mut matching = Vec::new();
        for u in 0..self.left_count {
            if let Some(v) = self.pair_left[u] {
                matching.push((u, v));
            }
        }
        
        matching.sort_unstable();
        matching
    }
    
    fn validate_matching(&self, matching: &[(usize, usize)]) {
        let mut left_degree = vec![0; self.left_count];
        let mut right_degree = vec![0; self.right_count];
        let mut errors = 0;
        
        eprintln!("\n=== Validation Report ===");
        eprintln!("Matching size (claimed): {}", matching.len());
        
        for &(u, v) in matching {
            // Check if edge exists in original graph
            if !self.graph[u].contains(&v) {
                eprintln!("ERROR: Edge ({}, {}) in matching but NOT in original graph!", u, v);
                errors += 1;
            }
            
            left_degree[u] += 1;
            right_degree[v] += 1;
        }
        
        for i in 0..self.left_count {
            if left_degree[i] > 1 {
                eprintln!("ERROR: Left node {} appears in {} edges (should be at most 1)!", i, left_degree[i]);
                errors += 1;
            }
        }
        
        for i in 0..self.right_count {
            if right_degree[i] > 1 {
                eprintln!("ERROR: Right node {} appears in {} edges (should be at most 1)!", i, right_degree[i]);
                errors += 1;
            }
        }
        
        let unique_left = left_degree.iter().filter(|&&d| d > 0).count();
        let unique_right = right_degree.iter().filter(|&&d| d > 0).count();
        
        eprintln!("Number of edges in matching: {}", matching.len());
        eprintln!("Left nodes matched: {}", unique_left);
        eprintln!("Right nodes matched: {}", unique_right);
        
        if errors > 0 {
            eprintln!("VALIDATION FAILED: {} errors found", errors);
        } else {
            eprintln!("VALIDATION PASSED: Matching is valid");
        }
        eprintln!("=========================\n");
    }
}

fn load_graph_from_file(filename: &str) -> io::Result<(usize, usize, Vec<(usize, usize)>)> {
    let file = File::open(filename)?;
    let mut lines = BufReader::new(file).lines();
    
    let first_line = lines.next().ok_or(io::Error::new(
        io::ErrorKind::InvalidData,
        "Empty file"
    ))??;
    
    let parts: Vec<usize> = first_line
        .split_whitespace()
        .map(|s| s.parse().map_err(|_| io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid number format"
        )))
        .collect::<Result<Vec<_>, _>>()?;
    
    if parts.len() != 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "First line must have 3 numbers"
        ));
    }
    
    let (left_count, right_count, edge_count) = (parts[0], parts[1], parts[2]);
    
    let mut edges = Vec::new();
    for line in lines.take(edge_count) {
        let line = line?;
        let nums: Vec<usize> = line
            .split_whitespace()
            .map(|s| s.parse().map_err(|_| io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid edge format"
            )))
            .collect::<Result<Vec<_>, _>>()?;
        
        if nums.len() != 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Each edge line must have 2 numbers"
            ));
        }
        
        edges.push((nums[0], nums[1]));
    }
    
    Ok((left_count, right_count, edges))
}

fn main() {
    println!("Hopcroft-Karp Algorithm - Rust Implementation");
    println!("==============================================\n");
    
    let args: Vec<String> = env::args().collect();
    
    if args.len() > 1 {
        let filename = &args[1];
        println!("Loading graph from: {}", filename);
        
        match load_graph_from_file(filename) {
            Ok((left_count, right_count, edges)) => {
                println!("File: {}", filename);
                println!("Graph: {} left nodes, {} right nodes, {} edges", 
                        left_count, right_count, edges.len());
                
                let start = Instant::now();
                let mut hk = HopcroftKarp::new(left_count, right_count, &edges);
                let matching = hk.maximum_matching();
                let duration = start.elapsed();
                
                hk.validate_matching(&matching);
                
                println!("Matching size: {}", matching.len());
                println!("Execution time: {} ms", duration.as_millis());
            }
            Err(e) => {
                eprintln!("Error loading file: {}", e);
                std::process::exit(1);
            }
        }
    } else {
        println!("Usage: {} <filename>", args[0]);
        std::process::exit(1);
    }
}

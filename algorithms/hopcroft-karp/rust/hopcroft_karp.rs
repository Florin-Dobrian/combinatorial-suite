/*
 * Hopcroft-Karp Algorithm for Maximum Bipartite Matching (Rust Implementation)
 * Time complexity: O(E * sqrt(V))
 */

use std::collections::{HashMap, HashSet, VecDeque};
use std::time::Instant;

#[derive(Clone)]
struct HopcroftKarp {
    left: HashSet<String>,
    graph: HashMap<String, Vec<String>>,
    pair_left: HashMap<String, Option<String>>,
    pair_right: HashMap<String, Option<String>>,
    dist: HashMap<String, i32>,
}

impl HopcroftKarp {
    fn new(
        left_nodes: &[String],
        right_nodes: &[String],
        edges: &[(String, String)],
    ) -> Self {
        let left: HashSet<String> = left_nodes.iter().cloned().collect();
        let right: HashSet<String> = right_nodes.iter().cloned().collect();
        
        let mut graph: HashMap<String, Vec<String>> = HashMap::new();
        
        // Build adjacency list
        for (u, v) in edges {
            if left.contains(u) && right.contains(v) {
                graph.entry(u.clone()).or_insert_with(Vec::new).push(v.clone());
            }
        }
        
        // Initialize pairs
        let mut pair_left = HashMap::new();
        let mut pair_right = HashMap::new();
        
        for u in &left {
            pair_left.insert(u.clone(), None);
        }
        for v in &right {
            pair_right.insert(v.clone(), None);
        }
        
        HopcroftKarp {
            left,
            graph,
            pair_left,
            pair_right,
            dist: HashMap::new(),
        }
    }
    
    fn bfs(&mut self) -> bool {
        let mut queue = VecDeque::new();
        
        // Initialize distances and queue with unmatched left nodes
        for u in &self.left {
            if self.pair_left.get(u).unwrap().is_none() {
                self.dist.insert(u.clone(), 0);
                queue.push_back(u.clone());
            } else {
                self.dist.insert(u.clone(), i32::MAX);
            }
        }
        
        self.dist.insert(String::new(), i32::MAX);
        
        // BFS
        while let Some(u) = queue.pop_front() {
            let u_dist = *self.dist.get(&u).unwrap();
            
            if u_dist < *self.dist.get("").unwrap() {
                if let Some(neighbors) = self.graph.get(&u) {
                    for v in neighbors {
                        let paired_node = self.pair_right.get(v).unwrap().clone();
                        
                        let paired_key = paired_node.clone().unwrap_or_default();
                        let paired_dist = *self.dist.get(&paired_key).unwrap_or(&i32::MAX);
                        
                        if paired_dist == i32::MAX {
                            self.dist.insert(paired_key.clone(), u_dist + 1);
                            if paired_node.is_some() {
                                queue.push_back(paired_key);
                            }
                        }
                    }
                }
            }
        }
        
        *self.dist.get("").unwrap() != i32::MAX
    }
    
    fn dfs(&mut self, u: Option<String>) -> bool {
        if u.is_none() {
            return true;
        }
        
        let u_str = u.unwrap();
        let u_dist = *self.dist.get(&u_str).unwrap();
        
        if let Some(neighbors) = self.graph.get(&u_str).cloned() {
            for v in neighbors {
                let paired_node = self.pair_right.get(&v).unwrap().clone();
                
                let paired_key = paired_node.clone().unwrap_or_default();
                let paired_dist = *self.dist.get(&paired_key).unwrap_or(&i32::MAX);
                
                if paired_dist == u_dist + 1 {
                    if self.dfs(paired_node) {
                        self.pair_right.insert(v.clone(), Some(u_str.clone()));
                        self.pair_left.insert(u_str.clone(), Some(v));
                        return true;
                    }
                }
            }
        }
        
        self.dist.insert(u_str, i32::MAX);
        false
    }
    
    fn maximum_matching(&mut self) -> Vec<(String, String)> {
        // While there exist augmenting paths
        while self.bfs() {
            let left_nodes: Vec<String> = self.left.iter().cloned().collect();
            for u in left_nodes {
                if self.pair_left.get(&u).unwrap().is_none() {
                    self.dfs(Some(u));
                }
            }
        }
        
        // Build the matching vector
        let mut matching = Vec::new();
        for u in &self.left {
            if let Some(Some(v)) = self.pair_left.get(u) {
                matching.push((u.clone(), v.clone()));
            }
        }
        
        matching
    }
}

// Generate a large random bipartite graph for benchmarking
fn generate_large_graph(
    left_size: usize,
    right_size: usize,
    edges_per_left_node: usize,
) -> Vec<(String, String)> {
    let mut edges = Vec::new();
    
    for i in 0..left_size {
        for j in 0..edges_per_left_node {
            let right_idx = (i * edges_per_left_node + j) % right_size;
            edges.push((
                format!("L{}", i),
                format!("R{}", right_idx),
            ));
        }
    }
    
    edges
}

fn main() {
    println!("Rust Hopcroft-Karp Implementation");
    println!("==================================\n");
    
    // Small example
    let left: Vec<String> = vec!["A", "B", "C", "D"]
        .into_iter()
        .map(|s| s.to_string())
        .collect();
    let right: Vec<String> = vec!["1", "2", "3", "4"]
        .into_iter()
        .map(|s| s.to_string())
        .collect();
    let edges: Vec<(String, String)> = vec![
        ("A", "1"), ("A", "2"),
        ("B", "2"), ("B", "3"),
        ("C", "3"), ("C", "4"),
        ("D", "4"),
    ]
    .into_iter()
    .map(|(a, b)| (a.to_string(), b.to_string()))
    .collect();
    
    let mut hk = HopcroftKarp::new(&left, &right, &edges);
    let matching = hk.maximum_matching();
    
    println!("Small example:");
    println!("Matching size: {}", matching.len());
    print!("Matching: ");
    for (u, v) in &matching {
        print!("({},{}) ", u, v);
    }
    println!("\n");
    
    // Benchmark with larger graph
    println!("Benchmarking with larger graph...");
    let left_size = 1000;
    let right_size = 1000;
    let edges_per_node = 10;
    
    let large_left: Vec<String> = (0..left_size)
        .map(|i| format!("L{}", i))
        .collect();
    let large_right: Vec<String> = (0..right_size)
        .map(|i| format!("R{}", i))
        .collect();
    
    let large_edges = generate_large_graph(left_size, right_size, edges_per_node);
    
    println!(
        "Graph size: {} left nodes, {} right nodes, {} edges",
        left_size,
        right_size,
        large_edges.len()
    );
    
    let start = Instant::now();
    let mut large_hk = HopcroftKarp::new(&large_left, &large_right, &large_edges);
    let large_matching = large_hk.maximum_matching();
    let duration = start.elapsed();
    
    println!("Matching size: {}", large_matching.len());
    println!("Execution time: {} ms", duration.as_millis());
}

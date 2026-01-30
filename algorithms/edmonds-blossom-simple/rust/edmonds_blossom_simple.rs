/*
 * Edmonds' Blossom Algorithm (Simple Version) for Maximum Cardinality Matching
 * Time complexity: O(V⁴)
 * 
 * This is the straightforward implementation that finds one augmenting path per iteration.
 */

use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::time::Instant;

struct EdmondsBlossomSimple {
    vertices: HashSet<String>,
    graph: HashMap<String, HashSet<String>>,
    mate: HashMap<String, Option<String>>,
}

impl EdmondsBlossomSimple {
    fn new(vertex_list: &[String], edges: &[(String, String)]) -> Self {
        let vertices: HashSet<String> = vertex_list.iter().cloned().collect();
        let mut graph: HashMap<String, HashSet<String>> = HashMap::new();
        
        // Build adjacency list (undirected graph)
        for (u, v) in edges {
            if vertices.contains(u) && vertices.contains(v) && u != v {
                graph.entry(u.clone()).or_insert_with(HashSet::new).insert(v.clone());
                graph.entry(v.clone()).or_insert_with(HashSet::new).insert(u.clone());
            }
        }
        
        // Initialize matching
        let mut mate = HashMap::new();
        for v in &vertices {
            mate.insert(v.clone(), None);
        }
        
        EdmondsBlossomSimple {
            vertices,
            graph,
            mate,
        }
    }
    
    fn is_outer(&self, v: &str, parent: &HashMap<String, String>) -> bool {
        if !parent.contains_key(v) {
            return false;
        }
        
        let mut current = v.to_string();
        let mut distance = 0;
        let mut visited = HashSet::new();
        
        while let Some(p) = parent.get(&current) {
            if &current == p || visited.contains(&current) {
                break;
            }
            visited.insert(current.clone());
            current = p.clone();
            distance += 1;
        }
        
        distance % 2 == 0
    }
    
    fn find_blossom(&self, v: &str, w: &str, 
                    parent: &HashMap<String, String>,
                    base: &HashMap<String, String>) -> String {
        // Mark all ancestors of v
        let mut path_v = HashSet::new();
        let mut current = v.to_string();
        
        while let Some(p) = parent.get(&current) {
            path_v.insert(base.get(&current).unwrap_or(&current).clone());
            if &current == p {
                break;
            }
            current = p.clone();
        }
        
        // Find first common ancestor from w
        current = w.to_string();
        while let Some(p) = parent.get(&current) {
            let base_current = base.get(&current).unwrap_or(&current);
            if path_v.contains(base_current) {
                return base_current.clone();
            }
            if &current == p {
                break;
            }
            current = p.clone();
        }
        
        v.to_string() // Fallback
    }
    
    fn contract_blossom(&self, blossom_base: &str, v: &str, w: &str,
                       base: &mut HashMap<String, String>,
                       parent: &HashMap<String, String>) {
        // Path from v to base
        let mut current = v.to_string();
        while base.get(&current).unwrap_or(&current) != blossom_base {
            base.insert(current.clone(), blossom_base.to_string());
            
            if let Some(Some(mate_current)) = self.mate.get(&current) {
                base.insert(mate_current.clone(), blossom_base.to_string());
                if let Some(next) = parent.get(mate_current) {
                    current = next.clone();
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        
        // Path from w to base
        current = w.to_string();
        while base.get(&current).unwrap_or(&current) != blossom_base {
            base.insert(current.clone(), blossom_base.to_string());
            
            if let Some(Some(mate_current)) = self.mate.get(&current) {
                base.insert(mate_current.clone(), blossom_base.to_string());
                if let Some(next) = parent.get(mate_current) {
                    current = next.clone();
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        
        base.insert(v.to_string(), blossom_base.to_string());
        base.insert(w.to_string(), blossom_base.to_string());
    }
    
    fn build_path(&self, v: &str, w: &str, parent: &HashMap<String, String>) -> Vec<String> {
        let mut path = vec![w.to_string(), v.to_string()];
        
        let mut current = v.to_string();
        while let Some(p) = parent.get(&current) {
            if &current == p {
                break;
            }
            path.push(p.clone());
            current = p.clone();
        }
        
        path
    }
    
    fn augment(&mut self, path: &[String]) {
        for i in (0..path.len().saturating_sub(1)).step_by(2) {
            let u = &path[i];
            let v = &path[i + 1];
            self.mate.insert(u.clone(), Some(v.clone()));
            self.mate.insert(v.clone(), Some(u.clone()));
        }
    }
    
    fn find_augmenting_path(&self, start: &str) -> Vec<String> {
        let mut parent: HashMap<String, String> = HashMap::new();
        let mut base: HashMap<String, String> = HashMap::new();
        
        parent.insert(start.to_string(), start.to_string());
        for v in &self.vertices {
            base.insert(v.clone(), v.clone());
        }
        
        let mut queue = VecDeque::new();
        queue.push_back(start.to_string());
        
        while let Some(v) = queue.pop_front() {
            if let Some(neighbors) = self.graph.get(&v) {
                for w in neighbors {
                    let base_v = base.get(&v).unwrap_or(&v);
                    let base_w = base.get(w).unwrap_or(w);
                    
                    // Skip if already in tree at same level
                    if base_v == base_w {
                        continue;
                    }
                    
                    // If w is unmatched, we found an augmenting path!
                    if self.mate.get(w).unwrap_or(&None).is_none() {
                        return self.build_path(&v, w, &parent);
                    }
                    
                    // If w is not in tree yet
                    if !parent.contains_key(w) {
                        parent.insert(w.clone(), v.clone());
                        if let Some(Some(mate_w)) = self.mate.get(w) {
                            parent.insert(mate_w.clone(), w.clone());
                            queue.push_back(mate_w.clone());
                        }
                    }
                    // If w is in tree at even distance (blossom found)
                    else if self.is_outer(w, &parent) {
                        let blossom_base = self.find_blossom(&v, w, &parent, &base);
                        self.contract_blossom(&blossom_base, &v, w, &mut base, &parent);
                        queue.push_back(blossom_base);
                    }
                }
            }
        }
        
        Vec::new()
    }
    
    fn maximum_matching(&mut self) -> Vec<(String, String)> {
        let mut improved = true;
        
        while improved {
            improved = false;
            
            let vertices: Vec<String> = self.vertices.iter().cloned().collect();
            for v in vertices {
                if self.mate.get(&v).unwrap_or(&None).is_none() {
                    let path = self.find_augmenting_path(&v);
                    if !path.is_empty() {
                        self.augment(&path);
                        improved = true;
                        break;
                    }
                }
            }
        }
        
        // Build the matching set
        let mut matching = Vec::new();
        let mut seen = HashSet::new();
        
        for u in &self.vertices {
            if let Some(Some(v)) = self.mate.get(u) {
                if !seen.contains(v) {
                    matching.push((u.clone(), v.clone()));
                    seen.insert(u.clone());
                }
            }
        }
        
        matching
    }
}

// Load graph from file
fn load_graph_from_file(filename: &str) -> io::Result<(Vec<String>, Vec<(String, String)>)> {
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
    
    if parts.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "First line must have 2 numbers"
        ));
    }
    
    let (vertex_count, edge_count) = (parts[0], parts[1]);
    
    let vertices: Vec<String> = (0..vertex_count).map(|i| format!("V{}", i)).collect();
    
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
        
        edges.push((format!("V{}", nums[0]), format!("V{}", nums[1])));
    }
    
    Ok((vertices, edges))
}

// Generate deterministic test graph
fn generate_test_graph(n: usize, edge_probability_percent: usize) -> (Vec<String>, Vec<(String, String)>) {
    let vertices: Vec<String> = (0..n).map(|i| format!("V{}", i)).collect();
    let mut edges = Vec::new();
    
    for i in 0..n {
        for j in (i + 1)..n {
            // Deterministic "random" based on i and j
            if ((i * 17 + j * 31) % 100) < edge_probability_percent {
                edges.push((vertices[i].clone(), vertices[j].clone()));
            }
        }
    }
    
    (vertices, edges)
}

fn run_example(vertices: &[String], edges: &[(String, String)], description: &str) {
    println!("{}", description);
    println!("Graph: {} vertices, {} edges", vertices.len(), edges.len());
    
    let start = Instant::now();
    let mut eb = EdmondsBlossomSimple::new(vertices, edges);
    let matching = eb.maximum_matching();
    let duration = start.elapsed();
    
    println!("Matching size: {}", matching.len());
    if matching.len() <= 10 {
        print!("Matching: ");
        for (u, v) in &matching {
            print!("({},{}) ", u, v);
        }
        println!();
    }
    println!("Execution time: {} ms", duration.as_millis());
    println!();
}

fn main() {
    println!("Edmonds' Blossom Algorithm (Simple) - Rust Implementation");
    println!("==========================================================\n");
    
    let args: Vec<String> = env::args().collect();
    
    // Check if a file was provided
    if args.len() > 1 {
        let filename = &args[1];
        println!("Loading graph from: {}", filename);
        
        match load_graph_from_file(filename) {
            Ok((vertices, edges)) => {
                run_example(&vertices, &edges, &format!("File: {}", filename));
            }
            Err(e) => {
                eprintln!("Error loading file: {}", e);
                std::process::exit(1);
            }
        }
    } else {
        // Run built-in examples
        println!("Running built-in examples (use: ./edmonds_blossom_simple_rust <filename> to load from file)\n");
        
        // Example 1: Simple graph with triangle
        let vertices: Vec<String> = vec!["A", "B", "C", "D", "E"]
            .into_iter()
            .map(|s| s.to_string())
            .collect();
        let edges: Vec<(String, String)> = vec![
            ("A", "B"), ("B", "C"), ("C", "A"),  // Triangle
            ("C", "D"), ("D", "E"),
        ]
        .into_iter()
        .map(|(a, b)| (a.to_string(), b.to_string()))
        .collect();
        
        run_example(&vertices, &edges, "Example 1: Graph with triangle (blossom)");
        
        // Example 2: Larger random graph
        let (vertices, edges) = generate_test_graph(20, 20);
        run_example(&vertices, &edges, "Example 2: Random graph (20 vertices)");
        
        // Benchmark
        println!("Benchmarking with larger graph...");
        let (vertices, edges) = generate_test_graph(100, 10);
        run_example(&vertices, &edges, "Benchmark: Random graph (100 vertices)");
    }
}

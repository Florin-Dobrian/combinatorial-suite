/*
 * Edmonds' Blossom Algorithm (Optimized Version) for Maximum Cardinality Matching
 * Time complexity: O(V²E)
 * 
 * This version uses efficient data structures and label tracking.
 * Uses indices instead of references to work with Rust's borrow checker.
 */

use std::collections::{HashSet, VecDeque};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::time::Instant;

struct EdmondsBlossomOptimized {
    vertex_count: usize,
    // Adjacency list: vertex index -> vector of neighbor indices
    graph: Vec<Vec<usize>>,
    // Matching: vertex index -> matched vertex index (or None)
    mate: Vec<Option<usize>>,
    // Vertex names for display
    vertex_names: Vec<String>,
}

impl EdmondsBlossomOptimized {
    fn new(vertices: &[String], edges: &[(usize, usize)]) -> Self {
        let vertex_count = vertices.len();
        let mut graph = vec![Vec::new(); vertex_count];
        
        // Build adjacency list
        for &(u, v) in edges {
            if u < vertex_count && v < vertex_count && u != v {
                graph[u].push(v);
                graph[v].push(u);
            }
        }
        
        // Sort adjacency lists for deterministic iteration
        for adj_list in &mut graph {
            adj_list.sort_unstable();
        }
        
        let mate = vec![None; vertex_count];
        let vertex_names = vertices.to_vec();
        
        EdmondsBlossomOptimized {
            vertex_count,
            graph,
            mate,
            vertex_names,
        }
    }
    
    fn find_base(v: usize, base: &[usize]) -> usize {
        let mut current = v;
        let mut visited = HashSet::new();
        
        // Follow base pointers until we reach a vertex that is its own base
        while base[current] != current && !visited.contains(&current) {
            visited.insert(current);
            current = base[current];
        }
        
        current
    }
    
    fn find_blossom_base(
        v: usize,
        w: usize,
        parent: &[Option<usize>],
        base: &[usize],
    ) -> usize {
        // Mark all ancestors of v
        let mut path_v = HashSet::new();
        let mut current = v;
        
        loop {
            let base_current = Self::find_base(current, base);
            path_v.insert(base_current);
            
            match parent[current] {
                Some(p) => current = p,
                None => break,
            }
        }
        
        // Find first common ancestor from w
        current = w;
        loop {
            let base_current = Self::find_base(current, base);
            if path_v.contains(&base_current) {
                return base_current;
            }
            
            match parent[current] {
                Some(p) => current = p,
                None => break,
            }
        }
        
        Self::find_base(v, base)
    }
    
    fn trace_and_update(
        start: usize,
        blossom_base: usize,
        base: &mut [usize],
        label: &mut [Option<u8>],
        parent: &[Option<usize>],
        mate: &[Option<usize>],
        queue: &mut VecDeque<usize>,
        in_queue: &mut HashSet<usize>,
    ) {
        let mut current = start;
        let mut visited = HashSet::new();
        
        loop {
            if visited.contains(&current) {
                break;
            }
            visited.insert(current);
            
            let current_base = Self::find_base(current, base);
            if current_base == blossom_base {
                break;
            }
            
            // Update base
            base[current] = blossom_base;
            
            // If this was an inner vertex, make it outer and add to queue
            if label[current] == Some(2) {
                label[current] = Some(1);
                if !in_queue.contains(&current) {
                    queue.push_back(current);
                    in_queue.insert(current);
                }
            }
            
            // Move to next vertex in the path
            if let Some(mate_current) = mate[current] {
                base[mate_current] = blossom_base;
                
                if let Some(next) = parent[mate_current] {
                    current = next;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
    }
    
    fn contract_blossom(
        blossom_base: usize,
        v: usize,
        w: usize,
        base: &mut [usize],
        label: &mut [Option<u8>],
        parent: &[Option<usize>],
        mate: &[Option<usize>],
        queue: &mut VecDeque<usize>,
        in_queue: &mut HashSet<usize>,
    ) {
        Self::trace_and_update(v, blossom_base, base, label, parent, mate, queue, in_queue);
        Self::trace_and_update(w, blossom_base, base, label, parent, mate, queue, in_queue);
    }
    
    fn build_path(v: usize, w: usize, parent: &[Option<usize>]) -> Vec<usize> {
        let mut path = vec![w, v];
        
        let mut current = v;
        while let Some(p) = parent[current] {
            path.push(p);
            current = p;
        }
        
        path
    }
    
    fn find_augmenting_path(&self, start: usize) -> Option<Vec<usize>> {
        let mut parent = vec![None; self.vertex_count];
        let mut base: Vec<usize> = (0..self.vertex_count).collect();
        let mut label = vec![None; self.vertex_count];
        let mut in_queue = HashSet::new();
        
        parent[start] = None;
        label[start] = Some(1); // Outer
        
        let mut queue = VecDeque::new();
        queue.push_back(start);
        in_queue.insert(start);
        
        while let Some(v) = queue.pop_front() {
            let v_base = Self::find_base(v, &base);
            
            for &w in &self.graph[v] {
                let w_base = Self::find_base(w, &base);
                
                // Skip if in same blossom
                if v_base == w_base {
                    continue;
                }
                
                // Case 1: w is unlabeled
                if label[w].is_none() {
                    if let Some(mate_w) = self.mate[w] {
                        // Add w (inner) and mate[w] (outer) to tree
                        label[w] = Some(2); // Inner
                        label[mate_w] = Some(1); // Outer
                        parent[w] = Some(v);
                        parent[mate_w] = Some(w);
                        
                        if !in_queue.contains(&mate_w) {
                            queue.push_back(mate_w);
                            in_queue.insert(mate_w);
                        }
                    } else {
                        // Found augmenting path!
                        return Some(Self::build_path(v, w, &parent));
                    }
                }
                // Case 2: w is outer (blossom detected)
                else if label[w] == Some(1) {
                    let blossom_base = Self::find_blossom_base(v, w, &parent, &base);
                    Self::contract_blossom(
                        blossom_base,
                        v,
                        w,
                        &mut base,
                        &mut label,
                        &parent,
                        &self.mate,
                        &mut queue,
                        &mut in_queue,
                    );
                }
            }
        }
        
        None
    }
    
    fn augment(&mut self, path: &[usize]) {
        for i in (0..path.len().saturating_sub(1)).step_by(2) {
            let u = path[i];
            let v = path[i + 1];
            self.mate[u] = Some(v);
            self.mate[v] = Some(u);
        }
    }
    
    fn maximum_matching(&mut self) -> Vec<(usize, usize)> {
        let mut improved = true;
        
        while improved {
            improved = false;
            
            for v in 0..self.vertex_count {
                if self.mate[v].is_none() {
                    if let Some(path) = self.find_augmenting_path(v) {
                        self.augment(&path);
                        improved = true;
                        break;
                    }
                }
            }
        }
        
        // Build matching set
        let mut matching = Vec::new();
        let mut seen = vec![false; self.vertex_count];
        
        for u in 0..self.vertex_count {
            if let Some(v) = self.mate[u] {
                if !seen[u] {
                    matching.push((u, v));
                    seen[u] = true;
                    seen[v] = true;
                }
            }
        }
        
        // Validate matching
        self.validate_matching(&matching);
        
        matching
    }
    
    fn validate_matching(&self, matching: &[(usize, usize)]) {
        let mut vertex_count_in_matching = vec![0; self.vertex_count];
        let mut errors = 0;
        
        eprintln!("\n=== Validation Report ===");
        eprintln!("Matching size (claimed): {}", matching.len());
        
        // Check 1: Each edge exists in the graph
        for &(u, v) in matching {
            let edge_exists = self.graph.get(u)
                .map(|neighbors| neighbors.contains(&v))
                .unwrap_or(false);
            
            if !edge_exists {
                eprintln!("ERROR: Edge ({}, {}) in matching but NOT in original graph!", u, v);
                errors += 1;
            }
            
            vertex_count_in_matching[u] += 1;
            vertex_count_in_matching[v] += 1;
        }
        
        // Check 2: No vertex appears in more than one edge
        for (vertex, &count) in vertex_count_in_matching.iter().enumerate() {
            if count > 1 {
                eprintln!("ERROR: Vertex {} appears in {} edges (should be at most 1)!", vertex, count);
                errors += 1;
            }
        }
        
        // Check 3: Explicit pairwise check - no two edges share a vertex
        for i in 0..matching.len() {
            for j in (i+1)..matching.len() {
                let (a, b) = matching[i];
                let (c, d) = matching[j];
                if a == c || a == d || b == c || b == d {
                    eprintln!("ERROR: Edges ({}, {}) and ({}, {}) share a vertex!", a, b, c, d);
                    errors += 1;
                }
            }
        }
        
        // Check 4: Count unique vertices
        let unique_vertices = vertex_count_in_matching.iter().filter(|&&c| c > 0).count();
        eprintln!("Number of edges in matching: {}", matching.len());
        eprintln!("Number of unique vertices: {}", unique_vertices);
        
        if errors > 0 {
            eprintln!("VALIDATION FAILED: {} errors found", errors);
        } else {
            eprintln!("VALIDATION PASSED: Matching is valid");
        }
        eprintln!("=========================\n");
    }
    
    #[allow(dead_code)]
    fn save_matching_to_file(&self, matching: &[(usize, usize)], filename: &str) {
        use std::io::Write;
        
        if let Ok(mut file) = std::fs::File::create(filename) {
            writeln!(file, "{}", matching.len()).ok();
            
            // Normalize edges (smaller vertex first) and sort
            let mut normalized_matching: Vec<(usize, usize)> = matching
                .iter()
                .map(|(u, v)| if u < v { (*u, *v) } else { (*v, *u) })
                .collect();
            normalized_matching.sort();
            
            for (u, v) in normalized_matching {
                writeln!(file, "{} {}", u, v).ok();
            }
        }
    }
}

fn load_graph_from_file(filename: &str) -> io::Result<(Vec<String>, Vec<(usize, usize)>)> {
    let file = File::open(filename)?;
    let mut lines = BufReader::new(file).lines();
    
    let first_line = lines.next().ok_or(io::Error::new(
        io::ErrorKind::InvalidData,
        "Empty file",
    ))??;
    
    let parts: Vec<usize> = first_line
        .split_whitespace()
        .map(|s| {
            s.parse().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "Invalid number format")
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    
    if parts.len() != 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "First line must have 2 numbers",
        ));
    }
    
    let (vertex_count, edge_count) = (parts[0], parts[1]);
    
    let vertices: Vec<String> = (0..vertex_count).map(|i| format!("V{:04}", i)).collect();
    
    let mut edges = Vec::new();
    for line in lines.take(edge_count) {
        let line = line?;
        let nums: Vec<usize> = line
            .split_whitespace()
            .map(|s| {
                s.parse().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "Invalid edge format")
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        
        if nums.len() != 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Each edge line must have 2 numbers",
            ));
        }
        
        edges.push((nums[0], nums[1]));
    }
    
    Ok((vertices, edges))
}

fn generate_test_graph(n: usize, edge_probability_percent: usize) -> (Vec<String>, Vec<(usize, usize)>) {
    let vertices: Vec<String> = (0..n).map(|i| format!("V{:04}", i)).collect();
    let mut edges = Vec::new();
    
    for i in 0..n {
        for j in (i + 1)..n {
            // Deterministic "random" based on i and j
            if ((i * 17 + j * 31) % 100) < edge_probability_percent {
                edges.push((i, j));
            }
        }
    }
    
    (vertices, edges)
}

fn run_example(vertices: &[String], edges: &[(usize, usize)], description: &str) {
    println!("{}", description);
    println!("Graph: {} vertices, {} edges", vertices.len(), edges.len());
    
    let start = Instant::now();
    let mut eb = EdmondsBlossomOptimized::new(vertices, edges);
    let matching = eb.maximum_matching();
    let duration = start.elapsed();
    
    println!("Matching size: {}", matching.len());
    if matching.len() <= 10 {
        print!("Matching: ");
        for (u, v) in &matching {
            print!("({},{}) ", eb.vertex_names[*u], eb.vertex_names[*v]);
        }
        println!();
    }
    println!("Execution time: {} ms", duration.as_millis());
    println!();
}

fn main() {
    println!("Edmonds' Blossom Algorithm (Optimized) - Rust Implementation");
    println!("==============================================================\n");
    
    let args: Vec<String> = env::args().collect();
    
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
        println!("Running built-in examples (use: ./edmonds_blossom_optimized_rust <filename> to load from file)\n");
        
        // Example 1: Triangle
        let vertices: Vec<String> = vec!["A", "B", "C", "D", "E"]
            .into_iter()
            .map(|s| s.to_string())
            .collect();
        let edges: Vec<(usize, usize)> = vec![
            (0, 1), (1, 2), (2, 0), // Triangle: A-B-C
            (2, 3), (3, 4),         // C-D-E
        ];
        
        run_example(&vertices, &edges, "Example 1: Graph with triangle (blossom)");
        
        // Example 2: Random graph
        let (vertices, edges) = generate_test_graph(20, 20);
        run_example(&vertices, &edges, "Example 2: Random graph (20 vertices)");
        
        // Benchmark
        println!("Benchmarking with larger graph...");
        let (vertices, edges) = generate_test_graph(100, 10);
        run_example(&vertices, &edges, "Benchmark: Random graph (100 vertices)");
        
        // Larger benchmark
        println!("Benchmarking with much larger graph...");
        let (vertices, edges) = generate_test_graph(500, 5);
        run_example(&vertices, &edges, "Benchmark: Random graph (500 vertices)");
    }
}

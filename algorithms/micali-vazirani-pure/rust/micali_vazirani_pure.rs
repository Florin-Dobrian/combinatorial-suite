/*
 * Micali-Vazirani Pure Algorithm - O(E√V) Maximum Matching - Rust Implementation
 *
 * True MV with DDFS, tenacity, regular + hanging bridges, petal contraction.
 * Ported from production Jorants MV-Matching-V2 (via C++ port).
 *
 * All integers, no hash containers, fully deterministic.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;

/* DDFS result codes */
const DDFS_EMPTY: i32 = 0;
const DDFS_PETAL: i32 = 1;
const DDFS_PATH: i32 = 2;

/* =========================================================================
 * Node
 * ========================================================================= */
struct Node {
    preds: Vec<i32>,
    pred_to: Vec<(usize, usize)>, /* (target, index in target's preds) */
    hanging_bridges: Vec<usize>,

    min_level: i32,
    max_level: i32,
    even_level: i32,
    odd_level: i32,
    match_: i32,
    bud: i32,
    above: i32,
    below: i32,
    ddfs_green: i32,
    ddfs_red: i32,
    number_preds: i32,
    deleted: bool,
    visited: bool,
}

impl Node {
    fn new() -> Self {
        Node {
            preds: Vec::new(),
            pred_to: Vec::new(),
            hanging_bridges: Vec::new(),
            min_level: NIL,
            max_level: NIL,
            even_level: NIL,
            odd_level: NIL,
            match_: NIL,
            bud: NIL,
            above: NIL,
            below: NIL,
            ddfs_green: NIL,
            ddfs_red: NIL,
            number_preds: 0,
            deleted: false,
            visited: false,
        }
    }

    fn set_min_level(&mut self, level: i32) {
        self.min_level = level;
        if level % 2 != 0 { self.odd_level = level; }
        else { self.even_level = level; }
    }

    fn set_max_level(&mut self, level: i32) {
        self.max_level = level;
        if level % 2 != 0 { self.odd_level = level; }
        else { self.even_level = level; }
    }

    fn outer(&self) -> bool {
        self.even_level != NIL && (self.odd_level == NIL || self.even_level < self.odd_level)
    }

    fn reset(&mut self) {
        self.preds.clear();
        self.pred_to.clear();
        self.hanging_bridges.clear();
        self.min_level = NIL;
        self.max_level = NIL;
        self.even_level = NIL;
        self.odd_level = NIL;
        self.bud = NIL;
        self.above = NIL;
        self.below = NIL;
        self.ddfs_green = NIL;
        self.ddfs_red = NIL;
        self.number_preds = 0;
        self.deleted = false;
        self.visited = false;
    }
}

/* =========================================================================
 * MVGraph — the full algorithm
 * ========================================================================= */
struct MVGraph {
    nodes: Vec<Node>,
    edges: Vec<usize>,       /* flat adjacency array (CSR) */
    adj_start: Vec<usize>,
    deg: Vec<usize>,

    levels: Vec<Vec<usize>>,
    bridges: Vec<Vec<(usize, usize)>>,

    green_stack: Vec<(i32, usize)>,
    red_stack: Vec<(i32, usize)>,
    path_found: Vec<usize>,
    ddfs_nodes_seen: Vec<usize>,
    ddfs_bottleneck: i32,

    matchnum: usize,
    bridgenum: i32,
    todonum: i32,
}

impl MVGraph {
    fn new() -> Self {
        MVGraph {
            nodes: Vec::new(),
            edges: Vec::new(),
            adj_start: Vec::new(),
            deg: Vec::new(),
            levels: Vec::new(),
            bridges: Vec::new(),
            green_stack: Vec::new(),
            red_stack: Vec::new(),
            path_found: Vec::new(),
            ddfs_nodes_seen: Vec::new(),
            ddfs_bottleneck: NIL,
            matchnum: 0,
            bridgenum: 0,
            todonum: 0,
        }
    }

    /* ---- construction ---- */
    fn build(&mut self, n: usize, edge_list: &[(usize, usize)]) {
        self.nodes = (0..n).map(|_| Node::new()).collect();
        let mut adj = vec![Vec::new(); n];
        for &(u, v) in edge_list {
            if u < n && v < n && u != v {
                adj[u].push(v);
                adj[v].push(u);
            }
        }
        for a in &mut adj { a.sort_unstable(); }

        self.adj_start = vec![0; n];
        self.deg = vec![0; n];
        self.edges.clear();
        for i in 0..n {
            self.adj_start[i] = self.edges.len();
            self.deg[i] = adj[i].len();
            for &nb in &adj[i] { self.edges.push(nb); }
        }
    }

    /* ---- greedy initialization ---- */
    fn greedy_init(&mut self) {
        let n = self.nodes.len();
        for j in 0..n {
            if self.nodes[j].match_ == NIL {
                for k in 0..self.deg[j] {
                    let i = self.edges[self.adj_start[j] + k];
                    if self.nodes[i].match_ == NIL {
                        self.nodes[j].match_ = i as i32;
                        self.nodes[i].match_ = j as i32;
                        self.matchnum += 1;
                        break;
                    }
                }
            }
        }
    }

    /* ---- helpers ---- */
    fn add_to_level(&mut self, level: usize, node: usize) {
        while self.levels.len() <= level { self.levels.push(Vec::new()); }
        self.levels[level].push(node);
        self.todonum += 1;
    }

    fn add_to_bridges(&mut self, level: usize, n1: usize, n2: usize) {
        while self.bridges.len() <= level { self.bridges.push(Vec::new()); }
        self.bridges[level].push((n1, n2));
        self.bridgenum += 1;
    }

    fn tenacity(&self, n1: usize, n2: usize) -> i32 {
        if self.nodes[n1].match_ == n2 as i32 {
            /* matched bridge */
            if self.nodes[n1].odd_level != NIL && self.nodes[n2].odd_level != NIL {
                return self.nodes[n1].odd_level + self.nodes[n2].odd_level + 1;
            }
        } else {
            /* unmatched bridge */
            if self.nodes[n1].even_level != NIL && self.nodes[n2].even_level != NIL {
                return self.nodes[n1].even_level + self.nodes[n2].even_level + 1;
            }
        }
        NIL
    }

    fn bud_star(&self, c: usize) -> usize {
        let b = self.nodes[c].bud;
        if b == NIL { c } else { self.bud_star(b as usize) }
    }

    fn bud_star_includes(&self, c: usize, goal: usize) -> bool {
        if c == goal { return true; }
        let b = self.nodes[c].bud;
        if b == NIL { return false; }
        self.bud_star_includes(b as usize, goal)
    }

    /* ---- reset between phases ---- */
    fn phase_reset(&mut self) {
        for lv in &mut self.levels { lv.clear(); }
        for br in &mut self.bridges { br.clear(); }
        self.bridgenum = 0;
        self.todonum = 0;
        let n = self.nodes.len();
        for i in 0..n {
            self.nodes[i].reset();
            if self.nodes[i].match_ == NIL {
                self.add_to_level(0, i);
                self.nodes[i].set_min_level(0);
            }
        }
    }

    /* ---- step_to: core level-building step ---- */
    fn step_to(&mut self, to: usize, from: usize, level: i32) {
        let next = level + 1;
        let tl = self.nodes[to].min_level;
        if tl == NIL || tl >= next {
            if tl != next {
                self.add_to_level(next as usize, to);
                self.nodes[to].set_min_level(next);
            }
            self.nodes[to].preds.push(from as i32);
            self.nodes[to].number_preds += 1;
            let idx = self.nodes[to].preds.len() - 1;
            self.nodes[from].pred_to.push((to, idx));
        } else {
            /* found a bridge */
            let ten = self.tenacity(to, from);
            if ten == NIL {
                self.nodes[to].hanging_bridges.push(from);
                self.nodes[from].hanging_bridges.push(to);
            } else {
                self.add_to_bridges(((ten - 1) / 2) as usize, to, from);
            }
        }
    }

    /* ---- MIN phase ---- */
    fn min_phase(&mut self, i: usize) {
        if i >= self.levels.len() { return; }
        let level_snap = self.levels[i].clone();
        for &current in &level_snap {
            self.todonum -= 1;
            let match_ = self.nodes[current].match_;
            if i % 2 == 0 {
                /* even level — explore non-matching edges */
                let start = self.adj_start[current];
                let d = self.deg[current];
                for j in 0..d {
                    let edge = self.edges[start + j];
                    if edge as i32 != match_ {
                        self.step_to(edge, current, i as i32);
                    }
                }
            } else {
                /* odd level — follow matching edge only */
                if match_ != NIL {
                    self.step_to(match_ as usize, current, i as i32);
                }
            }
        }
    }

    /* ---- MAX phase ---- */
    fn max_phase(&mut self, i: usize) -> bool {
        let mut found = false;
        if i >= self.bridges.len() { return false; }

        let bridge_snap = self.bridges[i].clone();
        for &(n1, n2) in &bridge_snap {
            self.bridgenum -= 1;
            if self.nodes[n1].deleted || self.nodes[n2].deleted { continue; }

            let result = self.ddfs(n1, n2);
            if result == DDFS_EMPTY { continue; }

            if result == DDFS_PATH {
                self.find_path(n1, n2);
                self.augment_path();
                if self.nodes.len() / 2 <= self.matchnum { return true; }
                self.remove_path();
                found = true;
            } else {
                /* DDFS_PETAL */
                let b = self.ddfs_bottleneck as usize;
                let current_ten = (i * 2 + 1) as i32;
                let seen = self.ddfs_nodes_seen.clone();
                for &itt in &seen {
                    self.nodes[itt].bud = b as i32;
                    let ml = self.nodes[itt].min_level;
                    self.nodes[itt].set_max_level(current_ten - ml);
                    let max_lv = self.nodes[itt].max_level as usize;
                    self.add_to_level(max_lv, itt);
                    let hangs = self.nodes[itt].hanging_bridges.clone();
                    for &hanging in &hangs {
                        let hanging_ten = self.tenacity(itt, hanging);
                        if hanging_ten != NIL {
                            self.add_to_bridges(((hanging_ten - 1) / 2) as usize, itt, hanging);
                        }
                    }
                }
            }
        }
        found
    }

    /* ==================================================================
     * DDFS — Double Depth-First Search
     * ================================================================== */

    fn add_pred_to_stack(preds: &[i32], cur_node: usize, stack: &mut Vec<(i32, usize)>) {
        for &pred in preds {
            if pred != NIL {
                stack.push((cur_node as i32, pred as usize));
            }
        }
    }

    fn edge_valid(e: (i32, i32)) -> bool {
        !(e.0 == NIL && e.1 == NIL)
    }

    fn stack_pop(stack: &mut Vec<(i32, usize)>) -> (i32, i32) {
        if let Some((a, b)) = stack.pop() {
            (a, b as i32)
        } else {
            (NIL, NIL)
        }
    }

    fn prepare_next(&mut self, nx: &mut (i32, i32)) {
        if nx.0 != NIL {
            self.nodes[nx.0 as usize].below = nx.1;
        }
        nx.1 = self.bud_star(nx.1 as usize) as i32;
    }

    fn level_of(&self, nx: &(i32, i32)) -> i32 {
        let n = self.bud_star(nx.1 as usize);
        self.nodes[n].min_level
    }

    fn step_into(&mut self, c: &mut i32, nx: &mut (i32, i32),
                 stack: &mut Vec<(i32, usize)>, green_top: usize, red_top: usize) {
        self.prepare_next(nx);
        let nx_second = nx.1 as usize;
        if !self.nodes[nx_second].visited {
            self.nodes[nx_second].above = nx.0;
            *c = nx_second as i32;
            self.nodes[nx_second].visited = true;
            self.nodes[nx_second].ddfs_green = green_top as i32;
            self.nodes[nx_second].ddfs_red = red_top as i32;
            self.ddfs_nodes_seen.push(nx_second);
            let preds = self.nodes[nx_second].preds.clone();
            Self::add_pred_to_stack(&preds, nx_second, stack);
        }
        let popped = Self::stack_pop(stack);
        nx.0 = popped.0;
        nx.1 = popped.1;
    }

    fn ddfs(&mut self, green_top: usize, red_top: usize) -> i32 {
        self.ddfs_nodes_seen.clear();
        self.ddfs_bottleneck = NIL;

        self.red_stack.clear();
        self.green_stack.clear();

        let mut g: i32 = NIL;
        let mut r: i32 = NIL;

        if self.bud_star(red_top) == self.bud_star(green_top) { return DDFS_EMPTY; }
        if self.nodes[green_top].min_level == 0 && self.nodes[red_top].min_level == 0 {
            return DDFS_PATH;
        }

        let mut ng: (i32, i32) = (NIL, green_top as i32);
        let mut nr: (i32, i32) = (NIL, red_top as i32);
        let mut red_before: (i32, i32) = (NIL, NIL);
        let mut green_before: (i32, i32) = (NIL, NIL);

        loop {
            /* check termination */
            if !(r == NIL || g == NIL ||
                 self.nodes[r as usize].min_level > 0 ||
                 self.nodes[g as usize].min_level > 0) {
                break;
            }

            /* balance levels */
            while Self::edge_valid((nr.0, nr.1)) && Self::edge_valid((ng.0, ng.1))
                  && self.level_of(&nr) != self.level_of(&ng) {

                while Self::edge_valid((nr.0, nr.1))
                      && self.level_of(&nr) > self.level_of(&ng) {
                    let mut sr = std::mem::take(&mut self.red_stack);
                    self.step_into(&mut r, &mut nr, &mut sr, green_top, red_top);
                    self.red_stack = sr;
                }

                if !Self::edge_valid((nr.0, nr.1)) {
                    nr = red_before;
                    let mut tmp = red_before.0;
                    while self.nodes[tmp as usize].above != NIL {
                        let rc = self.nodes[tmp as usize].above as usize;
                        let preds = self.nodes[rc].preds.clone();
                        for &ri in &preds {
                            if ri == NIL { continue; }
                            if self.bud_star(ri as usize) == tmp as usize {
                                self.nodes[rc].below = ri;
                                break;
                            }
                        }
                        tmp = self.nodes[tmp as usize].above;
                    }
                }

                while Self::edge_valid((ng.0, ng.1))
                      && self.level_of(&nr) < self.level_of(&ng) {
                    let mut sg = std::mem::take(&mut self.green_stack);
                    self.step_into(&mut g, &mut ng, &mut sg, green_top, red_top);
                    self.green_stack = sg;
                }

                if !Self::edge_valid((ng.0, ng.1)) {
                    ng = green_before;
                    let mut tmp = green_before.0;
                    while self.nodes[tmp as usize].above != NIL {
                        let rc = self.nodes[tmp as usize].above as usize;
                        let preds = self.nodes[rc].preds.clone();
                        for &ri in &preds {
                            if ri == NIL { continue; }
                            if self.bud_star(ri as usize) == tmp as usize {
                                self.nodes[rc].below = ri;
                                break;
                            }
                        }
                        tmp = self.nodes[tmp as usize].above;
                    }
                }
            }

            /* check collision */
            if self.bud_star(nr.1 as usize) == self.bud_star(ng.1 as usize) {
                if !self.red_stack.is_empty() {
                    red_before = nr;
                    self.prepare_next(&mut nr);
                    let popped = Self::stack_pop(&mut self.red_stack);
                    nr.0 = popped.0;
                    nr.1 = popped.1;
                    if Self::edge_valid((nr.0, nr.1)) { r = nr.0; }
                    else { nr = red_before; }
                } else if !self.green_stack.is_empty() {
                    green_before = ng;
                    self.prepare_next(&mut ng);
                    let popped = Self::stack_pop(&mut self.green_stack);
                    ng.0 = popped.0;
                    ng.1 = popped.1;
                    if Self::edge_valid((ng.0, ng.1)) { g = ng.0; }
                    else { ng = green_before; }
                } else {
                    self.prepare_next(&mut nr);
                    self.prepare_next(&mut ng);
                    self.ddfs_bottleneck = nr.1;
                    return DDFS_PETAL;
                }
            } else {
                /* step both sides */
                let mut sr = std::mem::take(&mut self.red_stack);
                self.step_into(&mut r, &mut nr, &mut sr, green_top, red_top);
                self.red_stack = sr;

                let mut sg = std::mem::take(&mut self.green_stack);
                self.step_into(&mut g, &mut ng, &mut sg, green_top, red_top);
                self.green_stack = sg;
            }
        }
        DDFS_PATH
    }

    /* ==================================================================
     * Path finding and augmentation
     * ================================================================== */

    fn find_path(&mut self, n1: usize, n2: usize) {
        self.path_found.clear();
        self.walk_down_path(n1);
        self.path_found.reverse();
        self.walk_down_path(n2);
    }

    fn walk_down_path(&mut self, start: usize) {
        let mut cur = start as i32;
        while cur != NIL {
            let c = cur as usize;
            if self.nodes[c].bud != NIL {
                cur = self.walk_blossom(c) as i32;
            } else {
                self.path_found.push(c);
                cur = self.nodes[c].below;
            }
        }
    }

    fn jump_bridge(&mut self, cur: usize) -> usize {
        let dg = self.nodes[cur].ddfs_green;
        let dr = self.nodes[cur].ddfs_red;

        if dg == cur as i32 { return dr as usize; }
        if dr == cur as i32 { return dg as usize; }

        if self.bud_star_includes(dg as usize, cur) {
            let before = self.path_found.len();
            let mut b = dg as usize;
            while b != cur { b = self.walk_blossom(b); }
            self.path_found[before..].reverse();
            return dr as usize;
        } else {
            let before = self.path_found.len();
            let mut b = dr as usize;
            while b != cur { b = self.walk_blossom(b); }
            self.path_found[before..].reverse();
            return dg as usize;
        }
    }

    fn walk_blossom(&mut self, cur: usize) -> usize {
        if self.nodes[cur].outer() {
            self.walk_blossom_down(cur, NIL)
        } else {
            let top = self.walk_blossom_up(cur);
            let before = top;
            let jumped = self.jump_bridge(top);
            self.walk_blossom_down(jumped, before as i32)
        }
    }

    fn walk_blossom_down(&mut self, cur: usize, before: i32) -> usize {
        let before_node = if before == NIL { cur } else { before as usize };
        let b = self.nodes[cur].bud;
        let mut c = cur as i32;
        while c != NIL && c != b {
            let cu = c as usize;
            if self.nodes[cu].ddfs_green != self.nodes[before_node].ddfs_green ||
               self.nodes[cu].ddfs_red != self.nodes[before_node].ddfs_red {
                c = self.walk_blossom(cu) as i32;
            } else {
                self.path_found.push(cu);
                c = self.nodes[cu].below;
            }
        }
        c as usize
    }

    fn walk_blossom_up(&mut self, cur: usize) -> usize {
        let mut c = cur;
        loop {
            self.path_found.push(c);
            if self.nodes[c].above == NIL { break; }
            let b_node = self.nodes[c].above as usize;
            let b = self.nodes[b_node].below;
            if b != c as i32 && self.bud_star_includes(b as usize, c) {
                let before = self.path_found.len();
                let mut bb = b as usize;
                while bb != c { bb = self.walk_blossom(bb); }
                self.path_found[before..].reverse();
            }
            c = self.nodes[c].above as usize;
        }
        c
    }

    fn augment_path(&mut self) {
        let mut i = 0;
        while i + 1 < self.path_found.len() {
            let n1 = self.path_found[i];
            let n2 = self.path_found[i + 1];
            self.nodes[n1].match_ = n2 as i32;
            self.nodes[n2].match_ = n1 as i32;
            i += 2;
        }
        self.matchnum += 1;
    }

    fn remove_path(&mut self) {
        while let Some(current) = self.path_found.pop() {
            if !self.nodes[current].deleted {
                self.nodes[current].deleted = true;
                let pred_to = self.nodes[current].pred_to.clone();
                for &(target, idx) in &pred_to {
                    if !self.nodes[target].deleted {
                        self.nodes[target].preds[idx] = NIL;
                        self.nodes[target].number_preds -= 1;
                        if self.nodes[target].number_preds <= 0 {
                            self.path_found.push(target);
                        }
                    }
                }
            }
        }
    }

    /* ---- main matching driver ---- */
    fn max_match(&mut self) {
        let n = self.nodes.len();
        for i in 0..n {
            if self.nodes[i].match_ == NIL {
                self.add_to_level(0, i);
                self.nodes[i].set_min_level(0);
            }
        }
        let mut found = self.max_match_phase();
        while self.nodes.len() / 2 > self.matchnum && found {
            self.phase_reset();
            found = self.max_match_phase();
        }
    }

    fn max_match_phase(&mut self) -> bool {
        let n = self.nodes.len();
        let mut found = false;
        for i in 0..(n / 2 + 1) {
            if !found {
                if self.todonum <= 0 && self.bridgenum <= 0 { return false; }
                self.min_phase(i);
                found = self.max_phase(i);
            }
        }
        found
    }

    fn get_matching(&self) -> Vec<(usize, usize)> {
        let mut result = Vec::new();
        for i in 0..self.nodes.len() {
            if self.nodes[i].match_ != NIL && (self.nodes[i].match_ as usize) > i {
                result.push((i, self.nodes[i].match_ as usize));
            }
        }
        result
    }
}

/* =========================================================================
 * File I/O, validation, and main
 * ========================================================================= */

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

fn validate_matching(n: usize, matching: &[(usize, usize)]) {
    let mut deg = vec![0i32; n];
    let mut errors = 0;

    for &(u, v) in matching {
        deg[u] += 1;
        deg[v] += 1;
    }
    for i in 0..n {
        if deg[i] > 1 {
            eprintln!("ERROR: Vertex {} in {} edges!", i, deg[i]);
            errors += 1;
        }
    }
    let matched = deg.iter().filter(|&&d| d > 0).count();

    println!("\n=== Validation Report ===");
    println!("Matching size: {}", matching.len());
    println!("Matched vertices: {}", matched);
    println!("{}", if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" });
    println!("=========================\n");
}

fn main() {
    println!("Micali-Vazirani Pure Algorithm - Rust Implementation");
    println!("=====================================================");
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
            let mut mv = MVGraph::new();
            mv.build(n, &edges);
            mv.greedy_init();
            mv.max_match();
            let duration = start.elapsed();

            let matching = mv.get_matching();
            validate_matching(n, &matching);

            println!("Matching size: {}", matching.len());
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

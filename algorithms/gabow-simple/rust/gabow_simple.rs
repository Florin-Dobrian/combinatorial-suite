/*
 * Gabow's Algorithm (Simple) - O(V * E) Maximum Matching
 *
 * Faithful to Gabow 1976: forest BFS with blossom contraction via
 * union-find. No physical contraction -- bases are tracked virtually.
 * Epoch-based interleaved LCA, path-only contraction, bridge recording
 * for augmentation through blossoms.
 *
 * Forest search: each iteration labels ALL free vertices as EVEN roots
 * simultaneously and grows a search forest. An augmenting path is found
 * when two different trees meet (EVEN-EVEN edge across trees, detected
 * by find_lca returning NIL). One augmentation per iteration, then full
 * reset and repeat until no augmenting path exists.
 *
 * Complexity: O(V * E) -- each iteration does O(E) work, at most V/2
 * augmentations total.
 *
 * Rust implementation -- fully deterministic, no hash containers.
 */

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const NIL: i32 = -1;
const UNLABELED: i32 = 0;
const EVEN: i32 = 1;
const ODD: i32 = 2;

struct GabowSimple {
    n: usize,
    greedy_size: usize,
    graph: Vec<Vec<usize>>,
    mate: Vec<i32>,
    base: Vec<usize>,
    parent: Vec<i32>,
    label: Vec<i32>,

    /* Bridge recording for ODD vertices absorbed into blossoms */
    bridge_src: Vec<i32>,
    bridge_tgt: Vec<i32>,

    /* Epoch-based interleaved LCA */
    lca_tag1: Vec<usize>,
    lca_tag2: Vec<usize>,
    lca_epoch: usize,
}

impl GabowSimple {
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
            adj.dedup();
        }

        GabowSimple {
            n,
            graph,
            mate: vec![NIL; n],
            base: vec![0; n],
            parent: vec![NIL; n],
            label: vec![UNLABELED; n],
            bridge_src: vec![NIL; n],
            bridge_tgt: vec![NIL; n],
            lca_tag1: vec![0; n],
            lca_tag2: vec![0; n],
            lca_epoch: 0,
            greedy_size: 0,
        }
    }

    fn greedy_init(&mut self) -> usize {
        let mut cnt = 0usize;
        for u in 0..self.n {
            if self.mate[u] != NIL { continue; }
            let neighbors: Vec<usize> = self.graph[u].clone();
            for &v in &neighbors {
                if self.mate[v] == NIL {
                    self.mate[u] = v as i32;
                    self.mate[v] = u as i32;
                    cnt += 1;
                    break;
                }
            }
        }
        cnt
    }

    fn greedy_init_md(&mut self) -> usize {
        let mut cnt = 0usize;
        let mut deg = vec![0usize; self.n];
        for u in 0..self.n {
            for &v in &self.graph[u] {
                deg[v] += 1;
            }
        }
        let mut order: Vec<usize> = (0..self.n).collect();
        order.sort_unstable_by(|&a, &b| deg[a].cmp(&deg[b]).then(a.cmp(&b)));
        for u in order {
            if self.mate[u] != NIL { continue; }
            let mut best: i32 = -1;
            let mut best_deg = usize::MAX;
            let neighbors: Vec<usize> = self.graph[u].clone();
            for &v in &neighbors {
                if self.mate[v] == NIL && deg[v] < best_deg {
                    best = v as i32;
                    best_deg = deg[v];
                }
            }
            if best >= 0 {
                self.mate[u] = best;
                self.mate[best as usize] = u as i32;
                cnt += 1;
            }
        }
        cnt
    }

    /* Path-halving find for union-find base */
    fn find_base(&mut self, mut v: usize) -> usize {
        while self.base[v] != v {
            self.base[v] = self.base[self.base[v]];
            v = self.base[v];
        }
        v
    }

    /* Interleaved LCA using epoch tags.
     * Returns the LCA base if u and v are in the same tree, or NIL if
     * they are in different trees (= augmenting path). */
    fn find_lca(&mut self, u: usize, v: usize) -> i32 {
        self.lca_epoch += 1;
        let ep = self.lca_epoch;
        let mut hx = self.find_base(u);
        let mut hy = self.find_base(v);
        self.lca_tag1[hx] = ep;
        self.lca_tag2[hy] = ep;
        loop {
            if self.lca_tag1[hy] == ep { return hy as i32; }
            if self.lca_tag2[hx] == ep { return hx as i32; }
            let hxr = self.mate[hx] == NIL;
            let hyr = self.mate[hy] == NIL;
            if hxr && hyr { return NIL; }
            if !hxr {
                hx = self.find_base(self.parent[self.mate[hx] as usize] as usize);
                self.lca_tag1[hx] = ep;
            }
            if !hyr {
                hy = self.find_base(self.parent[self.mate[hy] as usize] as usize);
                self.lca_tag2[hy] = ep;
            }
        }
    }

    /* Path-only contraction: walk from x back to lca, union bases.
     * For each ODD vertex mv on the path, record the bridge (x, y)
     * and enqueue mv as newly-EVEN if it wasn't already. */
    fn shrink_path(
        &mut self,
        lca: usize,
        x: usize,
        y: usize,
        queue: &mut Vec<usize>,
    ) {
        let mut v = self.find_base(x);
        while v != lca {
            let mv = self.mate[v] as usize;
            /* Union both v and mv into lca's component */
            let fv = self.find_base(v);
            self.base[fv] = lca;
            let fmv = self.find_base(mv);
            self.base[fmv] = lca;
            self.base[lca] = lca;

            /* Record bridge for mv */
            self.bridge_src[mv] = x as i32;
            self.bridge_tgt[mv] = y as i32;

            /* If mv was ODD and not yet enqueued as EVEN, enqueue it */
            if self.label[mv] != EVEN {
                self.label[mv] = EVEN;
                queue.push(mv);
            }

            /* Walk up */
            v = self.find_base(self.parent[mv] as usize);
        }
    }

    /* Trace from vertex v to vertex u (or to a root if u==NIL),
     * collecting edge pairs for augmentation.
     *   - No bridge -> "originally EVEN": step mate -> parent
     *   - Has bridge -> "originally ODD, absorbed into blossom":
     *     recurse through bridge */
    fn trace_path(&self, v: i32, u: i32, pairs: &mut Vec<(i32, i32)>) {
        struct Frame {
            v: i32,
            u: i32,
            phase: i32,
            sb: i32,
            tb: i32,
        }
        let mut stk: Vec<Frame> = vec![Frame { v, u, phase: 0, sb: 0, tb: 0 }];

        while !stk.is_empty() {
            let len = stk.len();
            let f = &mut stk[len - 1];

            if f.v == f.u {
                stk.pop();
                continue;
            }

            if f.phase == 0 {
                let fv = f.v as usize;
                if self.bridge_src[fv] == NIL {
                    /* Originally EVEN vertex (no bridge) */
                    if self.mate[fv] == NIL {
                        /* Root (free vertex) -- done */
                        stk.pop();
                        continue;
                    }
                    let mv = self.mate[fv];
                    let pmv = self.parent[mv as usize];
                    pairs.push((mv, pmv));
                    stk[len - 1].v = pmv;
                    continue;
                }
                /* Has bridge -- originally ODD, absorbed into blossom */
                let sb = self.bridge_src[fv];
                let tb = self.bridge_tgt[fv];
                let mate_fv = self.mate[fv];
                stk[len - 1].sb = sb;
                stk[len - 1].tb = tb;
                stk[len - 1].phase = 1;
                stk.push(Frame { v: sb, u: mate_fv, phase: 0, sb: 0, tb: 0 });
                continue;
            }
            if f.phase == 1 {
                let sb = f.sb;
                let tb = f.tb;
                let fu = f.u;
                pairs.push((sb, tb));
                stk[len - 1].phase = 2;
                stk.push(Frame { v: tb, u: fu, phase: 0, sb: 0, tb: 0 });
                continue;
            }
            stk.pop();
        }
    }

    /* Augment along the path:
     *   root_u ~~~ u -- v ~~~ root_v
     * Collect all edge pairs, then flip mate for all of them. */
    fn augment_two_sides(&mut self, u: usize, v: usize) {
        let mut pairs: Vec<(i32, i32)> = Vec::new();
        /* The cross-tree edge */
        pairs.push((u as i32, v as i32));
        /* Trace from u to its root */
        self.trace_path(u as i32, NIL, &mut pairs);
        /* Trace from v to its root */
        self.trace_path(v as i32, NIL, &mut pairs);
        /* Flip all */
        for &(a, b) in &pairs {
            self.mate[a as usize] = b;
            self.mate[b as usize] = a;
        }
    }

    /* Find one augmenting path in the forest and augment.
     * Returns true if an augmentation was performed. */
    fn find_and_augment(&mut self) -> bool {
        /* Reset per-iteration state */
        for i in 0..self.n {
            self.base[i] = i;
            self.parent[i] = NIL;
            self.label[i] = UNLABELED;
            self.bridge_src[i] = NIL;
            self.bridge_tgt[i] = NIL;
        }

        let mut queue: Vec<usize> = Vec::with_capacity(self.n);

        /* All free vertices become EVEN roots */
        for v in 0..self.n {
            if self.mate[v] == NIL {
                self.label[v] = EVEN;
                queue.push(v);
            }
        }

        let mut qi = 0;
        while qi < queue.len() {
            let u = queue[qi];
            qi += 1;

            /* Check that u is still effectively EVEN */
            let bu = self.find_base(u);
            if self.label[bu] != EVEN { continue; }

            let neighbors = self.graph[u].clone();
            for &v in &neighbors {
                let bu2 = self.find_base(u);
                let bv = self.find_base(v);
                if bu2 == bv { continue; }
                if v as i32 == self.mate[u] { continue; }

                if self.label[bv] == UNLABELED {
                    /* v is matched and unlabeled -> grow step */
                    self.label[v] = ODD;
                    self.parent[v] = u as i32;
                    let w = self.mate[v] as usize;
                    self.label[w] = EVEN;
                    queue.push(w);

                } else if self.label[bv] == EVEN {
                    /* EVEN-EVEN edge: blossom or augmenting path */
                    let lca = self.find_lca(u, v);
                    if lca != NIL {
                        /* Same tree -> blossom contraction */
                        let lca_u = lca as usize;
                        self.shrink_path(lca_u, u, v, &mut queue);
                        self.shrink_path(lca_u, v, u, &mut queue);
                    } else {
                        /* Different trees -> augmenting path! */
                        self.augment_two_sides(u, v);
                        return true;
                    }
                }
                /* label[bv] == ODD: ignore */
            }
        }
        false
    }

    fn maximum_matching(&mut self, greedy_mode: i32) -> Vec<(usize, usize)> {
        self.greedy_size = match greedy_mode {
            1 => self.greedy_init(),
            2 => self.greedy_init_md(),
            _ => 0,
        };

        while self.find_and_augment() {}

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
        if graph[u].binary_search(&v).is_err() {
            eprintln!("ERROR: Edge ({}, {}) not in graph!", u, v);
            errors += 1;
        }
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
    println!(
        "{}",
        if errors > 0 { "VALIDATION FAILED" } else { "VALIDATION PASSED" }
    );
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
    println!("Gabow's Algorithm (Simple) - Rust Implementation");
    println!("==================================================\n");

    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <filename> [--greedy|--greedy-md]", args[0]);
        std::process::exit(1);
    }

    let greedy_mode: i32 = if args.iter().any(|a| a == "--greedy-md") {
        2
    } else if args.iter().any(|a| a == "--greedy") {
        1
    } else {
        0
    };

    match load_graph(&args[1]) {
        Ok((n, edges)) => {
            println!("Graph: {} vertices, {} edges", n, edges.len());
            let start = Instant::now();
            let mut gabow = GabowSimple::new(n, &edges);
            let matching = gabow.maximum_matching(greedy_mode);
            let duration = start.elapsed();
            validate_matching(n, &gabow.graph, &matching);
            println!("Matching size: {}", matching.len());
            if greedy_mode > 0 {
                let gs = gabow.greedy_size;
                let fs = matching.len();
                println!("Greedy init size: {}", gs);
                if fs > 0 {
                    println!("Greedy/Final: {:.2}%", 100.0 * gs as f64 / fs as f64);
                } else {
                    println!("Greedy/Final: NA");
                }
            }
            println!("Time: {} ms", duration.as_millis());
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

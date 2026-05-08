// Tentative: single-source general-graph MWM (Edmonds-style with weighted blossoms).
// Sparse, flat CSR, single-source SDAP-with-potentials.
// Real-valued weights, virtual blossom contraction (G1-style), single PQ.
//
// This is the analog of hungarian_mwm_csr.cpp but for general (not bipartite)
// graphs.  The structural differences from V1:
//
//   1. No S/T partition.  All vertices are equally "source-side."
//      mate[] is a single array; reduced-cost duals y[v] are per-vertex.
//   2. Blossoms (odd cycles in the alternating tree) are handled via virtual
//      contraction: a union-find structure base_par[] maps each vertex to the
//      base of its top-level blossom.  Sub-blossoms stay where they are in
//      the graph; find_base() resolves labels and duals.
//   3. The PQ holds edge IDs (not vertex IDs).  See the work log subsection
//      "Why bipartite uses vertices and general uses edges."
//   4. Three event types must be tracked, all encoded as edge events plus a
//      blossom event:
//        - GROW: edge from outer (S-blossom) to UNLABELED becomes tight.
//        - AUGMENT or BLOSSOM: edge between two outer vertices becomes tight.
//          Augment if in different trees (we have only one tree here, so it
//          becomes the same as the path-to-free-vertex check) or blossom
//          if in the same tree.
//        - EXPAND: ODD blossom's z(B) reaches 0.  ODD-only because EVEN
//          blossoms have z growing (no upper bound, no event); only ODD
//          blossom z shrinks toward 0 and forces an expand to keep z >= 0.
//      Edges in PQ for GROW/AUGMENT/BLOSSOM; ODD blossoms in PQ for EXPAND.
//      EVEN blossoms never enter the PQ.
//   5. Single binary-heap PQ, heterogeneous entries with type tag.
//
// TRUST MODE: this implementation is "trust mode" — we return only M, no
// LP-witness verification structures.  Blossoms are per-iteration scratch
// state; they form during a search, may go through expand events, and
// dissolve at iteration end with z(B) redistributed to sub-blossom y
// values.  No surviving blossoms across iterations.  This trades the LP-
// witness machinery (BT laminar storage, b membership map) for code
// simplicity, while preserving correctness of the matching itself.

#include <vector>
#include <queue>
#include <limits>

using namespace std;

// Flat CSR for general graph.  Each undirected edge has a unique ID and
// appears in the adjacency list of both endpoints.  esrc[e]/etgt[e] give
// stable endpoints; opposite(v, e) returns the other endpoint.

struct Graph {
    int numVtxs;
    vector<int> off;          // size numVtxs + 1
    vector<int> adj;          // size 2*numEdges, edge IDs (each edge appears twice)
    vector<int> esrc, etgt;   // size numEdges, edge endpoints
    vector<double> wght;      // size numEdges
};

struct Matching {
    vector<int> mate;         // mate[v] = matched partner, or -1 if unmatched
    vector<double> yDual;     // y(v), per-vertex dual (real-valued)
    // zDual[B] is the blossom's dual (0 for trivial blossoms = vertices).
    // In TRUST MODE (no LP-witness verification needed), zDual is per-
    // iteration scratch state: blossoms form during a search, may
    // accumulate z values, then dissolve at iteration end with z
    // redistributed back to sub-blossom y values.  Across iterations,
    // zDual is empty.  We carry it in the struct here for in-iteration
    // state, but in a stricter trust-mode design it could be local
    // scratch in `augment` rather than a returned field.
    //
    // For verification mode (LP-witness retention), this would also need
    // to be paired with a laminar membership map and survive across
    // iterations.  See work log: "LP-witness for general-graph MWM" and
    // "Trust mode: skip the LP-witness machinery."
    vector<double> zDual;
};

static constexpr double INF = numeric_limits<double>::infinity();

enum Label { UNLABELED = 0, EVEN = 1, ODD = 2 };

// PQ entry types.
enum EntryType { EDGE_EVENT = 0, BLOSSOM_EVENT = 1 };

struct PQEntry {
    double eventTime;     // global Δ at which this event fires
    int id;               // edge ID for EDGE_EVENT, blossom ID for BLOSSOM_EVENT
    EntryType type;
    bool operator>(const PQEntry& o) const { return eventTime > o.eventTime; }
};

bool augment(const Graph& graph, int sFirst, Matching& matching) {
    // Per-iteration scratch state.  Allocated fresh; production would reuse.
    int n = graph.numVtxs;

    vector<int> baseParent(n, 0);             // union-find for blossom membership
    for (int v = 0; v < n; ++v) baseParent[v] = v;
    auto findBase = [&](int v) -> int {
        while (baseParent[v] != v) {
            baseParent[v] = baseParent[baseParent[v]];   // path compression
            v = baseParent[v];
        }
        return v;
    };

    vector<Label> label(n, UNLABELED);         // label of each top-level blossom
    vector<int> bRoot(n, -1);                  // for ODD: source-side vertex of edge that labeled this blossom
    vector<int> bPair(n, -1);                  // for EVEN: source-side vertex's partner that brought us here

    // Per-vertex (baseY, baseGblDualShift) representation: y_now(v) computed
    // on demand depending on current label and gblDualShift.
    vector<double> baseY(n, 0.0);
    vector<double> baseGblDualShift(n, 0.0);

    double gblDualShift = 0.0;

    // Bailout tracking (delta1-style direct counter, analog of V1's
    // allowedGblDualShift / sLast):
    //
    //   `allowedGblDualShift` bounds how far gblDualShift can advance before
    //   some EVEN vertex's y_now would go below 0.  Initialized to
    //   yDual[sFirst]; updated as new EVEN entities join the tree.
    //
    //   `evenLast` tracks which EVEN entity has the smallest projected
    //   y_now=0 time (i.e., contributes the binding bound).  When the
    //   bailout fires, evenLast is the vertex/blossom whose match-status
    //   gets flipped (free in the non-trivial case) to capture the
    //   weight improvement.
    //
    // When we'd pop a PQ entry with eventTime >= allowedGblDualShift, the
    // bailout has won; we break out and dispatch to the bailout branch.
    double allowedGblDualShift = matching.yDual[sFirst];
    int evenLast = sFirst;

    // Update bailout bound when a new EVEN entity joins the tree.  Computes
    // its projected y_now=0 time as baseGblDualShift[v] + baseY[v]; takes
    // the min with the current allowedGblDualShift; updates evenLast if a
    // tighter bound is found.
    //
    // Call this from GROW (when an EVEN-label vertex is freshly assigned)
    // and from BLOSSOM events (when a new blossom super-structure is
    // created with EVEN label, taking the min over its sub-blossoms'
    // bounds).
    auto updateAllowedShift = [&](int v) {
        double bound = baseGblDualShift[v] + baseY[v];
        if (bound < allowedGblDualShift) {
            allowedGblDualShift = bound;
            evenLast = v;
        }
    };

    // Initialize: sFirst becomes EVEN at gblDualShift=0.
    label[sFirst] = EVEN;
    baseY[sFirst] = matching.yDual[sFirst];
    baseGblDualShift[sFirst] = 0.0;

    // y_now(v) computed from current label, baseY, gblDualShift:
    //   EVEN  → yDual[v] = baseY - (gblDualShift - baseGblDualShift)
    //   ODD   → yDual[v] = baseY + (gblDualShift - baseGblDualShift)
    //   UNLABELED → yDual[v] = baseY (unchanged)
    auto yNow = [&](int v) -> double {
        int b = findBase(v);
        double Δ_v = gblDualShift - baseGblDualShift[v];
        if (label[b] == EVEN) return baseY[v] - Δ_v;
        if (label[b] == ODD)  return baseY[v] + Δ_v;
        return baseY[v];
    };

    // Edge slack at current gblDualShift:
    //   slack(u, v) = yNow(u) + yNow(v) + Σ z(B for B containing both u and v) - wght
    // For trivial blossoms (just vertices), no z contribution.  Non-trivial
    // blossoms contribute z(B) for each blossom containing both endpoints.
    // Computed on demand.
    auto edgeSlack = [&](int e) -> double {
        int u = graph.esrc[e], v = graph.etgt[e];
        double s = yNow(u) + yNow(v) - graph.wght[e];
        // TODO: add Σ z(B) for shared blossoms.  Skipped in this sketch.
        return s;
    };

    priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>> pq;

    // Schedule edge `e` in the PQ.  Compute its event time based on labels of
    // its endpoints' blossoms.  Stale entries (from earlier label states) get
    // filtered at pop time by re-checking edgeSlack.
    auto scanEdge = [&](int e) {
        int u = graph.esrc[e], v = graph.etgt[e];
        int bu = findBase(u), bv = findBase(v);
        if (bu == bv) return;                           // intra-blossom edge, ignore
        if (label[bu] == ODD || label[bv] == ODD) return; // ODD endpoint can't be source

        double slack = edgeSlack(e);
        double tightAt;
        if (label[bu] == EVEN && label[bv] == UNLABELED) {
            // GROW: outer to unlabeled.  Slack decreases by 1·gblDualShift per step.
            tightAt = gblDualShift + slack;
        } else if (label[bu] == UNLABELED && label[bv] == EVEN) {
            tightAt = gblDualShift + slack;
        } else if (label[bu] == EVEN && label[bv] == EVEN) {
            // BLOSSOM or AUGMENT: outer-outer.  Both endpoints' duals
            // decrease, so slack decreases by 2·gblDualShift per step.
            tightAt = gblDualShift + slack / 2.0;
        } else {
            return;  // both UNLABELED: no event
        }
        pq.push({tightAt, e, EDGE_EVENT});
    };

    // Schedule blossom `B` for expand event.  Fires when z(B) reaches 0.
    // ODD-blossom-only: EVEN blossoms have z growing with Δ (no upper-bound
    // event), so they never enter the PQ.  Only ODD blossoms have z
    // shrinking with Δ; expand must fire when z hits 0 to prevent
    // dual non-negativity violation.  Trivial blossoms (single
    // vertices) have no z, so skip them too.
    auto scanBlossom = [&](int B) {
        if (label[B] != ODD) return;
        if (B < graph.numVtxs) return;  // trivial blossom (single vertex), no z
        double tightAt = gblDualShift + matching.zDual[B] / 2.0;
        pq.push({tightAt, B, BLOSSOM_EVENT});
    };

    // Initial scan: outgoing edges of sFirst.
    for (int k = graph.off[sFirst]; k < graph.off[sFirst + 1]; ++k) {
        scanEdge(graph.adj[k]);
    }

    int augmentEndpoint = -1;  // free vertex reached, ending iteration successfully

    while (!pq.empty()) {
        PQEntry top = pq.top(); pq.pop();

        // Bailout check: if this event would advance gblDualShift past
        // allowedGblDualShift, the bailout has won — dispatch via the
        // post-loop branch with augmentEndpoint < 0 and current evenLast.
        if (top.eventTime >= allowedGblDualShift) break;

        gblDualShift = top.eventTime;

        if (top.type == EDGE_EVENT) {
            int e = top.id;
            // Stale-entry guard: predicted tightness assumed certain label
            // states; if labels have changed, the prediction is invalid.
            if (edgeSlack(e) > 1e-12) continue;  // not actually tight (real-valued tolerance)

            int u = graph.esrc[e], v = graph.etgt[e];
            int bu = findBase(u), bv = findBase(v);
            if (bu == bv) continue;

            if (label[bu] == EVEN && label[bv] == UNLABELED) {
                // GROW: extend the tree.
                // bv becomes ODD (T-blossom), bv's mate becomes EVEN.
                // TODO: full GROW logic, including pulling in the matched
                // partner as the next EVEN vertex, scanning its edges.
                // After labeling the new EVEN vertex w with baseY=yDual[w]
                // and baseGblDualShift=gblDualShift, call
                // updateAllowedShift(w) to incorporate w's bound into
                // the bailout tracking.
            } else if (label[bu] == UNLABELED && label[bv] == EVEN) {
                // GROW (symmetric).  See above; call updateAllowedShift
                // on the freshly-labeled EVEN vertex.
                // TODO.
            } else if (label[bu] == EVEN && label[bv] == EVEN) {
                // BLOSSOM or AUGMENT.  Walk back along parent pointers from
                // both endpoints; if they meet at an LCA in the tree, it's
                // a blossom; if they don't, it's an augmenting path to a
                // free vertex, but in single-source we always have one tree
                // so this is always a blossom unless one endpoint reaches
                // the free root.
                // TODO: LCA detection, blossom shrink (assigning a fresh
                // blossom ID, calling baseParent assignments via union-find,
                // setting z(B)=0, scheduling expand event, rescanning
                // edges incident to new blossom).  Also call
                // updateAllowedShift on the new EVEN super-blossom (its
                // bound is the min over sub-blossoms' bounds).
            }
        } else {
            // BLOSSOM_EVENT: expand the blossom.
            int B = top.id;
            if (label[B] != ODD) continue;  // stale
            // TODO: expand B.  Sub-blossoms get exposed, with appropriate
            // labels (alternating around the blossom cycle).  Rescan edges.
        }

        if (augmentEndpoint >= 0) break;
    }

    if (augmentEndpoint < 0) {
        // No augmenting path found.  Bailout fires.
        //
        // Set gblDualShift to allowedGblDualShift to drive evenLast's y_now
        // to 0.  Two cases:
        //
        //   Trivial bailout (evenLast == sFirst): the iteration's source
        //   vertex itself has the smallest dual; no other EVEN vertex
        //   provides a tighter bound.  No splice, no weight improvement;
        //   we just apply the dual update so that yDual[sFirst] = 0,
        //   marking sFirst as effectively "given up" for this iteration.
        //
        //   Non-trivial bailout (evenLast != sFirst): some EVEN vertex
        //   other than sFirst hits y_now=0 first.  Splice the matching
        //   along the partial alternating tree from evenLast back to
        //   sFirst, freeing evenLast and matching sFirst along the way.
        //   Cardinality unchanged; weight strictly increases by
        //   allowedGblDualShift.
        //
        // Then apply the dual update for all visited vertices/blossoms.
        gblDualShift = allowedGblDualShift;

        if (evenLast != sFirst) {
            // TODO: walk the partial alternating tree from evenLast back
            // to sFirst, lifting through nested blossoms as needed.
            // For each tree edge along the path, flip its M-status.
            // Free evenLast (matching.mate[evenLast] = -1, after picking
            // up its current mate to flip the connecting M-edge).
            // Mirror of V1's splice walk in lines ~150-165 of
            // hungarian_mwm_csr.cpp, with the additional complication of
            // walking through (and possibly expanding) nested blossoms.
        }

        // Apply the dual update: persist gblDualShift adjustments to baseY
        // values for all visited EVEN/ODD vertices and blossoms.
        // TODO: traverse visited set, set baseY[v] = yNow(v) and
        // baseGblDualShift[v] = 0; similar for blossom z values.

        return false;
    }

    // Augmenting path found.  Walk back from augmentEndpoint along parent
    // pointers, flipping matching status.
    // TODO: augmentation walk through nested blossoms (recursive lifting
    // when the path goes through a blossom's interior).

    // Apply the dual update: persist gblDualShift adjustments to baseY values
    // for all visited vertices/blossoms.
    // TODO: traverse visited set, set baseY[v] = yNow(v) and baseGblDualShift[v] = 0.

    return true;
}

Matching mwmGeneral(const Graph& graph) {
    Matching matching{
        vector<int>(graph.numVtxs, -1),
        vector<double>(graph.numVtxs, 0.0),
        vector<double>()  // zDual grows on demand
    };

    // Initialize y(v) so that for every edge (u,v), y(u) + y(v) >= wght(u,v).
    // Standard choice: y(v) = max incident edge weight / 2 (each edge gets
    // half from each endpoint).  Could also be y(v) = max incident weight,
    // matching V1's bipartite init.
    for (int v = 0; v < graph.numVtxs; ++v) {
        double maxW = 0.0;
        for (int k = graph.off[v]; k < graph.off[v + 1]; ++k) {
            int e = graph.adj[k];
            if (graph.wght[e] > maxW) maxW = graph.wght[e];
        }
        matching.yDual[v] = maxW / 2.0;
    }

    for (int v = 0; v < graph.numVtxs; ++v) {
        if (matching.mate[v] < 0) augment(graph, v, matching);
    }

    return matching;
}

// What's missing from this sketch (TODOs):
//
//   - GROW step: when an edge (s, v) becomes tight with s EVEN and v UNLABELED,
//     label v as ODD, label v's mate as EVEN, and rescan edges from the new
//     EVEN vertex.  About 30 lines.
//
//   - LCA detection in BLOSSOM step: walk parent pointers from both endpoints
//     of the tight edge, find common ancestor.  Use a marking scheme
//     (similar to G2's path1/path2 = strue trick).  About 25 lines.
//
//   - Blossom shrink: assign fresh blossom ID B, set baseParent[v] = B for
//     all sub-blossoms in the cycle, set z(B) = 0, label[B] = EVEN, rescan
//     edges incident to B (which means all edges from any sub-blossom
//     to outside).  About 40 lines.
//
//   - Blossom expand: when z(B) reaches 0 for an ODD blossom B, expose
//     sub-blossoms with alternating labels around the cycle.  Rescan their
//     edges.  About 50 lines.
//
//   - Augmentation walk through nested blossoms: when the path goes through
//     a blossom B, the in-edge and out-edge of B in the path identify two
//     vertices on B's cycle; the matching flip walks the shorter arc of the
//     cycle (the one of correct parity).  Recursive for nested blossoms.
//     About 40 lines.
//
//   - Bailout: track allowedGblDualShift across visited EVEN
//     blossoms (analog of V1's allowedGblDualShift across visited S).  When
//     no GROW or AUGMENT can fire, the free root either has dual driven to
//     0 (trivial bailout) or some other EVEN vertex does and we flip along
//     the partial tree (non-trivial bailout).  Detection mechanism and
//     dispatch branch implemented; the splice walk through nested
//     blossoms remains TODO.  About 50 lines.
//
//   - Final dual update: after augmenting or bailing, persist
//     yNow values back into baseY for the next iteration's initial
//     state.  About 20 lines.
//
//   - Blossom z(B) terms in edgeSlack: when both endpoints are inside a
//     non-trivial blossom B (or share enclosing blossoms), the slack must
//     include Σ z(B).  About 10 lines, careful indexing.
//
// Total expected size: ~500-600 lines once filled in.  This sketch shows the
// scaffolding (data structures, scan_edge, scan_blossom, main event loop,
// PQ structure, virtual contraction setup).  The combinatorial mechanics
// (shrink, expand, augmentation walk, bailout) are stubbed as TODOs.

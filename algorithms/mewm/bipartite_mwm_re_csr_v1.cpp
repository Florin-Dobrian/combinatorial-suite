// Hungarian-style maximum-weight bipartite matching.
// Sparse, flat CSR, single-source SAP-with-potentials.
// Lazy binary-heap PQ (Strategy A from the Dijkstra doc).

#include <vector>
#include <queue>
#include <limits>

using namespace std;

// Partition-local vertex IDs: s in [0, sNumVtxs), t in [0, tNumVtxs).
// Flat CSR for the S side. Edges go S -> T.
//   off[s..s+1]  indexes into adj[] and wght[].
//   adj[k] = T-vertex index (in [0, tNumVtxs))
//   wght[k] = w(s, adj[k])

struct Graph {
    int sNumVtxs, tNumVtxs;
    vector<int> off;   // size sNumVtxs + 1
    vector<int> adj;       // size m, T-vertex indices
    vector<double> wght;    // size m
};

// Output: sMate[s] = matched t, or -1 if s unmatched.
//         tMate[t] = matched s, or -1 if t unmatched.
//         sDual, tDual: LP-duality witness potentials.
//   Reduced cost of edge (s, t) is  sDual[s] + tDual[t] - w(s, t),  >= 0 invariant.
//   Matching edges are tight: sDual[s] + tDual[t] == w(s, t).

struct Matching {
    vector<int> sMate, tMate;
    vector<double> sDual, tDual;
};

static constexpr double INF = numeric_limits<double>::infinity();

// One Hungarian iteration: find a shortest augmenting path from `sFirst` (free in S)
// to some free T-vertex, augment along it, update duals selectively.
// Returns true on augmentation, false on the no-improvement bailout
// (it is cheaper to leave sFirst unmatched than to augment).

bool augment(const Graph& graph, int sFirst, Matching& matching) {
    // Per-iteration scratch.  In production these would be reused across
    // iterations with selective cleanup; here, allocated fresh for clarity.
    vector<double> sDist(graph.sNumVtxs, INF), tDist(graph.tNumVtxs, INF);
    vector<int> tPrev(graph.tNumVtxs, -1);    // for t in tree, the s before it
    vector<int> sPrev(graph.sNumVtxs, -1);    // for s in tree, the t before it (-1 if s == sFirst)
    vector<char> sVstd(graph.sNumVtxs, 0), tVstd(graph.tNumVtxs, 0);
    vector<int> sPrcdQue, tPrcdQue;   // for selective dual update

    sDist[sFirst] = 0.0;
    sPrcdQue.push_back(sFirst);  sVstd[sFirst] = 1;

    // Best-so-far cost-to-reach of any S-vertex in this Dijkstra.
    // The nonperfect-MWBM bailout (no-improvement termination): if no
    // free T is cheaper than allowedGblDualShift, it is better to leave sFirst unmatched
    // and free sLast instead.
    double allowedGblDualShift = matching.sDual[sFirst];
    int sLast = sFirst;

    // Lazy binary-heap PQ keyed by tentative distance to a T-vertex.
    using P = pair<double, int>;
    priority_queue<P, vector<P>, greater<P>> tPrcbQue;

    // Edge relaxation from a known S-vertex `s`.  For each outgoing edge
    // (s, t), compute the reduced-cost slack and the candidate distance
    // tCrtDist = sCrtDist + slack.  If tCrtDist improves the best known
    // distance to t, record it, set tPrev[t] = s, mark t processed (so
    // it appears once in tPrcdQue for the dual update), and push a fresh
    // entry onto the PQ.  Stale entries are filtered at pop time, so we
    // never call decrease-key — this is the lazy-PQ idiom from the
    // Dijkstra doc (Strategy A).
    //
    // Implemented as a lambda capturing the enclosing scratch arrays
    // by reference (`[&]`) so the call site reads cleanly as
    // `relaxFromS(s)` rather than passing eight parameters every time.

    auto relaxFromS = [&](int s) {
        double sCrtDist = sDist[s];
        for (int k = graph.off[s]; k < graph.off[s + 1]; ++k) {
            int t = graph.adj[k];
            double stEdgWght = graph.wght[k];
            double slack = (matching.sDual[s] + matching.tDual[t]) - stEdgWght;   // >= 0
            double tCrtDist = sCrtDist + slack;
            if (tCrtDist < tDist[t]) {
                tDist[t] = tCrtDist;
                tPrev[t] = s;
                if (!tVstd[t]) { tPrcdQue.push_back(t); tVstd[t] = 1; }
                tPrcbQue.push({tCrtDist, t});                                  // insert, don't update
            }
        }
    };

    relaxFromS(sFirst);

    int tLast = -1;
    double gblDualShift;

    while (!tPrcbQue.empty()) {
        auto [tCrtDist, t] = tPrcbQue.top(); tPrcbQue.pop();

        // Stale-entry guard: a later push gave t a shorter distance.
        if (tCrtDist > tDist[t]) continue;

        // Bailout: the cheapest reachable T is no better than
        // unmatching some already-processed S-vertex.  Two one-sided
        // constraints on the iteration's global dual shift gblDualShift
        // clash: closing the slack on this T's incoming edge needs
        // gblDualShift at least tCrtDist; keeping all visited sDual
        // values non-negative bounds gblDualShift at allowedGblDualShift
        // from above.  When tCrtDist >= allowedGblDualShift no legal
        // gblDualShift exists, since any global dual shift large enough
        // to make this T's incoming edge tight would drive at least one
        // visited S-vertex's dual negative.  Formally tDual[t] would
        // need to be adjusted by gblDualShift - tCrtDist which must be
        // >= 0, thus gblDualShift >= tCrtDist.  Therefore if
        // tCrtDist >= allowedGblDualShift we would either make the
        // dual of at least one S-vertex negative or we would not make
        // any difference in terms of weight.
        if (tCrtDist >= allowedGblDualShift) break;

        if (matching.tMate[t] < 0) {
            // Free T reached.  This is the SAP endpoint.
            tLast = t;
            gblDualShift = tCrtDist;
            break;
        }

        // T is matched.  Follow the unique matching edge to its S-mate
        // and relax outgoing edges from there.
        int sNext = matching.tMate[t];
        sDist[sNext] = tCrtDist;
        sPrev[sNext] = t;
        if (!sVstd[sNext]) { sPrcdQue.push_back(sNext); sVstd[sNext] = 1; }

        if (tCrtDist + matching.sDual[sNext] < allowedGblDualShift) {
            allowedGblDualShift = tCrtDist + matching.sDual[sNext];
            sLast = sNext;
        }

        relaxFromS(sNext);
    }

    if (tLast < 0) {
        // No free T was reached.  Either the comparison test
        // tCrtDist >= allowedGblDualShift fired mid-loop (forward progress structurally
        // blocked), or the PQ exhausted before any free T appeared.
        // Both cases: no cardinality-augmenting path; possibly a
        // weight-augmenting path if sLast != sFirst.
        gblDualShift = allowedGblDualShift;
        if (sLast != sFirst) {
            // Walk the weight-augmenting path: predecessor chain from sLast
            // back to sFirst, flipping alternation along the way.
            // sLast (currently matched) becomes free; sFirst (currently free)
            // becomes matched along the partial tree.  Cardinality unchanged,
            // weight strictly increases.
            matching.sMate[sLast] = -1;
            int s = sLast;
            while (s != sFirst) {
                int t = sPrev[s];                // T preceding s in the tree
                int ss = tPrev[t];           // S preceding that t
                matching.sMate[ss] = t;
                matching.tMate[t] = ss;
                s = ss;
            }
        }
    } else {
        // Augment along the predecessor chain from tLast back to sFirst.
        int t = tLast;
        while (t >= 0) {
            int s = tPrev[t];
            int tt = sPrev[s];         // T before s in the tree, -1 if s == sFirst
            matching.sMate[s] = t;
            matching.tMate[t] = s;
            t = tt;
        }
    }

    // Selective dual update: shift duals of processed vertices by
    // (gblDualShift - dist[v]).  Vertices not in sPrcdQue / tPrcdQue are untouched.
    // Reduced cost remains non-negative for all edges by construction.
    for (int s : sPrcdQue) {
        double lclDualShift = gblDualShift - sDist[s];
        if (lclDualShift > 0) matching.sDual[s] -= lclDualShift;
    }
    for (int t : tPrcdQue) {
        double lclDualShift = gblDualShift - tDist[t];
        if (lclDualShift > 0) matching.tDual[t] += lclDualShift;
    }

    return tLast >= 0;
}

// Outer loop: iterate over free S-vertices, augmenting from each.
// Initialization: sDual[s] = max incident edge weight; tDual[t] = 0.
// Reduced costs are non-negative at start.

Matching mwm_hungarian(const Graph& graph) {
    Matching matching{
        vector<int>(graph.sNumVtxs, -1),
        vector<int>(graph.tNumVtxs, -1),
        vector<double>(graph.sNumVtxs, 0.0),
        vector<double>(graph.tNumVtxs, 0.0)
    };

    // Init sDual[s] = max w(s, *), so every s has at least one tight edge.
    for (int s = 0; s < graph.sNumVtxs; ++s) {
        double maxWght = 0.0;
        for (int k = graph.off[s]; k < graph.off[s + 1]; ++k) {
            if (graph.wght[k] > maxWght) maxWght = graph.wght[k];
        }
        matching.sDual[s] = maxWght;
    }

    for (int s = 0; s < graph.sNumVtxs; ++s) {
        if (matching.sMate[s] < 0) augment(graph, s, matching);
    }

    return matching;
}

// Floating-point handling, LEDA-style: pre-scale doubles to integers-in-doubles
// so all algorithm arithmetic is exact.
//
// The Hungarian algorithm performs additions and subtractions on weights and
// duals (slack = sDual + tDual - wght, dual updates by gblDualShift - dist,
// path arithmetic accumulating slacks).  With arbitrary doubles, roundoff
// errors break the invariants the algorithm relies on:
//
//   - Slack should be exactly 0 on M-edges (tightness).  Roundoff makes
//     slack drift to small positive/negative values, breaking complementary
//     slackness checks and causing the bailout test (tCrtDist >= allowed)
//     to fire spuriously or fail to fire when it should.
//   - Dual non-negativity (sDual >= 0, tDual >= 0).  A dual that should land
//     at exactly 0 may end up at -1e-15, causing the per-vertex
//     non-negativity check to skip dual updates that should fire.
//   - Distance comparisons in the PQ.  Tied distances may sort
//     non-deterministically under roundoff, producing different traversal
//     orders on different runs of the same input.
//
// The fix: pre-scale the weights to a representation where all algorithm
// arithmetic stays in the exact-integer regime of double precision.
//
//   w1[e] = sign(w[e]) * floor(|w[e]| * S) / S
//
// where S is the largest power of two with S < 2^53 / (3 * C),
// C = max |w[e]|.  This guarantees that:
//
//   1. Each w1[e] is an integer-valued double (the floor produces an exact
//      integer; division by S=2^k preserves it as an integer-valued double).
//   2. All sums of up to 3*C*S of these integers stay within 2^53, the
//      exact-integer range of double precision.
//   3. The 3 factor accounts for the worst-case path arithmetic
//      (sDual + tDual + slack-style accumulations along an alternating path).
//
// Optimal matching weight differs from the unscaled version by at most
// n * 3 * C * 2^-52, which is the rounding error from the floor() truncation
// summed over n M-edges.  For typical n=10^6, C=10^4, this gives error
// bound roughly 10^-9 — far below any practical sensitivity threshold.
//
// Usage: compute w1 once at startup, pass it to mwm_hungarian as the
// weight array.  The output matching is exact for the rounded weights,
// approximate (within the bound above) for the originals.  Free vertices
// have dual exactly 0; M-edges are exactly tight.  All algorithm invariants
// hold bit-exactly.
//
// Reference: LEDA's MWBM_SCALE_WEIGHTS and MWBM_check_weight_bound.

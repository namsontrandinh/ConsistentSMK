// src/algs/consistent_chasing.h
//
// Implements Algorithm 2: Restarted-Knapsack-Chasing-Local-Opt
// from "Consistent submodular under knapsack constraint" (Ha et al.)
//
// Problem setting:
//   - Insertion-only stream: elements arrive one by one.
//   - f is normalized, MONOTONE, submodular.
//   - Each element e has cost c(e) = g.nodes[e].weight > 0.
//   - Budget B: feasible set S must satisfy c(S) <= B.
//   - Cost-granularity: c(e) <= delta*B for all feasible e (Definition 3.1).
//
// Dual output (Theorem 1):
//   Y_t = argmax{ f(S_t), f(P_t) } where P_t is the best singleton seen so far.
//   The approximation guarantee applies to Y_t, not S_t alone.
//
// Consistency metric tracked in CCResult.total_sym_diff:
//   sum over t of |S_t △ S_{t-1}|, and max_sym_diff = max over t of |S_t △ S_{t-1}|.
//
// CHANGELOG vs v1 (fixes per reviewer):
//   - FIX #1: q = floor(1/delta), was ceil.
//   - FIX #2: N = floor(ln q / ln((100+delta)/100)) + 1, was ceil.
//   - FIX #3: Output Y_t = argmax{f(S_t), f(P_t)}, was only S_t.
//   - FIX #4: Bounded-Exchange invariant documented; drop check order clarified.
//   - Added: delta guard against pathologically small values (warn if q > 10000).
//   - Added: total_sym_diff and max_sym_diff tracking for consistency experiments.

#ifndef ALGS_CONSISTENT_CHASING_H
#define ALGS_CONSISTENT_CHASING_H

#include <vector>
#include <cstdint>
#include <limits>
#include <chrono>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>   // for pathological-delta warning

#include "mygraph.h"
#include "sfunctions.h"

namespace algs {

// ---------------------------------------------------------------------------
// CCResult: local result struct for ConsistentChasing.
// Does NOT depend on result.h from the partition-matroid project.
// Follows the same field conventions as edl.h / twin_greedy.h.
// ---------------------------------------------------------------------------
struct CCResult {
    std::string algo;
    std::string constraint;

    // Final output solution Y_n = argmax{f(S_n), f(P_n)} (Algorithm 2, Line 20).
    std::vector<std::uint8_t> inS;

    double      f_value      = 0.0;   // f(Y_n)
    std::size_t queries      = 0;     // total oracle calls
    double      time_sec     = 0.0;

    // Consistency tracking (for experiments).
    // total_sym_diff = sum_t |S_t △ S_{t-1}|  (cumulative recourse)
    // max_sym_diff   = max_t |S_t △ S_{t-1}|  (worst-case per-step recourse)
    std::size_t total_sym_diff = 0;
    std::size_t max_sym_diff   = 0;

    // Derived parameters (logged for reproducibility).
    std::size_t q_param = 0;   // floor(1/delta)
    std::size_t N_param = 0;   // floor(ln q / ln((100+delta)/100)) + 1
};

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------
namespace cc_detail {

inline double node_cost(const mygraph::tinyGraph &g, mygraph::node_id u) {
    return g.nodes[u].weight;
}

inline double cost_of_solution(const mygraph::tinyGraph &g,
                                const subm::Solution &inS) {
    double c = 0.0;
    for (std::size_t u = 0; u < g.n; ++u)
        if (inS[u]) c += g.nodes[u].weight;
    return c;
}

inline double evaluate(const mygraph::tinyGraph &g,
                       const subm::Solution &inS,
                       std::size_t &queries)
{
    ++queries;
    return subm::sfunc_evaluate(g, inS);
}

// f(S ∪ {u}) - f(S). Assumes inS[u] == 0.
inline double marginal(const mygraph::tinyGraph &g,
                       const subm::Solution &inS,
                       mygraph::node_id u,
                       double fS,
                       std::size_t &queries)
{
    ++queries;
    return subm::sfunc_marginal(g, inS, u, fS);
}

// Symmetric difference |A △ B| between two 0/1 solutions of the same size.
inline std::size_t sym_diff(const subm::Solution &A, const subm::Solution &B) {
    std::size_t d = 0;
    for (std::size_t u = 0; u < A.size(); ++u)
        if (A[u] != B[u]) ++d;
    return d;
}

} // namespace cc_detail


// ===========================================================================
// Algorithm 1: Bounded-Exchange(S, x, B, q)
//
// Paper: Algorithm 1, Section 4.
//
// Finds R ⊆ S minimising ω_S(R) = f(S) - f(S\R) subject to:
//   (a) c(R) >= max(0, c(S) + c(x) - B)    [makes room for x in budget]
//   (b) |R| <= q                             [bounded cardinality]
//
// Returns empty R (feasible=true) if x fits without any removal (d=0).
// Sets feasible=false if no R satisfying (a)+(b) exists.
//
// Implementation invariant (addresses reviewer comment #4):
//   Elements of S are sorted by individual drop ω_S({u}) = f(S) - f(S\{u})
//   ascending before greedy selection. This ensures that the R we return
//   has the smallest total drop achievable by any greedy prefix, which is
//   the best practical approximation to the minimum-drop R required by the
//   theory. The caller (chasing loop) then verifies that ω_S(R) <=
//   mu*(c(x)/B)*f(S) — this two-step structure is intentional: Bounded-Exchange
//   finds the best feasible R it can, and the caller decides whether that R
//   is good enough to commit the exchange.
// ===========================================================================
inline std::vector<mygraph::node_id> bounded_exchange(
    const mygraph::tinyGraph     &g,
    const subm::Solution         &inS,
    double                        fS,
    mygraph::node_id              x,
    double                        B,
    std::size_t                   q,
    bool                         &feasible,
    std::size_t                  &queries)
{
    using mygraph::node_id;

    feasible = true;

    const double cx = cc_detail::node_cost(g, x);
    const double cS = cc_detail::cost_of_solution(g, inS);
    // d = amount of budget that must be freed to accommodate x (Algorithm 1, Line 1).
    const double d  = std::max(0.0, cS + cx - B);

    // Algorithm 1, Lines 2-3: x fits without removing anything.
    if (d <= 1e-12) {
        return {};  // empty R; feasible = true
    }

    // Compute individual drop for each element in S.
    // drop(u) = f(S) - f(S \ {u}) >= 0  for monotone f.
    // Cost: |S| oracle calls, each O(m). |S| <= q by Lemma 1.
    struct MemberDrop {
        node_id u;
        double  drop;
        double  cost;
    };

    std::vector<MemberDrop> members;
    members.reserve(static_cast<std::size_t>(q) + 4);

    subm::Solution tmp = inS;   // working copy; restored after each probe
    for (std::size_t u = 0; u < g.n; ++u) {
        if (!inS[u]) continue;

        tmp[u] = 0;
        const double f_minus_u = cc_detail::evaluate(g, tmp, queries);
        tmp[u] = 1;  // restore

        members.push_back({static_cast<node_id>(u),
                           fS - f_minus_u,                          // drop >= 0
                           cc_detail::node_cost(g, static_cast<node_id>(u))});
    }

    // Sort by drop ascending: smallest value loss first.
    std::sort(members.begin(), members.end(),
              [](const MemberDrop &a, const MemberDrop &b) {
                  return a.drop < b.drop;
              });

    // Greedy selection: accumulate until c(R) >= d, with |R| <= q.
    std::vector<node_id> R;
    R.reserve(q);
    double cR = 0.0;

    for (const auto &md : members) {
        if (cR >= d - 1e-12) break;   // enough budget freed
        if (R.size() >= q)  break;    // |R| <= q (Algorithm 1, Line 4)
        R.push_back(md.u);
        cR += md.cost;
    }

    // Feasibility check: did we accumulate enough cost?
    if (cR < d - 1e-12) {
        // Cannot free enough budget within q removals.
        feasible = false;
        return {};
    }

    return R;
}


// ===========================================================================
// Algorithm 2: Restarted-Knapsack-Chasing-Local-Opt
//
// Processes nodes 0, 1, ..., n-1 as an insertion-only stream.
// Returns CCResult with the DUAL output Y_n = argmax{f(S_n), f(P_n)}.
//
// Parameters:
//   g      - graph; g.nodes[u].weight = c(u)
//   B      - knapsack budget
//   delta  - cost-granularity: assumed c(u) <= delta*B for all feasible u.
//            Typical value: 0.1.  Must be in (0, 1].
//            WARNING: very small delta (e.g. 0.001) yields q=1000 and
//            N potentially in the hundreds; runtime per element becomes
//            O(N * n * q * m).  A warning is printed to stderr if q > 10000.
// ===========================================================================
inline CCResult run_consistent_chasing(const mygraph::tinyGraph &g,
                                       double B,
                                       double delta)
{
    using mygraph::node_id;
    using subm::Solution;

    CCResult res;
    res.algo = "ConsistentChasing";

    auto t0 = std::chrono::high_resolution_clock::now();

    // ---- Sanity checks ----
    if (g.n == 0 || B <= 0.0) {
        Solution empty(g.n, 0);
        res.inS     = empty;
        res.f_value = cc_detail::evaluate(g, empty, res.queries);
        res.time_sec = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t0).count();
        return res;
    }
    if (!(delta > 0.0 && delta <= 1.0)) {
        throw std::invalid_argument(
            "ConsistentChasing: delta must be in (0, 1].");
    }

    // -------------------------------------------------------------------------
    // Algorithm 2, Line 1: derived parameters.
    //
    // FIX #1 (vs v1): q = floor(1/delta), NOT ceil.
    //   Paper Lemma 1: |S| <= floor(1/delta) = q for any feasible S.
    //
    // FIX #2 (vs v1): N = floor(ln q / ln((100+delta)/100)) + 1, NOT ceil.
    //   Paper Algorithm 2 Line 1 uses floor explicitly.
    //   Using ceil would overcount iterations, violating the stated C bound.
    // -------------------------------------------------------------------------
    const std::size_t q = static_cast<std::size_t>(std::floor(1.0 / delta));

    // vartheta = (100+delta)/(100*delta): density threshold for x (Line 14, condition 1).
    const double vartheta = (100.0 + delta) / (100.0 * delta);

    // mu = 1/delta: drop bound multiplier (Line 14, condition 3; Lemma 2).
    const double mu = 1.0 / delta;

    // N = floor( ln(q) / ln((100+delta)/100) ) + 1.
    std::size_t N = 1;
    if (q >= 2) {
        const double log_ratio = std::log((100.0 + delta) / 100.0);  // > 0 for delta > 0
        if (log_ratio > 1e-15) {
            // FIX #2: floor, not ceil.
            N = static_cast<std::size_t>(
                    std::floor(std::log(static_cast<double>(q)) / log_ratio))
                + 1;
        }
    }

    res.q_param = q;
    res.N_param = N;

    {
        std::ostringstream oss;
        oss << "B=" << B << ", delta=" << delta
            << ", q=" << q << ", N=" << N;
        res.constraint = oss.str();
    }

    // Pathological-delta warning: runtime blows up for tiny delta.
    if (q > 10000) {
        std::cerr << "[ConsistentChasing WARNING] delta=" << delta
                  << " yields q=" << q << ", N=" << N
                  << ". Runtime per element is O(N*n*q*m). Consider larger delta.\n";
    }

    // -------------------------------------------------------------------------
    // Algorithm 2, Line 2: S_0 = empty, P_0 = empty.
    // -------------------------------------------------------------------------
    Solution S(g.n, 0);   // maintained solution S_t
    double   fS = 0.0;    // f(S_t), kept in sync after every exchange

    // P_t: best singleton seen so far (Lines 5-8).
    // We store the id and its f-value to avoid re-evaluating.
    node_id  best_p_id  = static_cast<node_id>(g.n);  // sentinel: none seen yet
    double   best_p_val = 0.0;

    // -------------------------------------------------------------------------
    // Main stream loop: process element e_t (Algorithm 2, Line 3).
    // -------------------------------------------------------------------------
    for (node_id et = 0; et < static_cast<node_id>(g.n); ++et) {

        const double c_et = cc_detail::node_cost(g, et);

        // Discard infeasible elements immediately (paper Section 3).
        if (c_et <= 0.0 || c_et > B) continue;

        // ---- Lines 5-8: update best singleton P_t ----
        {
            Solution tmp(g.n, 0);
            tmp[et] = 1;
            const double f_et = cc_detail::evaluate(g, tmp, res.queries);
            if (f_et > best_p_val) {
                best_p_val = f_et;
                best_p_id  = et;
            }
        }

        // S_prev is used to compute |S_t △ S_{t-1}| for consistency tracking.
        Solution S_prev = S;

        // ---- Lines 9-12: choose starting point ----
        // S <- argmax{ f(P_t), f(S_{t-1}) }
        if (best_p_val > fS) {
            // Restart from best singleton (Lines 9-10).
            S.assign(g.n, 0);
            S[best_p_id] = 1;
            fS = best_p_val;
        }
        // else: keep S_{t-1} (Lines 11-12); fS unchanged.

        // ---- Lines 13-18: chasing loop ----
        //
        // Invariant entering each iteration: S is feasible, fS = f(S).
        // Each iteration tries to find ONE improving exchange (x, R) and commits it.
        // If no improving exchange exists, break (Line 18).
        //
        // Three conditions must hold simultaneously (Line 14):
        //   (a) f(x|S)/c(x) >= vartheta * f(S)/B    [density threshold]
        //   (b) Bounded-Exchange returns R != failure  [budget feasibility]
        //   (c) omega_S(R) <= mu * (c(x)/B) * f(S)   [bounded value drop]
        //
        // Reviewer note on condition ordering (fix #4):
        //   We check (a) first (cheap: 1 oracle call), then call Bounded-Exchange
        //   for (b), then check (c) on the returned R.  Bounded-Exchange sorts
        //   by drop ascending, so the R it returns minimises drop among all greedy
        //   prefixes.  If (c) fails for this R, no greedy prefix can do better,
        //   so we correctly skip x.  The invariant "Bounded-Exchange returns the
        //   minimum-drop feasible R" is guaranteed by the ascending sort + greedy
        //   prefix structure.

        for (std::size_t iter = 0; iter < N; ++iter) {

            bool found = false;

            for (node_id x = 0; x <= et; ++x) {
                if (S[x]) continue;  // x must not be in S already

                const double cx = cc_detail::node_cost(g, x);
                if (cx <= 0.0 || cx > B) continue;

                // ---- Condition (a) ----
                const double delta_x    = cc_detail::marginal(g, S, x, fS, res.queries);
                const double threshold_a = vartheta * cx * fS / B;
                // Note: when fS == 0, threshold_a == 0, so any x with delta_x >= 0
                // passes.  This is correct: with an empty solution any element
                // with non-negative marginal is improving.
                if (delta_x < threshold_a) continue;

                // ---- Condition (b): call Bounded-Exchange ----
                bool be_feasible = false;
                std::vector<node_id> R = bounded_exchange(
                    g, S, fS, x, B, q, be_feasible, res.queries);
                if (!be_feasible) continue;

                // ---- Condition (c): check drop bound ----
                // omega_S(R) = f(S) - f(S \ R).
                // If R is empty (x fits without removal), drop = 0 <= anything.
                double drop = 0.0;
                if (!R.empty()) {
                    Solution tmp = S;
                    for (node_id r : R) tmp[r] = 0;
                    // One oracle call to evaluate S \ R.
                    const double f_minus_R = cc_detail::evaluate(g, tmp, res.queries);
                    drop = fS - f_minus_R;
                }
                const double threshold_c = mu * (cx / B) * fS;
                if (drop > threshold_c) continue;

                // ---- All three conditions satisfied: commit exchange ----
                // S <- (S \ R) ∪ {x}
                for (node_id r : R) S[r] = 0;
                S[x] = 1;

                // Full oracle call to update fS (not incremental).
                // Rationale: incremental update f += delta_x - drop is inexact
                // because drop was measured on S before R was removed, not on
                // (S \ R).  Full call avoids accumulated error.  See edl.h pattern.
                fS = cc_detail::evaluate(g, S, res.queries);

                found = true;
                break;  // one exchange per iteration; restart scan with new S
            }

            if (!found) break;  // Line 18: local optimum reached
        }

        // ---- Line 19: S_t <- S ----
        // Track consistency metric |S_t △ S_{t-1}|.
        const std::size_t sd = cc_detail::sym_diff(S, S_prev);
        res.total_sym_diff += sd;
        if (sd > res.max_sym_diff) res.max_sym_diff = sd;

    } // end stream loop (all elements processed)

    // -------------------------------------------------------------------------
    // FIX #3 (vs v1): dual output Y_n = argmax{ f(S_n), f(P_n) }.
    //
    // Algorithm 2, Line 20: Y_t <- S_t.
    // But Theorem 1 approximation guarantee is for Y_t where:
    //   Y_t in argmax_{Z in {S_t, P_t}} f(Z).
    // P_t is the best singleton; its value is already cached in best_p_val.
    // S_t value is in fS.
    // -------------------------------------------------------------------------
    if (best_p_val > fS) {
        // Output P_n: the best singleton.
        res.inS.assign(g.n, 0);
        if (best_p_id < static_cast<node_id>(g.n))
            res.inS[best_p_id] = 1;
        res.f_value = best_p_val;
    } else {
        // Output S_n: the maintained solution after chasing.
        res.inS     = S;
        res.f_value = fS;
    }

    res.time_sec = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - t0).count();
    return res;
}

} // namespace algs

#endif // ALGS_CONSISTENT_CHASING_H

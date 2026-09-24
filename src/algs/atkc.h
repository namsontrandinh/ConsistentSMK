// src/algs/atkc.h
//
// Implements:
//   Algorithm 3: Threshold-Size-Prioritized-Exchange(S, x, tau, B, q)
//   Algorithm 4: Adaptive-Threshold-Knapsack-Chasing (ATKC)
// from "Consistent Monotone Submodular Maximization under Budget
// Constraints" (Ha et al.), Section 5.
//
// Problem setting (identical to RKCL, Section 4):
//   - Insertion-only stream: elements arrive one by one, V_t = V_{t-1} u {e_t}.
//   - f is normalized (f(empty)=0), MONOTONE, submodular.
//   - Each element e has cost c(e) = subm_obj_common::node_cost(g, e) > 0.
//   - Budget B: feasible set S must satisfy c(S) <= B.
//   - Cost-granularity: c(e) >= delta*B for all feasible e -> |S| <= q = floor(1/delta).
//
// ATKC reuses two facts already established for RKCL (NOT re-derived here):
//   - q = floor(1/delta) bounds the size of any feasible solution.
//   - M_t <= f(O_t) <= q * M_t, where M_t = f(P_t) is the best singleton value.
//
// Output convention (paper, Line 30):
//   Y_t = argmax over { P_t } u { S_{j,t}, G_{j,t} : 0 <= j <= J } of f(.)
//
// Consistency metric tracked (same convention as rkcl_im.h / csmk.h):
//   |Y_t triangle Y_{t-1}|  (symmetric difference of consecutive OUTPUTS).
//
// -----------------------------------------------------------------------
// KNOWN IMPLEMENTATION GAP (documented, same spirit as Bounded-Exchange):
//
// Algorithm 3, Line 3 requires, for each removal size k, the BEST subset
// R' of S with |R'| = k (i.e. argmax over C(|S|, k) subsets). Enumerating
// all subsets is combinatorially infeasible once |S| grows past a handful
// of elements. This implementation instead uses the same greedy heuristic
// already used for Bounded-Exchange in this project: elements of S are
// sorted once by ascending individual drop
//     drop(u) = f(S) - f(S \ {u})
// and, for a given k, R_k is taken to be the length-k PREFIX of that
// sorted order. This does not guarantee the theoretical R is found, but
// it is the standard practical relaxation for this class of exchange
// subroutines (see consistent_chasing analysis notes). If R_k does not
// meet the budget/threshold test, k is increased by one (grow the
// prefix), following the same size-priority order the paper mandates.
// -----------------------------------------------------------------------
//
// COMPUTATIONAL COST WARNING:
// Per element arrival, ATKC runs (J+1) independent lanes, each performing
// up to N successful-exchange iterations, each iteration scanning up to
// |V_t| candidates, each candidate call to Algorithm 3 costing up to
// O(|S|) oracle calls just to rank members by drop, plus up to q more
// oracle calls to test growing prefixes. With the IC objective (Monte
// Carlo spread oracle), a single oracle call is itself expensive. Test
// with a SMALL prefix (--active_n or a truncated graph) and a reduced
// KIC_MC before running on the full email.bin / fb.bin streams. Use
// ATKCParams::max_x_scan_per_iter and ATKCParams::chase_cap (both 0 =
// unlimited / theoretical) to bound worst-case runtime during debugging.

#ifndef ALGS_ATKC_H
#define ALGS_ATKC_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "mygraph.h"
#include "sfunctions_impl.h"
#include "objectvalue/objectvalue_common.h"
#include "algs/result.h"

namespace algs {

using ATKCClock = std::chrono::high_resolution_clock;
using mygraph::node_id;
using subm::Solution;

// =====================================================
// Parameters (paper inputs: delta, epsilon).
//
// max_x_scan_per_iter == 0:
//   Scan all candidates x in V_t \ S during one chase iteration
//   (theoretical behaviour).
// max_x_scan_per_iter > 0:
//   Stop scanning after this many candidates have been TRIED (not
//   necessarily accepted) in a single chase iteration. Useful safety
//   valve for early experiments; when capped, the certified guarantee
//   no longer strictly applies (same caveat style as RKCLParams.chase_cap
//   in rkcl_im.h).
//
// chase_cap == 0:
//   Use the theoretical N = floor(3q/delta) + 1 for every lane.
// chase_cap > 0:
//   Use min(N, chase_cap) successful-exchange iterations per lane.
// =====================================================
struct ATKCParams {
    double delta = 0.1;
    double eps   = 0.5;

    std::uint64_t max_x_scan_per_iter = 0;

    // Practical deviation for Algorithm 4 Lines 27-29 (one-element
    // augmentation), which the manuscript specifies as an unbounded
    // scan of Vt \ Sj,t at EVERY lane, EVERY timestep -> O(active_n)
    // per lane per timestep, O(active_n^2 * lanes) overall.
    // 0 = scan every candidate (theoretical).
    // > 0 = scan only the this many most-recently-arrived candidates
    // each augmentation pass. Not certified by the theorem when capped.
    std::uint64_t aug_scan_cap = 0;
    std::uint64_t chase_cap           = 0;
};

// =====================================================
// Derived parameters (Line 1): q, J, N.
// =====================================================
struct ATKCDerivedParams {
    double delta = 0.0;
    double eps   = 0.0;

    std::uint64_t q = 0;   // floor(1/delta)
    std::uint64_t J = 0;   // ceil( log_{1+eps}(q/eps) )   -> lanes 0..J
    std::uint64_t theoretical_N = 0;  // floor(3q/delta) + 1
    std::uint64_t used_N        = 0;  // min(theoretical_N, chase_cap) if capped

    bool chase_was_capped = false;
};

inline ATKCDerivedParams atkc_derive_params(const ATKCParams& p) {
    if (!(p.delta > 0.0 && p.delta <= 1.0)) {
        throw std::invalid_argument("ATKC: delta must be in (0,1].");
    }
    if (!(p.eps > 0.0 && p.eps < 1.0)) {
        throw std::invalid_argument("ATKC: eps must be in (0,1).");
    }

    ATKCDerivedParams out;
    out.delta = p.delta;
    out.eps   = p.eps;

    // q = floor(1/delta)  (Line 1, shared with RKCL).
    out.q = static_cast<std::uint64_t>(std::floor(1.0 / p.delta + 1e-12));
    out.q = std::max<std::uint64_t>(1, out.q);

    // J = ceil( log_{1+eps}(q/eps) ) = ceil( ln(q/eps) / ln(1+eps) ).
    const double q_over_eps = static_cast<double>(out.q) / p.eps;
    if (q_over_eps <= 1.0) {
        out.J = 0;
    } else {
        const double ln_ratio = std::log(q_over_eps);
        const double ln_base  = std::log(1.0 + p.eps);
        out.J = static_cast<std::uint64_t>(std::ceil(ln_ratio / ln_base - 1e-9));
    }

    // N = floor(3q/delta) + 1  (Line 1).
    out.theoretical_N = static_cast<std::uint64_t>(
        std::floor(3.0 * static_cast<double>(out.q) / p.delta + 1e-9)) + 1;

    out.used_N = out.theoretical_N;
    // chase_cap is applied by the caller once it knows theoretical_N,
    // stored here for convenience/logging.
    return out;
}

// =====================================================
// Basic helpers (same conventions as edl.h / multistream.h / rkcl_im.h).
// =====================================================
inline double atkc_node_cost(const mygraph::tinyGraph& g, node_id u) {
    return subm_obj_common::node_cost(g, u);
}

inline double atkc_cost_of(const mygraph::tinyGraph& g, const Solution& inS) {
    subm::solution_check_size(g, inS);
    double c = 0.0;
    for (std::size_t u = 0; u < g.n; ++u) {
        if (inS[u]) c += atkc_node_cost(g, static_cast<node_id>(u));
    }
    return c;
}

inline std::uint64_t atkc_symmetric_difference(const Solution& A, const Solution& B) {
    if (A.size() != B.size()) {
        throw std::invalid_argument("ATKC: solution sizes differ.");
    }
    std::uint64_t d = 0;
    for (std::size_t u = 0; u < A.size(); ++u) {
        if (static_cast<bool>(A[u]) != static_cast<bool>(B[u])) ++d;
    }
    return d;
}

inline double atkc_eval_singleton(
    const mygraph::tinyGraph& g, Solution& tmp, node_id e, std::uint64_t& queries)
{
    tmp[e] = 1;
    const double val = subm::sfunc_evaluate(g, tmp);
    ++queries;
    tmp[e] = 0;
    return val;
}

// =====================================================
// Algorithm 3: Threshold-Size-Prioritized-Exchange(S, x, tau, B, q)
//
// Returns feasible=true and R (possibly empty) on success, or
// feasible=false (the paper's "bottom" symbol) on failure.
//
// See file-level comment for the documented greedy relaxation of Line 3.
// =====================================================
struct Algo3Result {
    bool feasible = false;
    std::vector<node_id> R;
};

inline Algo3Result threshold_size_prioritized_exchange(
    const mygraph::tinyGraph& g,
    const Solution& S,
    double fS,
    double cS,
    node_id x,
    double tau,
    double B,
    std::uint64_t q,
    std::uint64_t& queries)
{
    Algo3Result out;

    const double cx = atkc_node_cost(g, x);
    if (!(cx > 0.0) || cx > B) {
        return out; // x itself can never be feasible under this budget.
    }

    // ---- k = 0: try R = empty (Line 1, k=0 case) ----
    if (cS + cx <= B + 1e-9) {
        const double gain = subm::sfunc_marginal(g, S, x, fS);
        ++queries;
        if (gain >= tau * cx - 1e-12) {
            out.feasible = true;
            out.R.clear();
            return out;
        }
    }

    // ---- Rank members of S by ascending individual drop ----
    // drop(u) = f(S) - f(S \ {u}) >= 0 for monotone f.
    // This ordering is the same size-priority-compatible heuristic used
    // for Bounded-Exchange: it approximates "smallest total value loss
    // for a given removal size" without enumerating all C(|S|,k) subsets.
    struct MemberDrop {
        node_id u;
        double  drop;
        double  cost;
    };

    std::vector<MemberDrop> members;
    members.reserve(S.size());

    Solution tmp = S;
    for (std::size_t u = 0; u < g.n; ++u) {
        if (!S[u]) continue;

        tmp[u] = 0;
        const double f_minus_u = subm::sfunc_evaluate(g, tmp);
        ++queries;
        tmp[u] = 1;

        members.push_back({static_cast<node_id>(u),
                            fS - f_minus_u,
                            atkc_node_cost(g, static_cast<node_id>(u))});
    }

    std::sort(members.begin(), members.end(),
              [](const MemberDrop& a, const MemberDrop& b) { return a.drop < b.drop; });

    const std::size_t kmax = std::min(members.size(), static_cast<std::size_t>(q));

    // ---- k = 1 .. kmax: grow the prefix R_k one element at a time ----
    std::vector<node_id> Rk;
    Rk.reserve(kmax);
    double cRk = 0.0;

    for (std::size_t k = 1; k <= kmax; ++k) {
        Rk.push_back(members[k - 1].u);
        cRk += members[k - 1].cost;

        const double cost_after = cS - cRk + cx;
        if (cost_after > B + 1e-9) {
            continue; // not yet feasible; keep growing the prefix.
        }

        Solution afterS = S;
        for (node_id r : Rk) afterS[r] = 0;
        afterS[x] = 1;

        const double f_after = subm::sfunc_evaluate(g, afterS);
        ++queries;

        const double gain = f_after - fS;
        if (gain >= tau * cx - 1e-12) {
            out.feasible = true;
            out.R = Rk;
            return out;
        }
    }

    return out; // Line 8: bottom (failure).
}

// =====================================================
// Per-run diagnostics: which candidate produced Y_t.
// Encoding: -2 = P_t, 2*j = S_{j,t}, 2*j+1 = G_{j,t}.
// =====================================================
inline int atkc_encode_source_P() { return -2; }
inline int atkc_encode_source_S(std::uint64_t j) { return static_cast<int>(2 * j); }
inline int atkc_encode_source_G(std::uint64_t j) { return static_cast<int>(2 * j + 1); }

// =====================================================
// Consistency statistics from the direct (single-pass) run.
// Same field conventions as RKCLConsistencyStats in rkcl_im.h.
// =====================================================
struct ATKCConsistencyStats {
    Result final_result;

    ATKCDerivedParams params;

    std::vector<std::uint64_t> per_step_symmetric_difference;
    std::vector<int>           per_step_best_source; // see atkc_encode_source_*

    std::uint64_t C = 0;
    std::uint64_t cumulative_consistency = 0;
    std::uint64_t initial_changes = 0;

    std::uint64_t total_lane_restarts     = 0;
    std::uint64_t total_successful_exchanges = 0;

    std::uint64_t total_queries = 0;
    double total_time_sec = 0.0;
};

// =====================================================
// Algorithm 4: Adaptive-Threshold-Knapsack-Chasing.
//
// Single-pass, insertion-only run over V = {0, ..., active_n-1}.
// Unlike EDL/MultiStream (which are re-run from scratch on every
// prefix for consistency measurement), ATKC directly maintains
// S_{j,t-1} -> S_{j,t} for every lane, exactly as specified by the
// paper's pseudocode (this mirrors how rkcl_im.h maintains RKCL's
// S_{t-1} -> S_t; do not confuse this with the EDL/MultiStream
// prefix-rerun pattern).
// =====================================================
inline ATKCConsistencyStats run_ATKC_consistency_prefix(
    const mygraph::tinyGraph& g,
    double B,
    const ATKCParams& user_params,
    std::size_t active_n)
{
    ATKCConsistencyStats stats;
    const auto t0 = ATKCClock::now();

    active_n = std::min(active_n, g.n);

    stats.final_result.algo = "ATKC";
    {
        std::ostringstream oss;
        oss << "B=" << B << ", delta=" << user_params.delta
            << ", eps=" << user_params.eps;
        stats.final_result.constraint = oss.str();
    }
    stats.final_result.inS.assign(g.n, 0);

    if (active_n == 0 || B <= 0.0) {
        stats.final_result.f_value = subm::sfunc_evaluate(g, stats.final_result.inS);
        ++stats.total_queries;
        stats.final_result.queries = stats.total_queries;
        stats.total_time_sec = std::chrono::duration<double>(
            ATKCClock::now() - t0).count();
        stats.final_result.time_sec = stats.total_time_sec;
        return stats;
    }

    stats.params = atkc_derive_params(user_params);
    const std::uint64_t q = stats.params.q;
    const std::uint64_t J = stats.params.J;

    std::uint64_t N = stats.params.theoretical_N;
    if (user_params.chase_cap > 0 && user_params.chase_cap < N) {
        N = user_params.chase_cap;
        stats.params.chase_was_capped = true;
    }
    stats.params.used_N = N;

    stats.per_step_symmetric_difference.assign(active_n, 0);
    stats.per_step_best_source.assign(active_n, atkc_encode_source_P());

    // ---- Line 2: initial state ----
    // P_0 = empty; S_{j,0} = empty for all j = 0..J.
    Solution Pt(g.n, 0);
    double fP = 0.0; // f(empty) = 0 for normalized f; we still evaluate once below.
    {
        Solution empty(g.n, 0);
        fP = subm::sfunc_evaluate(g, empty);
        ++stats.total_queries;
    }
    bool has_singleton = false;
    Solution singleton_tmp(g.n, 0);

    const std::size_t num_lanes = static_cast<std::size_t>(J) + 1;

    std::vector<Solution> Sj(num_lanes, Solution(g.n, 0));
    std::vector<double>   fSj(num_lanes, fP);   // f(empty) initially
    std::vector<double>   cSj(num_lanes, 0.0);

    std::vector<Solution> Gj(num_lanes, Solution(g.n, 0));
    std::vector<double>   fGj(num_lanes, fP);
    std::vector<double>   cGj(num_lanes, 0.0);

    Solution Yt_prev(g.n, 0);
    double fYt_prev = fP;

    // ---- Line 3: main stream loop ----
    for (std::size_t t = 1; t <= active_n; ++t) {
        const node_id et = static_cast<node_id>(t - 1);

        // ---- Lines 5-8: update best feasible singleton P_t ----
        const double c_et = atkc_node_cost(g, et);
        if (c_et > 0.0 && c_et <= B) {
            const double val = atkc_eval_singleton(g, singleton_tmp, et, stats.total_queries);
            if (val > fP) {
                Pt.assign(g.n, 0);
                Pt[et] = 1;
                fP = val;
                has_singleton = true;
            }
        }

        // ---- Line 9: M_t = f(P_t) ----
        const double Mt = fP;

        Solution Yt;
        double fYt = 0.0;
        int best_source = atkc_encode_source_P();

        // ---- Lines 10-13: degenerate case M_t = 0 ----
        if (!(Mt > 1e-12)) {
            for (std::size_t j = 0; j < num_lanes; ++j) {
                Sj[j].assign(g.n, 0);
                fSj[j] = 0.0;
                cSj[j] = 0.0;
                Gj[j] = Sj[j];
                fGj[j] = 0.0;
                cGj[j] = 0.0;
            }
            Yt = Pt;       // Pt is empty here since Mt == 0.
            fYt = fP;      // == 0
            best_source = atkc_encode_source_P();
        }
        else {
            // ---- Lines 14-29: one lane per threshold tau_{j,t} ----
            for (std::size_t j = 0; j < num_lanes; ++j) {
                const double tau_j = (Mt / (3.0 * B)) *
                                      std::pow(1.0 + user_params.eps, static_cast<double>(j));

                Solution S;
                double fS = 0.0;
                double cS = 0.0;

                // ---- Lines 16-19: per-lane restart ----
                if (fSj[j] >= tau_j * cSj[j] - 1e-12) {
                    S = Sj[j];
                    fS = fSj[j];
                    cS = cSj[j];
                } else {
                    S.assign(g.n, 0);
                    fS = 0.0;
                    cS = 0.0;
                    ++stats.total_lane_restarts;
                }

                // ---- Lines 20-25: chase loop, up to N successful exchanges ----
                for (std::uint64_t iter = 0; iter < N; ++iter) {
                    bool found = false;

                    std::uint64_t scanned = 0;
                    for (node_id x = 0; x <= et; ++x) {
                        if (S[x]) continue;

                        if (user_params.max_x_scan_per_iter > 0 &&
                            scanned >= user_params.max_x_scan_per_iter)
                        {
                            break;
                        }
                        ++scanned;

                        Algo3Result ar = threshold_size_prioritized_exchange(
                            g, S, fS, cS, x, tau_j, B, q, stats.total_queries);

                        if (!ar.feasible) continue;

                        for (node_id r : ar.R) {
                            S[r] = 0;
                            cS -= atkc_node_cost(g, r);
                        }
                        S[x] = 1;
                        cS += atkc_node_cost(g, x);

                        // Full oracle call for the same reason documented
                        // in csmk.h: incremental accumulation across a
                        // removal-then-add step is inexact.
                        fS = subm::sfunc_evaluate(g, S);
                        ++stats.total_queries;

                        found = true;
                        ++stats.total_successful_exchanges;
                        break; // one exchange per iteration; rescan with new S.
                    }

                    if (!found) break; // Line 25: local optimum for this lane.
                }

                // ---- Line 26 ----
                Sj[j] = S;
                fSj[j] = fS;
                cSj[j] = cS;

                Gj[j] = S;
                fGj[j] = fS;
                cGj[j] = cS;

                // ---- Lines 27-29: one-element augmentation ----
                {
                    const bool aug_capped = (user_params.aug_scan_cap > 0);
                    const std::size_t total_candidates =
                        static_cast<std::size_t>(et) + 1;
                    const std::size_t aug_limit =
                        aug_capped
                            ? std::min<std::size_t>(user_params.aug_scan_cap, total_candidates)
                            : total_candidates;

                    for (std::size_t s = 0; s < aug_limit; ++s) {
                        const node_id x =
                            aug_capped
                                ? static_cast<node_id>(total_candidates - 1 - s)
                                : static_cast<node_id>(s);

                        if (Sj[j][x]) continue;

                        const double cx = atkc_node_cost(g, x);
                        if (!(cx > 0.0) || cSj[j] + cx > B + 1e-9) continue;

                        const double gain = subm::sfunc_marginal(g, Sj[j], x, fSj[j]);
                        ++stats.total_queries;

                        const double candidate_value = fSj[j] + gain;
                        if (candidate_value > fGj[j] + 1e-12) {
                            Gj[j] = Sj[j];
                            Gj[j][x] = 1;
                            fGj[j] = candidate_value;
                            cGj[j] = cSj[j] + cx;
                        }
                    }
                }
            }

            // ---- Line 30: Y_t = argmax over {P_t} u {S_j,t, G_j,t} ----
            Yt = Pt;
            fYt = fP;
            best_source = atkc_encode_source_P();

            for (std::size_t j = 0; j < num_lanes; ++j) {
                if (fSj[j] > fYt + 1e-12) {
                    Yt = Sj[j];
                    fYt = fSj[j];
                    best_source = atkc_encode_source_S(j);
                }
                if (fGj[j] > fYt + 1e-12) {
                    Yt = Gj[j];
                    fYt = fGj[j];
                    best_source = atkc_encode_source_G(j);
                }
            }
        }

        // ---- Consistency tracking on the OUTPUT Y_t (not on S_j,t) ----
        const std::uint64_t change = atkc_symmetric_difference(Yt_prev, Yt);

        if (t == 1) {
            stats.initial_changes = change;
            stats.per_step_symmetric_difference[0] = 0;
        } else {
            stats.per_step_symmetric_difference[t - 1] = change;
            stats.C = std::max(stats.C, change);
            stats.cumulative_consistency += change;
        }
        stats.per_step_best_source[t - 1] = best_source;

        Yt_prev = Yt;
        fYt_prev = fYt;
        (void)has_singleton; // kept for readability/debugging symmetry with rkcl_im.h
    }

    stats.final_result.inS = Yt_prev;
    stats.final_result.f_value = fYt_prev;
    stats.final_result.queries = stats.total_queries;

    stats.total_time_sec = std::chrono::duration<double>(ATKCClock::now() - t0).count();
    stats.final_result.time_sec = stats.total_time_sec;

    return stats;
}

inline ATKCConsistencyStats run_ATKC_consistency(
    const mygraph::tinyGraph& g, double B, const ATKCParams& params)
{
    return run_ATKC_consistency_prefix(g, B, params, g.n);
}

// Simple Result wrapper for experiments that do not need the detailed
// consistency trace (same convenience overload style as rkcl_im.h).
inline Result run_ATKC_prefix(
    const mygraph::tinyGraph& g, double B, const ATKCParams& params, std::size_t active_n)
{
    return run_ATKC_consistency_prefix(g, B, params, active_n).final_result;
}

inline Result run_ATKC(const mygraph::tinyGraph& g, double B, const ATKCParams& params) {
    return run_ATKC_prefix(g, B, params, g.n);
}

} // namespace algs

#endif // ALGS_ATKC_H
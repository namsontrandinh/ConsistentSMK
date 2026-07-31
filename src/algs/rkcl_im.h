// src/algs/rkcl_im.h
#ifndef ALGS_RKCL_IM_H
#define ALGS_RKCL_IM_H

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

using RKCLClock = std::chrono::high_resolution_clock;
using mygraph::node_id;
using subm::Solution;

// =====================================================
// Restarted-Knapsack-Chasing-Local-Opt (RKCL)
//
// Experimental implementation for monotone submodular
// maximization under a knapsack constraint.
//
// When the executable is compiled with the IC objective,
// sfunc_evaluate / sfunc_marginal become the influence
// spread oracle, so the same code runs IM on Facebook,
// Email, or any other graph loaded by tinyGraph.
//
// IMPORTANT:
// Consistency is measured by symmetric difference:
//     |S_t Delta S_{t-1}|
// exactly as in the ConsistentSMK manuscript.
// =====================================================

// -----------------------------------------------------
// Parameters.
//
// delta <= 0:
//   Automatically use the smallest observed ratio
//   min_{0<c(e)<=B} c(e)/B. This makes the empirical
//   granularity assumption valid for the loaded graph.
//
// theta:
//   lambda = 1/delta + theta.
//   The manuscript's explicit example uses theta = 0.01.
//
// chase_cap == 0:
//   Use the theoretical N.
// chase_cap > 0:
//   Use min(N, chase_cap). This is useful for large IM
//   experiments but the forced-local-optimum theorem is
//   then not certified when chase_cap < N.
// -----------------------------------------------------
struct RKCLParams {
    double delta = 0.0;
    double theta = 0.01;
    std::uint64_t chase_cap = 0;
    bool strict_granularity = true;
};

// -----------------------------------------------------
// Derived parameters.
// -----------------------------------------------------
struct RKCLDerivedParams {
    double delta = 0.0;
    double mu = 0.0;
    double lambda = 0.0;
    double gamma = 0.0;

    std::uint64_t q = 0;
    std::uint64_t theoretical_N = 0;
    std::uint64_t used_N = 0;

    std::uint64_t theoretical_C = 0;

    bool chase_was_capped = false;
};

// -----------------------------------------------------
// Statistics from the actual insertion-only run.
// Unlike EDL/MultiStream prefix reruns, RKCL directly
// maintains S_t from S_{t-1}.
// -----------------------------------------------------
struct RKCLConsistencyStats {
    Result final_result;

    RKCLDerivedParams params;

    std::vector<std::uint64_t>
        per_step_symmetric_difference;

    std::vector<std::uint64_t>
        per_step_exchanges;

    std::uint64_t C = 0;
    std::uint64_t cumulative_consistency = 0;
    std::uint64_t initial_changes = 0;

    std::uint64_t total_restarts = 0;
    std::uint64_t total_successful_exchanges = 0;

    std::uint64_t total_queries = 0;
    double total_time_sec = 0.0;
};

// =====================================================
// Basic helpers.
// =====================================================
inline double rkcl_node_cost(
    const mygraph::tinyGraph& g,
    node_id u)
{
    return subm_obj_common::node_cost(g, u);
}

inline double rkcl_cost_of(
    const mygraph::tinyGraph& g,
    const Solution& inS)
{
    subm::solution_check_size(g, inS);

    double c = 0.0;

    for (std::size_t u = 0;
         u < g.n;
         ++u)
    {
        if (inS[u]) {
            c += rkcl_node_cost(
                g,
                static_cast<node_id>(u));
        }
    }

    return c;
}

inline std::uint64_t
rkcl_symmetric_difference(
    const Solution& A,
    const Solution& B)
{
    if (A.size() != B.size()) {
        throw std::invalid_argument(
            "RKCL: solution sizes differ.");
    }

    std::uint64_t d = 0;

    for (std::size_t u = 0;
         u < A.size();
         ++u)
    {
        if (static_cast<bool>(A[u]) !=
            static_cast<bool>(B[u]))
        {
            ++d;
        }
    }

    return d;
}

inline double rkcl_eval_singleton(
    const mygraph::tinyGraph& g,
    Solution& tmp,
    node_id e,
    std::uint64_t& queries)
{
    tmp[e] = 1;

    const double val =
        subm::sfunc_evaluate(
            g,
            tmp);

    ++queries;

    tmp[e] = 0;

    return val;
}

// =====================================================
// Effective granularity of the loaded graph.
//
// delta_eff = min c(e)/B over feasible positive-cost
// elements. If every such element satisfies
// c(e) >= delta_eff * B, then q=floor(1/delta_eff)
// is a valid cardinality upper bound.
// =====================================================
inline double rkcl_effective_delta(
    const mygraph::tinyGraph& g,
    double B,
    std::size_t active_n)
{
    if (!(B > 0.0)) {
        return 0.0;
    }

    active_n =
        std::min(
            active_n,
            g.n);

    double delta_eff =
        std::numeric_limits<double>::infinity();

    for (std::size_t idx = 0;
         idx < active_n;
         ++idx)
    {
        const node_id e =
            static_cast<node_id>(idx);

        const double ce =
            rkcl_node_cost(
                g,
                e);

        if (ce > 0.0 &&
            ce <= B)
        {
            delta_eff =
                std::min(
                    delta_eff,
                    ce / B);
        }
    }

    if (!std::isfinite(delta_eff)) {
        return 0.0;
    }

    return std::min(
        1.0,
        delta_eff);
}

inline RKCLDerivedParams
rkcl_derive_params(
    const mygraph::tinyGraph& g,
    double B,
    const RKCLParams& p,
    std::size_t active_n)
{
    if (!(B > 0.0)) {
        throw std::invalid_argument(
            "RKCL: B must be positive.");
    }

    if (!(p.theta > 0.0)) {
        throw std::invalid_argument(
            "RKCL: theta must be positive.");
    }

    const double delta_eff =
        rkcl_effective_delta(
            g,
            B,
            active_n);

    if (!(delta_eff > 0.0)) {
        throw std::invalid_argument(
            "RKCL: no feasible positive-cost element.");
    }

    double delta =
        (p.delta > 0.0)
            ? p.delta
            : delta_eff;

    if (!(delta > 0.0 &&
          delta <= 1.0))
    {
        throw std::invalid_argument(
            "RKCL: delta must be in (0,1].");
    }

    // A user-provided delta larger than delta_eff violates
    // c(e) >= delta B for at least one feasible element.
    if (delta >
        delta_eff +
            1e-12)
    {
        if (p.strict_granularity) {
            throw std::invalid_argument(
                "RKCL: supplied delta violates graph cost granularity. "
                "Use --delta 0 for automatic delta or choose a smaller value.");
        }

        delta =
            delta_eff;
    }

    RKCLDerivedParams out;

    out.delta =
        delta;

    out.mu =
        1.0 /
        delta;

    out.lambda =
        out.mu +
        p.theta;

    // gamma = 1 + delta(lambda-mu)
    //       = 1 + delta*theta.
    out.gamma =
        1.0 +
        delta *
        p.theta;

    out.q =
        static_cast<std::uint64_t>(
            std::floor(
                1.0 /
                delta +
                1e-12));

    out.q =
        std::max<std::uint64_t>(
            1,
            out.q);

    if (!(out.gamma > 1.0)) {
        throw std::runtime_error(
            "RKCL: gamma must be greater than 1.");
    }

    const long double ln_q =
        std::log(
            static_cast<long double>(
                out.q));

    const long double ln_gamma =
        std::log(
            static_cast<long double>(
                out.gamma));

    long double raw_N =
        1.0L;

    if (out.q > 1) {
        raw_N =
            std::floor(
                ln_q /
                ln_gamma) +
            1.0L;
    }

    if (raw_N >
        static_cast<long double>(
            std::numeric_limits<std::uint64_t>::max()))
    {
        throw std::overflow_error(
            "RKCL: theoretical N overflows uint64_t.");
    }

    out.theoretical_N =
        static_cast<std::uint64_t>(
            raw_N);

    out.theoretical_N =
        std::max<std::uint64_t>(
            1,
            out.theoretical_N);

    out.used_N =
        out.theoretical_N;

    if (p.chase_cap > 0 &&
        p.chase_cap <
            out.used_N)
    {
        out.used_N =
            p.chase_cap;

        out.chase_was_capped =
            true;
    }

    // The manuscript proves C <= N(q+1)
    // for the uncapped theoretical run.
    if (out.theoretical_N >
        std::numeric_limits<std::uint64_t>::max() /
            (out.q + 1))
    {
        out.theoretical_C =
            std::numeric_limits<std::uint64_t>::max();
    }
    else {
        out.theoretical_C =
            out.theoretical_N *
            (out.q + 1);
    }

    return out;
}

// =====================================================
// Bounded-Exchange output.
// =====================================================
struct RKCLExchange {
    bool feasible = false;

    Solution after_removal;

    std::vector<node_id> removed;

    double released_cost = 0.0;
    double value_after_removal = 0.0;
    double loss = 0.0;
};

// =====================================================
// Practical Bounded-Exchange.
//
// The manuscript defines R as a minimum-loss deletion set.
// Computing that argmin exactly is exponential in |S|.
//
// For scalable Facebook/Email IM experiments we use:
//   1) iterative minimum loss-per-released-cost deletion;
//   2) verify the required bounded-loss inequality later;
//   3) if the greedy deletion violates the bound, use R=S.
//
// The fallback R=S is exactly the universal construction
// from Lemma 2 when mu=1/delta and c(x)>=delta B.
// Thus it preserves feasibility and the universal bounded
// exchange guarantee under the granularity assumption.
// =====================================================
inline RKCLExchange
rkcl_bounded_exchange(
    const mygraph::tinyGraph& g,
    const Solution& S,
    double fS,
    double cS,
    node_id x,
    double B,
    std::uint64_t& queries)
{
    RKCLExchange out;

    const double cx =
        rkcl_node_cost(
            g,
            x);

    if (!(cx > 0.0) ||
        cx > B)
    {
        return out;
    }

    const double deficit =
        std::max(
            0.0,
            cS +
                cx -
                B);

    out.after_removal =
        S;

    out.value_after_removal =
        fS;

    if (deficit <= 1e-12)
    {
        out.feasible = true;
        out.loss = 0.0;
        return out;
    }

    // Greedily remove the element with minimum
    // current loss / released cost.
    double released = 0.0;
    double cur_value = fS;

    while (released +
               1e-12 <
           deficit)
    {
        bool found = false;

        node_id best_u = 0;

        double best_score =
            std::numeric_limits<double>::infinity();

        double best_value_after =
            cur_value;

        for (std::size_t idx = 0;
             idx < g.n;
             ++idx)
        {
            if (!out.after_removal[idx]) {
                continue;
            }

            const node_id u =
                static_cast<node_id>(idx);

            const double cu =
                rkcl_node_cost(
                    g,
                    u);

            if (!(cu > 0.0)) {
                continue;
            }

            out.after_removal[u] = 0;

            const double val_minus =
                subm::sfunc_evaluate(
                    g,
                    out.after_removal);

            ++queries;

            out.after_removal[u] = 1;

            const double loss_u =
                std::max(
                    0.0,
                    cur_value -
                        val_minus);

            const double score =
                loss_u /
                cu;

            if (!found ||
                score <
                    best_score)
            {
                found = true;
                best_u = u;
                best_score = score;
                best_value_after =
                    val_minus;
            }
        }

        if (!found) {
            break;
        }

        out.after_removal[best_u] = 0;

        out.removed.push_back(
            best_u);

        released +=
            rkcl_node_cost(
                g,
                best_u);

        cur_value =
            best_value_after;
    }

    if (released +
            1e-12 <
        deficit)
    {
        // Should not happen for a feasible S and cx<=B,
        // but keep the routine robust.
        return RKCLExchange{};
    }

    out.feasible = true;
    out.released_cost = released;

    out.value_after_removal =
        cur_value;

    out.loss =
        std::max(
            0.0,
            fS -
                cur_value);

    return out;
}

// Universal fallback R=S from Lemma 2.
inline RKCLExchange
rkcl_full_eviction(
    const mygraph::tinyGraph& g,
    const Solution& S,
    double fS)
{
    RKCLExchange out;

    out.feasible = true;

    out.after_removal.assign(
        g.n,
        0);

    out.value_after_removal =
        0.0;

    out.loss =
        fS;

    for (std::size_t idx = 0;
         idx < g.n;
         ++idx)
    {
        if (!S[idx]) {
            continue;
        }

        const node_id u =
            static_cast<node_id>(idx);

        out.removed.push_back(
            u);

        out.released_cost +=
            rkcl_node_cost(
                g,
                u);
    }

    return out;
}

// =====================================================
// Find one (lambda,mu)-improving element.
//
// The theoretical algorithm allows any improving x.
// For deterministic and scalable experiments we scan
// V_t in arrival order and accept the first valid x.
// =====================================================
struct RKCLImprovingMove {
    bool found = false;

    node_id x = 0;

    RKCLExchange exchange;

    Solution next_S;

    double next_cost = 0.0;
    double next_value = 0.0;
};

inline RKCLImprovingMove
rkcl_find_improving_move(
    const mygraph::tinyGraph& g,
    const Solution& S,
    double fS,
    double cS,
    double B,
    const RKCLDerivedParams& params,
    std::size_t active_n,
    std::uint64_t& queries)
{
    RKCLImprovingMove move;

    active_n =
        std::min(
            active_n,
            g.n);

    for (std::size_t idx = 0;
         idx < active_n;
         ++idx)
    {
        if (S[idx]) {
            continue;
        }

        const node_id x =
            static_cast<node_id>(idx);

        const double cx =
            rkcl_node_cost(
                g,
                x);

        if (!(cx > 0.0) ||
            cx > B)
        {
            continue;
        }

        const double marginal =
            subm::sfunc_marginal(
                g,
                S,
                x,
                fS);

        ++queries;

        // Density condition:
        // f(x|S)/c(x) >= lambda f(S)/B.
        if (marginal /
                cx +
                1e-12 <
            params.lambda *
                fS /
                B)
        {
            continue;
        }

        RKCLExchange ex =
            rkcl_bounded_exchange(
                g,
                S,
                fS,
                cS,
                x,
                B,
                queries);

        if (!ex.feasible) {
            continue;
        }

        const double allowed_loss =
            params.mu *
            (cx / B) *
            fS;

        // If the practical low-loss deletion misses the
        // theoretical bound, use the universal R=S fallback.
        if (ex.loss >
            allowed_loss +
                1e-10)
        {
            ex =
                rkcl_full_eviction(
                    g,
                    S,
                    fS);

            // Under valid granularity and mu=1/delta,
            // this inequality must hold.
            if (ex.loss >
                allowed_loss +
                    1e-10)
            {
                continue;
            }
        }

        Solution next =
            ex.after_removal;

        next[x] = 1;

        const double next_cost =
            cS -
            ex.released_cost +
            cx;

        if (next_cost >
            B +
                1e-9)
        {
            continue;
        }

        // Evaluate the exchanged solution exactly through
        // the configured submodular oracle.
        const double next_value =
            subm::sfunc_evaluate(
                g,
                next);

        ++queries;

        // Protect against Monte-Carlo noise / numerical noise.
        if (next_value <=
            fS +
                1e-12)
        {
            continue;
        }

        move.found = true;
        move.x = x;
        move.exchange =
            std::move(ex);

        move.next_S =
            std::move(next);

        move.next_cost =
            next_cost;

        move.next_value =
            next_value;

        return move;
    }

    return move;
}

// =====================================================
// Core dynamic insertion-only RKCL run.
// =====================================================
inline RKCLConsistencyStats
run_RKCL_consistency_prefix(
    const mygraph::tinyGraph& g,
    double B,
    const RKCLParams& user_params,
    std::size_t active_n)
{
    RKCLConsistencyStats stats;

    const auto t0 =
        RKCLClock::now();

    active_n =
        std::min(
            active_n,
            g.n);

    stats.final_result.algo =
        "RKCL-IM";

    {
        std::ostringstream oss;

        oss
            << "B="
            << B
            << ", delta="
            << user_params.delta
            << ", theta="
            << user_params.theta
            << ", chase_cap="
            << user_params.chase_cap;

        stats.final_result.constraint =
            oss.str();
    }

    stats.final_result.inS.assign(
        g.n,
        0);

    if (active_n == 0 ||
        B <= 0.0)
    {
        stats.final_result.f_value =
            subm::sfunc_evaluate(
                g,
                stats.final_result.inS);

        ++stats.total_queries;

        stats.final_result.queries =
            stats.total_queries;

        stats.total_time_sec =
            std::chrono::duration<double>(
                RKCLClock::now() -
                t0)
                .count();

        stats.final_result.time_sec =
            stats.total_time_sec;

        return stats;
    }

    stats.params =
        rkcl_derive_params(
            g,
            B,
            user_params,
            active_n);

    stats.per_step_symmetric_difference.assign(
        active_n,
        0);

    stats.per_step_exchanges.assign(
        active_n,
        0);

    // S_0 = empty.
    Solution S(
        g.n,
        0);

    double fS =
        subm::sfunc_evaluate(
            g,
            S);

    ++stats.total_queries;

    double cS = 0.0;

    // P_0 = empty.
    Solution P(
        g.n,
        0);

    double fP = fS;
    double cP = 0.0;

    bool has_singleton = false;

    Solution singleton_tmp(
        g.n,
        0);

    for (std::size_t t = 1;
         t <= active_n;
         ++t)
    {
        const node_id et =
            static_cast<node_id>(
                t - 1);

        const Solution S_prev =
            S;

        // ---------------------------------------------
        // Update best feasible singleton P_t.
        // ---------------------------------------------
        const double cet =
            rkcl_node_cost(
                g,
                et);

        if (cet > 0.0 &&
            cet <= B)
        {
            const double singleton_value =
                rkcl_eval_singleton(
                    g,
                    singleton_tmp,
                    et,
                    stats.total_queries);

            if (!has_singleton ||
                singleton_value >
                    fP)
            {
                has_singleton = true;

                P.assign(
                    g.n,
                    0);

                P[et] = 1;

                fP =
                    singleton_value;

                cP =
                    cet;
            }
        }

        // ---------------------------------------------
        // Restart from better of S_{t-1} and P_t.
        // ---------------------------------------------
        if (has_singleton &&
            fP >
                fS +
                    1e-12)
        {
            S =
                P;

            fS =
                fP;

            cS =
                cP;

            ++stats.total_restarts;
        }

        // ---------------------------------------------
        // Chase exchange local optimum.
        // The theorem guarantees fewer than N successful
        // exchanges when using theoretical N.
        // ---------------------------------------------
        std::uint64_t step_exchanges = 0;

        for (std::uint64_t i = 0;
             i < stats.params.used_N;
             ++i)
        {
            RKCLImprovingMove move =
                rkcl_find_improving_move(
                    g,
                    S,
                    fS,
                    cS,
                    B,
                    stats.params,
                    t,
                    stats.total_queries);

            if (!move.found) {
                break;
            }

            S =
                std::move(
                    move.next_S);

            fS =
                move.next_value;

            cS =
                move.next_cost;

            ++step_exchanges;
            ++stats.total_successful_exchanges;
        }

        stats.per_step_exchanges[t - 1] =
            step_exchanges;

        // ---------------------------------------------
        // Actual consistency:
        // |S_t Delta S_{t-1}|.
        // ---------------------------------------------
        const std::uint64_t change =
            rkcl_symmetric_difference(
                S_prev,
                S);

        if (t == 1)
        {
            stats.initial_changes =
                change;

            // Definition normally starts at t >= 2.
            stats.per_step_symmetric_difference[0] =
                0;
        }
        else {
            stats.per_step_symmetric_difference[t - 1] =
                change;

            stats.C =
                std::max(
                    stats.C,
                    change);

            stats.cumulative_consistency +=
                change;
        }
    }

    stats.final_result.inS =
        S;

    stats.final_result.f_value =
        fS;

    stats.final_result.queries =
        stats.total_queries;

    stats.total_time_sec =
        std::chrono::duration<double>(
            RKCLClock::now() -
            t0)
            .count();

    stats.final_result.time_sec =
        stats.total_time_sec;

    return stats;
}

inline RKCLConsistencyStats
run_RKCL_consistency(
    const mygraph::tinyGraph& g,
    double B,
    const RKCLParams& params)
{
    return run_RKCL_consistency_prefix(
        g,
        B,
        params,
        g.n);
}

// Simple Result wrapper for experiments that do not need
// the detailed consistency trace.
inline Result
run_RKCL_prefix(
    const mygraph::tinyGraph& g,
    double B,
    const RKCLParams& params,
    std::size_t active_n)
{
    return run_RKCL_consistency_prefix(
        g,
        B,
        params,
        active_n)
        .final_result;
}

inline Result
run_RKCL(
    const mygraph::tinyGraph& g,
    double B,
    const RKCLParams& params)
{
    return run_RKCL_prefix(
        g,
        B,
        params,
        g.n);
}

} // namespace algs

#endif // ALGS_RKCL_IM_H

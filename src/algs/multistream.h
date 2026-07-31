// src/algs/multistream.h
#ifndef ALGS_MULTISTREAM_H
#define ALGS_MULTISTREAM_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
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

using MultiStreamClock = std::chrono::high_resolution_clock;
using mygraph::node_id;
using subm::Solution;

// =====================================================
// Common helpers.
// =====================================================
inline double multistream_node_cost(
    const mygraph::tinyGraph& g,
    node_id u)
{
    return subm_obj_common::node_cost(g, u);
}

inline double multistream_cost_of(
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
            c += multistream_node_cost(
                g,
                static_cast<node_id>(u));
        }
    }

    return c;
}

inline double multistream_eval_singleton(
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

// Default h used when the caller does not explicitly provide one.
// We keep h >= 2 because MultiStream uses
// 1 + 1 / (2^(h-1) - 1), which is undefined for h = 1.
inline int multistream_default_h(
    double eps)
{
    if (!(eps > 0.0 &&
          eps < 1.0))
    {
        throw std::invalid_argument(
            "MultiStream: eps must be in (0,1).");
    }

    const double x =
        1.0 / (8.0 * eps);

    double raw = 2.0;

    if (x > 1.0) {
        raw =
            std::log2(x) +
            1.0;
    }

    return std::max(
        2,
        static_cast<int>(
            std::ceil(raw)));
}

// =====================================================
// OneStream output used by MultiStream.
// M1 = Q*: feasible solution returned by OneStream.
// M2 = final cumulative set kept in memory.
// =====================================================
struct OneStreamOutput {
    Solution M1;
    Solution M2;

    double fM1 = 0.0;
    double fM2 = 0.0;

    double costM1 = 0.0;
    double costM2 = 0.0;

    std::uint64_t queries = 0;
};

// One nearly-feasible block S_i in Algorithm 1.
struct OneStreamBlock {
    std::vector<node_id> elems;
    double cost = 0.0;
};

// =====================================================
// OneStream on prefix V_t = {0,...,active_n-1}.
//
// This implementation follows Algorithm 1:
// 1) Maintain the cumulative union of a sliding window of blocks.
// 2) Accept e when
//      f(e | cumulative) / c(e) >= f(cumulative) / B.
// 3) Start a new block when current block cost reaches B.
// 4) When 2h blocks are present, delete the oldest h blocks.
// 5) Return Q* and the final cumulative set.
// =====================================================
inline OneStreamOutput
run_OneStream_prefix(
    const mygraph::tinyGraph& g,
    double B,
    int h,
    std::size_t active_n)
{
    if (B < 0.0) {
        throw std::invalid_argument(
            "OneStream: B must be non-negative.");
    }

    if (h < 1) {
        throw std::invalid_argument(
            "OneStream: h must be >= 1.");
    }

    OneStreamOutput out;

    const std::size_t n = g.n;

    active_n =
        std::min(
            active_n,
            n);

    out.M1.assign(n, 0);
    out.M2.assign(n, 0);

    if (active_n == 0 ||
        B <= 0.0)
    {
        out.fM1 =
            subm::sfunc_evaluate(
                g,
                out.M1);
        ++out.queries;

        out.fM2 = out.fM1;

        return out;
    }

    // Sliding window of candidate blocks S_j,...,S_i.
    std::deque<OneStreamBlock> blocks;
    blocks.emplace_back();

    // Explicit bitset for the current cumulative union.
    Solution cumulative(n, 0);

    double cumulative_cost = 0.0;

    double f_cumulative =
        subm::sfunc_evaluate(
            g,
            cumulative);
    ++out.queries;

    // Best feasible singleton e*.
    Solution singleton_tmp(n, 0);

    bool has_best_singleton = false;
    node_id best_singleton = 0;

    double best_singleton_value =
        -std::numeric_limits<double>::infinity();

    for (std::size_t idx = 0;
         idx < active_n;
         ++idx)
    {
        const node_id e =
            static_cast<node_id>(idx);

        const double ce =
            multistream_node_cost(
                g,
                e);

        // Elements with non-positive cost are ignored by this
        // experimental implementation. Elements with c(e) > B
        // can never be feasible and are discarded.
        if (ce <= 0.0 ||
            ce > B)
        {
            continue;
        }

        // Line 10: maintain the best feasible singleton.
        const double singleton_value =
            multistream_eval_singleton(
                g,
                singleton_tmp,
                e,
                out.queries);

        if (!has_best_singleton ||
            singleton_value >
                best_singleton_value)
        {
            has_best_singleton = true;
            best_singleton = e;
            best_singleton_value =
                singleton_value;
        }

        // Line 4: density test against the current cumulative set.
        const double delta =
            subm::sfunc_marginal(
                g,
                cumulative,
                e,
                f_cumulative);
        ++out.queries;

        const double lhs =
            delta / ce;

        const double rhs =
            f_cumulative / B;

        if (lhs < rhs) {
            continue;
        }

        // Line 5: add e to the current block S_i.
        blocks.back().elems.push_back(e);
        blocks.back().cost += ce;

        cumulative[e] = 1;
        cumulative_cost += ce;
        f_cumulative += delta;

        // Lines 6-9: current block becomes nearly feasible.
        if (blocks.back().cost >= B)
        {
            // If there are exactly 2h blocks, remove the oldest h.
            if (blocks.size() ==
                static_cast<std::size_t>(2 * h))
            {
                for (int r = 0;
                     r < h;
                     ++r)
                {
                    if (blocks.empty()) {
                        break;
                    }

                    const OneStreamBlock& old =
                        blocks.front();

                    for (node_id u :
                         old.elems)
                    {
                        cumulative[u] = 0;
                    }

                    cumulative_cost -=
                        old.cost;

                    blocks.pop_front();
                }

                // Deleting old blocks changes the cumulative
                // objective value. Re-evaluate it exactly.
                f_cumulative =
                    subm::sfunc_evaluate(
                        g,
                        cumulative);
                ++out.queries;
            }

            // Start S_{i+1} = empty.
            blocks.emplace_back();
        }
    }

    // M2 is the final cumulative set retained in memory.
    out.M2 =
        cumulative;

    out.fM2 =
        f_cumulative;

    out.costM2 =
        cumulative_cost;

    // Lines 12-15: extract Q*.
    if (cumulative_cost <= B)
    {
        out.M1 =
            cumulative;

        out.fM1 =
            f_cumulative;

        out.costM1 =
            cumulative_cost;

        return out;
    }

    // Reconstruct the chronological order of retained elements.
    std::vector<node_id> retained_order;

    std::size_t retained_size = 0;

    for (const OneStreamBlock& block :
         blocks)
    {
        retained_size +=
            block.elems.size();
    }

    retained_order.reserve(
        retained_size);

    for (const OneStreamBlock& block :
         blocks)
    {
        retained_order.insert(
            retained_order.end(),
            block.elems.begin(),
            block.elems.end());
    }

    // S(z): largest suffix whose total cost is <= B.
    Solution suffix(n, 0);
    double suffix_cost = 0.0;

    for (std::size_t pos =
             retained_order.size();
         pos > 0;
         --pos)
    {
        const node_id u =
            retained_order[pos - 1];

        const double cu =
            multistream_node_cost(
                g,
                u);

        if (suffix_cost + cu > B) {
            break;
        }

        suffix[u] = 1;
        suffix_cost += cu;
    }

    const double suffix_value =
        subm::sfunc_evaluate(
            g,
            suffix);
    ++out.queries;

    // Q* = better of S(z) and best singleton.
    if (has_best_singleton &&
        best_singleton_value >
            suffix_value)
    {
        out.M1.assign(n, 0);
        out.M1[best_singleton] = 1;

        out.fM1 =
            best_singleton_value;

        out.costM1 =
            multistream_node_cost(
                g,
                best_singleton);
    }
    else {
        out.M1 =
            std::move(suffix);

        out.fM1 =
            suffix_value;

        out.costM1 =
            suffix_cost;
    }

    return out;
}

inline OneStreamOutput
run_OneStream(
    const mygraph::tinyGraph& g,
    double B,
    int h)
{
    return run_OneStream_prefix(
        g,
        B,
        h,
        g.n);
}

// =====================================================
// MultiStream candidate A_rho.
// =====================================================
struct MultiStreamCandidate {
    double rho = 0.0;

    Solution S;

    double f_value = 0.0;
    double cost = 0.0;
};

// =====================================================
// Construct the threshold set
//
// P = {(1-eps)^(-z): z in Z and lower <= ... <= upper}
//
// We generate the exact geometric grid q^z,
// q = 1/(1-eps), rather than starting directly at lower.
// =====================================================
inline std::vector<double>
multistream_build_thresholds(
    double fM1,
    double fM2,
    double B,
    double eps,
    int h)
{
    if (!(eps > 0.0 &&
          eps < 1.0))
    {
        throw std::invalid_argument(
            "MultiStream: eps must be in (0,1).");
    }

    if (!(B > 0.0)) {
        return {};
    }

    if (h < 2) {
        throw std::invalid_argument(
            "MultiStream: h must be >= 2.");
    }

    const double denom =
        std::pow(
            2.0,
            static_cast<double>(h - 1)) -
        1.0;

    if (!(denom > 0.0)) {
        throw std::runtime_error(
            "MultiStream: invalid h in threshold bound.");
    }

    const double lower =
        (1.0 - eps) *
        fM1 /
        (2.0 * B);

    const double upper =
        (1.0 + 1.0 / denom) *
        fM2 /
        (eps * B);

    if (!(lower > 0.0) ||
        !(upper > 0.0) ||
        lower > upper)
    {
        return {};
    }

    const double q =
        1.0 /
        (1.0 - eps);

    const double logq =
        std::log(q);

    const double tol =
        1e-12;

    const long long z_min =
        static_cast<long long>(
            std::ceil(
                std::log(lower) /
                    logq -
                tol));

    const long long z_max =
        static_cast<long long>(
            std::floor(
                std::log(upper) /
                    logq +
                tol));

    if (z_min > z_max) {
        return {};
    }

    const long long count =
        z_max -
        z_min +
        1;

    // Safety guard against pathological numerical/input cases.
    if (count >
        1000000LL)
    {
        throw std::runtime_error(
            "MultiStream: threshold set is unexpectedly large.");
    }

    std::vector<double> P;

    P.reserve(
        static_cast<std::size_t>(
            count));

    for (long long z = z_min;
         z <= z_max;
         ++z)
    {
        const double rho =
            std::exp(
                static_cast<double>(z) *
                logq);

        if (rho + tol >= lower &&
            rho <= upper + tol)
        {
            P.push_back(rho);
        }
    }

    return P;
}

// =====================================================
// MultiStream on prefix V_t.
//
// Pass 1: OneStream -> M1, M2 and threshold set P.
// Pass 2: construct A_rho for every rho in P.
// Pass 3: test one-element augmentations A_rho U {e}.
// Final: greedily augment L* using elements already stored in A_rho.
// =====================================================
inline Result
run_MultiStream_prefix(
    const mygraph::tinyGraph& g,
    double B,
    double eps,
    int h,
    std::size_t active_n)
{
    Result res;

    res.algo =
        "MultiStream";

    {
        std::ostringstream oss;

        oss
            << "B="
            << B
            << ", eps="
            << eps
            << ", h="
            << h;

        res.constraint =
            oss.str();
    }

    res.inS.assign(
        g.n,
        0);

    res.queries = 0;

    const auto t0 =
        MultiStreamClock::now();

    active_n =
        std::min(
            active_n,
            g.n);

    if (!(eps > 0.0 &&
          eps < 1.0))
    {
        throw std::invalid_argument(
            "MultiStream: eps must be in (0,1).");
    }

    if (h < 2) {
        throw std::invalid_argument(
            "MultiStream: h must be >= 2.");
    }

    if (active_n == 0 ||
        B <= 0.0)
    {
        res.f_value =
            subm::sfunc_evaluate(
                g,
                res.inS);
        ++res.queries;

        res.time_sec =
            std::chrono::duration<double>(
                MultiStreamClock::now() -
                t0)
                .count();

        return res;
    }

    // =================================================
    // PASS 1: OneStream.
    // =================================================
    OneStreamOutput one =
        run_OneStream_prefix(
            g,
            B,
            h,
            active_n);

    res.queries +=
        one.queries;

    Solution Lstar =
        one.M1;

    double fLstar =
        one.fM1;

    double cLstar =
        one.costM1;

    // =================================================
    // Build P.
    // =================================================
    const std::vector<double> P =
        multistream_build_thresholds(
            one.fM1,
            one.fM2,
            B,
            eps,
            h);

    // Empty-set value for all A_rho.
    Solution empty(g.n, 0);

    const double f_empty =
        subm::sfunc_evaluate(
            g,
            empty);
    ++res.queries;

    std::vector<MultiStreamCandidate>
        candidates;

    candidates.reserve(
        P.size());

    for (double rho :
         P)
    {
        MultiStreamCandidate cand;

        cand.rho =
            rho;

        cand.S.assign(
            g.n,
            0);

        cand.f_value =
            f_empty;

        cand.cost =
            0.0;

        candidates.push_back(
            std::move(cand));
    }

    // =================================================
    // PASS 2: build every A_rho.
    // =================================================
    Solution singleton_tmp(
        g.n,
        0);

    for (std::size_t idx = 0;
         idx < active_n;
         ++idx)
    {
        const node_id e =
            static_cast<node_id>(idx);

        const double ce =
            multistream_node_cost(
                g,
                e);

        if (ce <= 0.0 ||
            ce > B)
        {
            continue;
        }

        // Line 5 prefilter: rho <= f(e)/c(e).
        const double singleton_value =
            multistream_eval_singleton(
                g,
                singleton_tmp,
                e,
                res.queries);

        const double singleton_density =
            singleton_value /
            ce;

        for (MultiStreamCandidate& cand :
             candidates)
        {
            if (cand.rho >
                singleton_density)
            {
                continue;
            }

            if (cand.cost + ce >
                B)
            {
                continue;
            }

            const double delta =
                subm::sfunc_marginal(
                    g,
                    cand.S,
                    e,
                    cand.f_value);
            ++res.queries;

            if (delta <
                cand.rho * ce)
            {
                continue;
            }

            // Lines 6-7.
            cand.S[e] = 1;
            cand.cost += ce;
            cand.f_value += delta;

            // Lines 8-9.
            if (cand.f_value >
                fLstar)
            {
                Lstar =
                    cand.S;

                fLstar =
                    cand.f_value;

                cLstar =
                    cand.cost;
            }
        }
    }

    // =================================================
    // PASS 3: try A_rho U {e}.
    // =================================================
    for (std::size_t idx = 0;
         idx < active_n;
         ++idx)
    {
        const node_id e =
            static_cast<node_id>(idx);

        const double ce =
            multistream_node_cost(
                g,
                e);

        if (ce <= 0.0 ||
            ce > B)
        {
            continue;
        }

        for (const MultiStreamCandidate& cand :
             candidates)
        {
            if (cand.S[e]) {
                continue;
            }

            if (cand.cost + ce >
                B)
            {
                continue;
            }

            const double delta =
                subm::sfunc_marginal(
                    g,
                    cand.S,
                    e,
                    cand.f_value);
            ++res.queries;

            const double candidate_value =
                cand.f_value +
                delta;

            // Algorithm 2 uses >= here.
            if (candidate_value >=
                fLstar)
            {
                Lstar =
                    cand.S;

                Lstar[e] = 1;

                fLstar =
                    candidate_value;

                cLstar =
                    cand.cost +
                    ce;
            }
        }
    }

    // =================================================
    // Final augmentation:
    // for each rho and each e in A_rho,
    // add e to L* whenever budget remains.
    //
    // This step assumes the objective is monotone,
    // as in the MultiStream paper.
    // =================================================
    for (const MultiStreamCandidate& cand :
         candidates)
    {
        for (std::size_t idx = 0;
             idx < active_n;
             ++idx)
        {
            if (!cand.S[idx] ||
                Lstar[idx])
            {
                continue;
            }

            const node_id e =
                static_cast<node_id>(idx);

            const double ce =
                multistream_node_cost(
                    g,
                    e);

            if (ce <= 0.0 ||
                cLstar + ce >
                    B)
            {
                continue;
            }

            const double delta =
                subm::sfunc_marginal(
                    g,
                    Lstar,
                    e,
                    fLstar);
            ++res.queries;

            Lstar[e] = 1;
            cLstar += ce;
            fLstar += delta;
        }
    }

    res.inS =
        std::move(Lstar);

    res.f_value =
        fLstar;

    res.time_sec =
        std::chrono::duration<double>(
            MultiStreamClock::now() -
            t0)
            .count();

    return res;
}

inline Result
run_MultiStream(
    const mygraph::tinyGraph& g,
    double B,
    double eps,
    int h)
{
    return run_MultiStream_prefix(
        g,
        B,
        eps,
        h,
        g.n);
}

// Convenience overload using the default h.
inline Result
run_MultiStream(
    const mygraph::tinyGraph& g,
    double B,
    double eps)
{
    return run_MultiStream(
        g,
        B,
        eps,
        multistream_default_h(
            eps));
}

// =====================================================
// Consistency statistics.
//
// IMPORTANT:
// To stay compatible with the current EDL experiment,
// the theory-style consistency metric is
//
//   change_t = |S_t \ S_{t-1}|
//
// and
//
//   C = max_{t>=2} change_t.
//
// We additionally record symmetric difference as a
// separate diagnostic:
//
//   symdiff_t = |S_t Delta S_{t-1}|.
//
// This routine is an empirical prefix-rerun measurement:
// S_t = MultiStream(V_t).
// =====================================================
struct MultiStreamConsistencyStats {
    Result final_result;

    std::vector<std::uint64_t>
        per_step;

    std::vector<std::uint64_t>
        per_step_symmetric_difference;

    std::uint64_t C = 0;
    std::uint64_t cumulative_consistency = 0;
    std::uint64_t initial_changes = 0;

    std::uint64_t max_symmetric_difference = 0;
    std::uint64_t cumulative_symmetric_difference = 0;

    std::uint64_t total_queries = 0;
    double total_time_sec = 0.0;
};

inline std::uint64_t
multistream_consistency_change(
    const Solution& prev,
    const Solution& cur)
{
    if (prev.size() !=
        cur.size())
    {
        throw std::invalid_argument(
            "MultiStream consistency: solution sizes differ.");
    }

    std::uint64_t changes = 0;

    for (std::size_t u = 0;
         u < cur.size();
         ++u)
    {
        if (cur[u] &&
            !prev[u])
        {
            ++changes;
        }
    }

    return changes;
}

inline std::uint64_t
multistream_symmetric_difference(
    const Solution& prev,
    const Solution& cur)
{
    if (prev.size() !=
        cur.size())
    {
        throw std::invalid_argument(
            "MultiStream symdiff: solution sizes differ.");
    }

    std::uint64_t changes = 0;

    for (std::size_t u = 0;
         u < cur.size();
         ++u)
    {
        if (static_cast<bool>(prev[u]) !=
            static_cast<bool>(cur[u]))
        {
            ++changes;
        }
    }

    return changes;
}

inline MultiStreamConsistencyStats
run_MultiStream_consistency(
    const mygraph::tinyGraph& g,
    double B,
    double eps,
    int h)
{
    MultiStreamConsistencyStats stats;

    stats.per_step.assign(
        g.n,
        0);

    stats.per_step_symmetric_difference.assign(
        g.n,
        0);

    Solution prev(
        g.n,
        0);

    for (std::size_t t = 1;
         t <= g.n;
         ++t)
    {
        Result cur =
            run_MultiStream_prefix(
                g,
                B,
                eps,
                h,
                t);

        const std::uint64_t change =
            multistream_consistency_change(
                prev,
                cur.inS);

        const std::uint64_t symdiff =
            multistream_symmetric_difference(
                prev,
                cur.inS);

        if (t == 1)
        {
            stats.initial_changes =
                change;

            stats.per_step[0] =
                0;

            stats.per_step_symmetric_difference[0] =
                0;
        }
        else {
            stats.per_step[t - 1] =
                change;

            stats.per_step_symmetric_difference[t - 1] =
                symdiff;

            stats.cumulative_consistency +=
                change;

            stats.C =
                std::max(
                    stats.C,
                    change);

            stats.cumulative_symmetric_difference +=
                symdiff;

            stats.max_symmetric_difference =
                std::max(
                    stats.max_symmetric_difference,
                    symdiff);
        }

        stats.total_queries +=
            cur.queries;

        stats.total_time_sec +=
            cur.time_sec;

        prev =
            cur.inS;

        stats.final_result =
            std::move(cur);
    }

    return stats;
}

inline MultiStreamConsistencyStats
run_MultiStream_consistency(
    const mygraph::tinyGraph& g,
    double B,
    double eps)
{
    return run_MultiStream_consistency(
        g,
        B,
        eps,
        multistream_default_h(
            eps));
}

} // namespace algs

#endif // ALGS_MULTISTREAM_H

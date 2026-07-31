// src/algs/edl.h
#ifndef ALGS_EDL_H
#define ALGS_EDL_H

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

using Clock = std::chrono::high_resolution_clock;
using mygraph::node_id;
using subm::Solution;

// =====================================================
// Scalar node cost for EDL.
// =====================================================
inline double edl_node_cost(
    const mygraph::tinyGraph& g,
    node_id u)
{
    return subm_obj_common::node_cost(g, u);
}

inline double edl_cost_of(
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
            c +=
                edl_node_cost(
                    g,
                    static_cast<node_id>(u));
        }
    }

    return c;
}

inline double edl_eval_singleton(
    const mygraph::tinyGraph& g,
    Solution& tmp,
    node_id e,
    std::uint64_t& queries)
{
    tmp[e] = 1;

    const double val =
        subm::sfunc_evaluate(g, tmp);

    ++queries;
    tmp[e] = 0;

    return val;
}

// =====================================================
// LA preprocessing restricted to the arriving prefix
// V_t = {0,...,active_n-1}.
// The Solution vector still has size g.n.
// =====================================================
inline Solution run_LA_preprocess_prefix(
    const mygraph::tinyGraph& g,
    double B,
    std::uint64_t& queries,
    double& out_fS,
    std::size_t active_n)
{
    if (B < 0.0) {
        throw std::invalid_argument(
            "LA: B must be non-negative.");
    }

    const std::size_t n = g.n;
    active_n = std::min(active_n, n);

    Solution X(n, 0);
    Solution Y(n, 0);

    std::vector<node_id> Xseq;
    std::vector<node_id> Yseq;

    Xseq.reserve(active_n);
    Yseq.reserve(active_n);

    double fX =
        subm::sfunc_evaluate(g, X);
    ++queries;

    double fY =
        subm::sfunc_evaluate(g, Y);
    ++queries;

    double cX = 0.0;
    double cY = 0.0;

    // Best feasible singleton.
    Solution tmp(n, 0);

    double best_single_val =
        -std::numeric_limits<double>::infinity();

    node_id best_single = 0;
    bool has_single = false;

    for (std::size_t idx = 0;
         idx < active_n;
         ++idx)
    {
        const node_id e =
            static_cast<node_id>(idx);

        const double ce =
            edl_node_cost(g, e);

        if (ce <= 0.0 || ce > B) {
            continue;
        }

        const double val =
            edl_eval_singleton(
                g,
                tmp,
                e,
                queries);

        if (!has_single ||
            val > best_single_val)
        {
            has_single = true;
            best_single_val = val;
            best_single = e;
        }
    }

    Solution S_single(n, 0);

    if (has_single) {
        S_single[best_single] = 1;
    }

    // Elements with c(e) <= B/2.
    const double halfB = B / 2.0;

    for (std::size_t idx = 0;
         idx < active_n;
         ++idx)
    {
        const node_id e =
            static_cast<node_id>(idx);

        const double ce =
            edl_node_cost(g, e);

        if (ce <= 0.0 ||
            ce > halfB)
        {
            continue;
        }

        const double densX =
            (cX > 0.0)
                ? (fX / cX)
                : 0.0;

        const double densY =
            (cY > 0.0)
                ? (fY / cY)
                : 0.0;

        const double dX =
            subm::sfunc_marginal(
                g,
                X,
                e,
                fX);
        ++queries;

        const double dY =
            subm::sfunc_marginal(
                g,
                Y,
                e,
                fY);
        ++queries;

        const double mdX = dX / ce;
        const double mdY = dY / ce;

        const bool canX =
            mdX >= densX;

        const bool canY =
            mdY >= densY;

        if (!canX &&
            !canY)
        {
            continue;
        }

        if (canX &&
            (!canY || mdX >= mdY))
        {
            X[e] = 1;
            Xseq.push_back(e);
            cX += ce;
            fX += dX;
        }
        else {
            Y[e] = 1;
            Yseq.push_back(e);
            cY += ce;
            fY += dY;
        }
    }

    // Best feasible suffix of a sequence.
    auto best_suffix =
        [&](const std::vector<node_id>& seq,
            double& best_val,
            std::size_t& best_j)
    {
        best_val =
            -std::numeric_limits<double>::infinity();

        best_j = 0;

        Solution cur(n, 0);
        double cur_cost = 0.0;

        best_val =
            subm::sfunc_evaluate(
                g,
                cur);
        ++queries;

        for (std::size_t t = 0;
             t < seq.size();
             ++t)
        {
            const node_id u =
                seq[seq.size() - 1 - t];

            const double cu =
                edl_node_cost(
                    g,
                    u);

            if (cu <= 0.0) {
                continue;
            }

            cur_cost += cu;

            if (cur_cost > B) {
                break;
            }

            cur[u] = 1;

            const double val =
                subm::sfunc_evaluate(
                    g,
                    cur);
            ++queries;

            const std::size_t j =
                t + 1;

            if (val > best_val) {
                best_val = val;
                best_j = j;
            }
        }
    };

    double bestX_val = 0.0;
    double bestY_val = 0.0;

    std::size_t bestX_j = 0;
    std::size_t bestY_j = 0;

    best_suffix(
        Xseq,
        bestX_val,
        bestX_j);

    best_suffix(
        Yseq,
        bestY_val,
        bestY_j);

    auto build_suffix_solution =
        [&](const std::vector<node_id>& seq,
            std::size_t j)
    {
        Solution sol(n, 0);

        if (j > seq.size()) {
            throw std::runtime_error(
                "LA: invalid suffix size.");
        }

        for (std::size_t t = 0;
             t < j;
             ++t)
        {
            const node_id u =
                seq[seq.size() - j + t];

            sol[u] = 1;
        }

        return sol;
    };

    Solution Xp =
        build_suffix_solution(
            Xseq,
            bestX_j);

    Solution Yp =
        build_suffix_solution(
            Yseq,
            bestY_j);

    Solution S = Xp;
    double fS = bestX_val;

    if (bestY_val > fS) {
        S = Yp;
        fS = bestY_val;
    }

    if (has_single &&
        best_single_val > fS)
    {
        S = S_single;
        fS = best_single_val;
    }

    out_fS = fS;
    return S;
}

inline Solution run_LA_preprocess(
    const mygraph::tinyGraph& g,
    double B,
    std::uint64_t& queries,
    double& out_fS)
{
    return run_LA_preprocess_prefix(
        g,
        B,
        queries,
        out_fS,
        g.n);
}

// =====================================================
// EDL on prefix V_t.
// =====================================================
inline Result run_EDL_prefix(
    const mygraph::tinyGraph& g,
    double B,
    double eps,
    std::size_t active_n)
{
    Result res;
    res.algo = "EDL";

    {
        std::ostringstream oss;
        oss
            << "B="
            << B
            << ", eps="
            << eps;

        res.constraint =
            oss.str();
    }

    res.inS.assign(g.n, 0);
    res.queries = 0;

    const auto t0 =
        Clock::now();

    active_n =
        std::min(
            active_n,
            g.n);

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
                Clock::now() - t0)
                .count();

        return res;
    }

    if (!(eps > 0.0 &&
          eps < 1.0))
    {
        throw std::invalid_argument(
            "EDL: eps must be in (0,1).");
    }

    // Step 1: LA preprocessing.
    double M = 0.0;

    Solution S_prime =
        run_LA_preprocess_prefix(
            g,
            B,
            res.queries,
            M,
            active_n);

    (void)S_prime;

    const double eps_p =
        eps / 14.0;

    if (!(eps_p > 0.0 &&
          eps_p < 1.0))
    {
        throw std::invalid_argument(
            "EDL: eps/14 must be in (0,1).");
    }

    // Step 2: EDL construction.
    Solution X(g.n, 0);
    Solution Y(g.n, 0);

    double fX =
        subm::sfunc_evaluate(
            g,
            X);
    ++res.queries;

    double fY =
        subm::sfunc_evaluate(
            g,
            Y);
    ++res.queries;

    double cX = 0.0;
    double cY = 0.0;

    const double base =
        1.0 /
        (1.0 - eps_p);

    const double target =
        19.0 /
        (eps_p * eps_p);

    int Imax = 0;

    if (target <= 1.0) {
        Imax = 1;
    }
    else {
        Imax =
            static_cast<int>(
                std::floor(
                    std::log(target) /
                    std::log(base))) +
            1;

        if (Imax < 0) {
            Imax = 0;
        }
    }

    for (int i = 0;
         i <= Imax;
         ++i)
    {
        const double theta =
            (19.0 *
             M *
             std::pow(
                 1.0 - eps_p,
                 static_cast<double>(i))) /
            (5.0 *
             eps_p *
             B);

        for (std::size_t idx = 0;
             idx < active_n;
             ++idx)
        {
            const node_id e =
                static_cast<node_id>(idx);

            if (X[e] ||
                Y[e])
            {
                continue;
            }

            const double ce =
                edl_node_cost(
                    g,
                    e);

            if (ce <= 0.0 ||
                ce > B)
            {
                continue;
            }

            bool chooseX = false;
            bool chooseY = false;

            double best_md =
                -std::numeric_limits<double>::infinity();

            double best_delta = 0.0;

            if (cX + ce <= B) {
                const double d =
                    subm::sfunc_marginal(
                        g,
                        X,
                        e,
                        fX);
                ++res.queries;

                const double md =
                    d / ce;

                if (md >= theta &&
                    md > best_md)
                {
                    best_md = md;
                    best_delta = d;
                    chooseX = true;
                    chooseY = false;
                }
            }

            if (cY + ce <= B) {
                const double d =
                    subm::sfunc_marginal(
                        g,
                        Y,
                        e,
                        fY);
                ++res.queries;

                const double md =
                    d / ce;

                if (md >= theta &&
                    md > best_md)
                {
                    best_md = md;
                    best_delta = d;
                    chooseX = false;
                    chooseY = true;
                }
            }

            if (!(chooseX ||
                  chooseY))
            {
                continue;
            }

            if (chooseX) {
                X[e] = 1;
                cX += ce;
                fX += best_delta;
            }
            else {
                Y[e] = 1;
                cY += ce;
                fY += best_delta;
            }
        }
    }

    const double valX =
        subm::sfunc_evaluate(
            g,
            X);
    ++res.queries;

    const double valY =
        subm::sfunc_evaluate(
            g,
            Y);
    ++res.queries;

    if (valY > valX) {
        res.inS =
            std::move(Y);

        res.f_value =
            valY;
    }
    else {
        res.inS =
            std::move(X);

        res.f_value =
            valX;
    }

    const auto t1 =
        Clock::now();

    res.time_sec =
        std::chrono::duration<double>(
            t1 - t0)
            .count();

    return res;
}

inline Result run_EDL(
    const mygraph::tinyGraph& g,
    double B,
    double eps)
{
    return run_EDL_prefix(
        g,
        B,
        eps,
        g.n);
}

// =====================================================
// Consistency statistics.
//
// change_t = |S_t \ S_{t-1}|
// C        = max_{t>=2} change_t
// cumulative_consistency = sum_{t>=2} change_t
//
// This is an empirical prefix-rerun measurement.
// =====================================================
struct EDLConsistencyStats {
    Result final_result;

    std::vector<std::uint64_t>
        per_step;

    std::uint64_t C = 0;
    std::uint64_t cumulative_consistency = 0;
    std::uint64_t initial_changes = 0;

    std::uint64_t total_queries = 0;
    double total_time_sec = 0.0;
};

inline std::uint64_t edl_consistency_change(
    const Solution& prev,
    const Solution& cur)
{
    if (prev.size() !=
        cur.size())
    {
        throw std::invalid_argument(
            "EDL consistency: solution sizes differ.");
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

inline EDLConsistencyStats
run_EDL_consistency(
    const mygraph::tinyGraph& g,
    double B,
    double eps)
{
    EDLConsistencyStats stats;

    stats.per_step.assign(
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
            run_EDL_prefix(
                g,
                B,
                eps,
                t);

        const std::uint64_t change =
            edl_consistency_change(
                prev,
                cur.inS);

        if (t == 1) {
            stats.initial_changes =
                change;

            stats.per_step[0] =
                0;
        }
        else {
            stats.per_step[t - 1] =
                change;

            stats.cumulative_consistency +=
                change;

            stats.C =
                std::max(
                    stats.C,
                    change);
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

} // namespace algs

#endif // ALGS_EDL_H

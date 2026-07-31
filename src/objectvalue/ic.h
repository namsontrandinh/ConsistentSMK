// src/objectvalue/ic.h
#ifndef OBJECTVALUE_IC_H
#define OBJECTVALUE_IC_H

#include "objectvalue_common.h"
#include "sfunctions.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <utility>
#include <vector>

#ifndef _OPENMP
#error "IC requires OpenMP. Compile and link with -fopenmp."
#endif

#include <omp.h>

namespace subm_obj_ic {

using mygraph::edge_id;
using mygraph::node_id;

namespace detail_ic {

inline double clamp01(double p) {
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return 1.0;
    return p;
}

inline bool env_exists(const char* name) {
    return std::getenv(name) != nullptr;
}

inline std::size_t read_env_size_t(
    const char* name,
    std::size_t defv)
{
    if (const char* s = std::getenv(name)) {
        char* end = nullptr;
        const unsigned long long v =
            std::strtoull(s, &end, 10);

        if (end && *end == '\0') {
            return static_cast<std::size_t>(v);
        }
    }

    return defv;
}

inline std::uint64_t read_env_u64(
    const char* name,
    std::uint64_t defv)
{
    if (const char* s = std::getenv(name)) {
        char* end = nullptr;
        const unsigned long long v =
            std::strtoull(s, &end, 10);

        if (end && *end == '\0') {
            return static_cast<std::uint64_t>(v);
        }
    }

    return defv;
}

inline double read_env_double(
    const char* name,
    double defv)
{
    if (const char* s = std::getenv(name)) {
        char* end = nullptr;
        const double v =
            std::strtod(s, &end);

        if (end && *end == '\0') {
            return v;
        }
    }

    return defv;
}

inline std::uint64_t splitmix64(
    std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x =
        (x ^ (x >> 30)) *
        0xbf58476d1ce4e5b9ULL;
    x =
        (x ^ (x >> 27)) *
        0x94d049bb133111ebULL;

    return x ^ (x >> 31);
}

inline std::uint64_t runtime_seed64() {
    std::random_device rd;

    const std::uint64_t r =
        (static_cast<std::uint64_t>(rd()) << 32) ^
        static_cast<std::uint64_t>(rd());

    const std::uint64_t t =
        static_cast<std::uint64_t>(
            std::chrono::high_resolution_clock::now()
                .time_since_epoch()
                .count());

    const std::uint64_t addr_mix =
        static_cast<std::uint64_t>(
            reinterpret_cast<std::uintptr_t>(&rd));

    return splitmix64(r ^ t ^ addr_mix);
}

// Cache scalar IC adjacency:
// adjs[u] = {(v, p(u,v)), ...}.
struct ICCache {
    const mygraph::tinyGraph* gptr = nullptr;
    std::size_t n = 0;
    std::size_t m = 0;
    bool undirected = true;

    std::vector<
        std::vector<
            std::pair<node_id, double>>> adjs;

    void rebuild(
        const mygraph::tinyGraph& g)
    {
        gptr = &g;
        n = g.n;
        m = g.m;
        undirected = g.undirected;

        adjs.assign(n, {});

        for (edge_id eid = 0;
             eid < g.edges.size();
             ++eid)
        {
            const auto& E = g.edges[eid];
            const node_id u = E.u;
            const node_id v = E.v;

            if (u >= n || v >= n) {
                continue;
            }

            const double p =
                clamp01(
                    subm_obj_common::edge_weight(
                        g,
                        eid));

            if (p <= 0.0) {
                continue;
            }

            // Directed edge u -> v.
            adjs[u].push_back({v, p});

            // For an undirected graph, also add v -> u.
            if (undirected) {
                adjs[v].push_back({u, p});
            }
        }
    }
};

inline const std::vector<
    std::vector<
        std::pair<node_id, double>>>&
get_adjs(
    const mygraph::tinyGraph& g)
{
    static ICCache cache;

    if (cache.gptr != &g ||
        cache.n != g.n ||
        cache.m != g.m ||
        cache.undirected != g.undirected)
    {
        cache.rebuild(g);
    }

    return cache.adjs;
}

} // namespace detail_ic

// IC objective:
//   f(S) = E[|A(S)|] - lambda |S|
//
// Environment variables:
//   KIC_MC      default 100
//   KIC_SEED    fixed seed if set; runtime-random otherwise
//   KIC_LAMBDA  default 1.0
inline double evaluate(
    const mygraph::tinyGraph& g,
    const subm::Solution& inS)
{
    const std::size_t n = g.n;

    if (n == 0) {
        return 0.0;
    }

    subm::solution_check_size(g, inS);

    std::vector<node_id> seeds;
    seeds.reserve(n);

    std::size_t S_card = 0;

    for (node_id u = 0;
         u < static_cast<node_id>(n);
         ++u)
    {
        if (inS[u]) {
            ++S_card;
            seeds.push_back(u);
        }
    }

    const std::size_t mc =
        detail_ic::read_env_size_t(
            "KIC_MC",
            100);

    const double lambda =
        detail_ic::read_env_double(
            "KIC_LAMBDA",
            1.0);

    if (mc == 0) {
        return
            -lambda *
            static_cast<double>(S_card);
    }

    const std::uint64_t base_seed =
        detail_ic::env_exists("KIC_SEED")
            ? detail_ic::read_env_u64(
                  "KIC_SEED",
                  42ULL)
            : detail_ic::runtime_seed64();

    const auto& adjs =
        detail_ic::get_adjs(g);

    double sum_spread = 0.0;

#pragma omp parallel
    {
        std::vector<std::uint32_t> seen(n, 0);
        std::uint32_t stamp = 1;

        auto bump_stamp =
            [](std::uint32_t& st,
               std::vector<std::uint32_t>& arr)
        {
            ++st;

            if (st == 0) {
                std::fill(
                    arr.begin(),
                    arr.end(),
                    0);
                st = 1;
            }
        };

        std::uniform_real_distribution<double>
            uni(0.0, 1.0);

#pragma omp for schedule(static) reduction(+:sum_spread)
        for (std::size_t it = 0;
             it < mc;
             ++it)
        {
            std::mt19937_64 rng(
                detail_ic::splitmix64(
                    base_seed ^
                    static_cast<std::uint64_t>(it)));

            bump_stamp(stamp, seen);

            std::size_t activated_cnt = 0;

            std::vector<node_id> frontier;
            frontier.reserve(seeds.size());

            for (node_id s : seeds) {
                if (s >= n) {
                    continue;
                }

                if (seen[s] == stamp) {
                    continue;
                }

                seen[s] = stamp;
                frontier.push_back(s);
                ++activated_cnt;
            }

            while (!frontier.empty()) {
                std::vector<node_id> next;
                next.reserve(frontier.size());

                for (node_id u : frontier) {
                    for (const auto& vp : adjs[u]) {
                        const node_id v = vp.first;
                        const double p = vp.second;

                        if (v >= n) {
                            continue;
                        }

                        if (seen[v] == stamp) {
                            continue;
                        }

                        if (uni(rng) < p) {
                            seen[v] = stamp;
                            next.push_back(v);
                            ++activated_cnt;
                        }
                    }
                }

                frontier.swap(next);
            }

            sum_spread +=
                static_cast<double>(
                    activated_cnt);
        }
    }

    const double expected_spread =
        sum_spread /
        static_cast<double>(mc);

    return
        expected_spread -
        lambda *
        static_cast<double>(S_card);
}

inline double marginal(
    const mygraph::tinyGraph& g,
    const subm::Solution& inS,
    mygraph::node_id x,
    double fS)
{
    subm::solution_check_size(g, inS);

    if (x >= g.n) {
        return 0.0;
    }

    if (inS[x]) {
        return 0.0;
    }

    subm::Solution inS_after = inS;
    inS_after[x] = 1;

    const double after =
        evaluate(g, inS_after);

    return after - fS;
}

} // namespace subm_obj_ic

#endif // OBJECTVALUE_IC_H

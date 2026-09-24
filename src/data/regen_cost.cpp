// src/data/regen_cost.cpp
//
// Reads an existing .bin graph and OVERWRITES node cost
// (nodes[i].weights[0], the field subm_obj_common::node_cost reads)
// with a fresh i.i.d. draw from Uniform[cost_min, cost_max].
// Everything else (edges, weights, alpha, part_id, K) is preserved
// unchanged. Written for the CSONET-2026 camera-ready experiment
// request: "cost sinh ngẫu nhiên trong khoảng 10-20".
//
// Usage:
//   ./regen_cost input.bin output.bin [cost_min=10] [cost_max=20] [seed=42]

#include <iostream>
#include <string>
#include <cstdlib>
#include <random>

#include "mygraph.h"

int main(int argc, char** argv) {
    using namespace mygraph;

    if (argc < 3) {
        std::cerr
            << "Usage:\n  " << argv[0]
            << " input.bin output.bin [cost_min=10] [cost_max=20] [seed=42]\n";
        return 1;
    }

    const std::string input_bin  = argv[1];
    const std::string output_bin = argv[2];

    double cost_min = (argc >= 4) ? std::strtod(argv[3], nullptr) : 10.0;
    double cost_max = (argc >= 5) ? std::strtod(argv[4], nullptr) : 20.0;
    unsigned seed    = (argc >= 6) ? static_cast<unsigned>(std::strtoul(argv[5], nullptr, 10)) : 42u;

    if (!(cost_min > 0.0) || !(cost_max > cost_min)) {
        std::cerr << "Error: require 0 < cost_min < cost_max.\n";
        return 1;
    }

    tinyGraph g;
    if (!g.read_binary(input_bin)) {
        std::cerr << "Error: cannot read graph: " << input_bin << "\n";
        return 1;
    }

    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(cost_min, cost_max);

    double sum = 0.0, mn = 1e300, mx = -1e300;

    for (std::size_t i = 0; i < g.n; ++i) {
        // K may be > 1 in general, but the cost read by
        // subm_obj_common::node_cost is always weights[0].
        // Only that slot needs the new draw; other topics (if any)
        // are left untouched.
        if (g.nodes[i].weights.empty()) {
            g.nodes[i].weights.assign(g.K > 0 ? g.K : 1, 1.0);
        }

        const double c = dist(gen);
        g.nodes[i].weights[0] = c;

        sum += c;
        mn = std::min(mn, c);
        mx = std::max(mx, c);
    }

    if (!g.write_binary(output_bin)) {
        std::cerr << "Error: write_binary failed: " << output_bin << "\n";
        return 1;
    }

    std::cout
        << "Done. Regenerated node cost for " << g.n << " nodes.\n"
        << "  input       = " << input_bin << "\n"
        << "  output      = " << output_bin << "\n"
        << "  cost range  = [" << cost_min << ", " << cost_max << "]\n"
        << "  seed        = " << seed << "\n"
        << "  actual min  = " << mn << "\n"
        << "  actual max  = " << mx << "\n"
        << "  actual mean = " << (sum / static_cast<double>(g.n)) << "\n";

    return 0;
}

// src/sfunctions.h
#ifndef SFUNCTIONS_H
#define SFUNCTIONS_H

#include <cstdint>
#include <stdexcept>
#include <vector>

#include "mygraph.h"

namespace subm {

using mygraph::node_id;

// Binary indicator vector of a selected set S.
// inS[u] = 1 iff u belongs to S; otherwise 0.
using Solution = std::vector<std::uint8_t>;

inline void solution_check_size(
    const mygraph::tinyGraph& g,
    const Solution& inS)
{
    if (inS.size() != g.n) {
        throw std::invalid_argument(
            "subm::Solution size does not match graph size.");
    }
}

// Implemented in sfunctions_impl.h after compile-time objective dispatch.
inline double sfunc_evaluate(
    const mygraph::tinyGraph& g,
    const Solution& inS);

inline double sfunc_marginal(
    const mygraph::tinyGraph& g,
    const Solution& inS,
    node_id u,
    double fS);

} // namespace subm

#endif // SFUNCTIONS_H

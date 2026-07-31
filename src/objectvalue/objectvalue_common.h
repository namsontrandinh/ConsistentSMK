// src/objectvalue/objectvalue_common.h
#ifndef OBJECTVALUE_COMMON_H
#define OBJECTVALUE_COMMON_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "mygraph.h"

namespace subm_obj_common {

using mygraph::edge_id;
using mygraph::node_id;

// Scalar node cost used by knapsack algorithms.
// Current tinyGraph stores node values in a vector;
// for the single-topic experiments, cost c(u)=weights[0].
inline double node_cost(
    const mygraph::tinyGraph& g,
    node_id u)
{
    if (u >= g.n) {
        return 0.0;
    }

    if (g.nodes[u].weights.empty()) {
        return 0.0;
    }

    return g.nodes[u].weights[0];
}

// Scalar edge weight / propagation probability.
// For the single-topic experiments, w(e)=weights[0].
inline double edge_weight(
    const mygraph::tinyGraph& g,
    edge_id eid)
{
    if (eid >= g.edges.size()) {
        return 0.0;
    }

    const auto& E = g.edges[eid];

    if (E.weights.empty()) {
        return 0.0;
    }

    return E.weights[0];
}

} // namespace subm_obj_common

#endif // OBJECTVALUE_COMMON_H

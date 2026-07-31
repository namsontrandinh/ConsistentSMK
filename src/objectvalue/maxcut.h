// src/objectvalue/maxcut.h  (submodular Max-Cut on a graph)
#ifndef OBJECTVALUE_MAXCUT_H
#define OBJECTVALUE_MAXCUT_H

#include "objectvalue_common.h"

namespace subm_obj_maxcut {

// f(S) = sum_{(u,v) in E} w_uv * [inS[u] != inS[v]]
inline double evaluate(const mygraph::tinyGraph &g, const subm::Solution &inS) {
    double val = 0.0;
    for (mygraph::edge_id e = 0; e < g.edges.size(); ++e) {
        const auto &E = g.edges[e];
        if ((inS[E.u] != 0) ^ (inS[E.v] != 0)) val += E.weight;
    }
    return val;
}

// Δ(u|S) for adding u (inS[u]==0)
// Undirected: scan incident[u] (contains all adjacent edges)
// Directed: scan incident[u] (out) and incoming[u] (in)
inline double marginal(const mygraph::tinyGraph &g,
                       const subm::Solution &inS,
                       mygraph::node_id u,
                       double /*fS*/)
{
    double delta = 0.0;

    auto contrib_one_edge = [&](const mygraph::Edge &E, mygraph::node_id other) {
        // before: inS[u]=0 => cut indicator = inS[other]
        // after : inS[u]=1 => cut indicator = 1 - inS[other]
        // delta = (1 - inS[other]) - inS[other] = 1 - 2*inS[other]
        const double s = (inS[other] ? 1.0 : 0.0);
        delta += E.weight * (1.0 - 2.0 * s);
    };

    // out / incident edges
    for (auto eid : g.incident[u]) {
        const auto &E = g.edges[eid];
        mygraph::node_id other;
        if (g.undirected) {
            other = (E.u == u) ? E.v : E.u;
            contrib_one_edge(E, other);
        } else {
            // directed: incident[u] are out-edges u->v
            other = E.v;
            contrib_one_edge(E, other);
        }
    }

    if (!g.undirected) {
        // in-edges ? -> u
        for (auto eid : g.incoming[u]) {
            const auto &E = g.edges[eid];
            // incoming[u] edges have v==u
            mygraph::node_id other = E.u;
            contrib_one_edge(E, other);
        }
    }

    return delta;
}

} // namespace subm_obj_maxcut

#endif // OBJECTVALUE_MAXCUT_H

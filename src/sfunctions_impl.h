// src/sfunctions_impl.h
#ifndef SFUNCTIONS_IMPL_H
#define SFUNCTIONS_IMPL_H

#include <stdexcept>

#include "sfunctions.h"

// Select exactly one scalar submodular objective at build time.
#if defined(SFUNC_MAXCUT)
  #include "objectvalue/maxcut.h"
  namespace obj = subm_obj_maxcut;
#elif defined(SFUNC_REVENUE)
  #include "objectvalue/revenue.h"
  namespace obj = subm_obj_revenue;
#elif defined(SFUNC_IC)
  #include "objectvalue/ic.h"
  namespace obj = subm_obj_ic;
#else
  #error "No objective selected. Define one of: -DSFUNC_MAXCUT, -DSFUNC_REVENUE, -DSFUNC_IC"
#endif

namespace subm {

inline double sfunc_evaluate(
    const mygraph::tinyGraph& g,
    const Solution& inS)
{
    solution_check_size(g, inS);
    return obj::evaluate(g, inS);
}

inline double sfunc_marginal(
    const mygraph::tinyGraph& g,
    const Solution& inS,
    node_id u,
    double fS)
{
    solution_check_size(g, inS);

    if (u >= g.n) {
        throw std::out_of_range(
            "subm::sfunc_marginal: node out of range.");
    }

    if (inS[u]) {
        return 0.0;
    }

    return obj::marginal(g, inS, u, fS);
}

} // namespace subm

#endif // SFUNCTIONS_IMPL_H

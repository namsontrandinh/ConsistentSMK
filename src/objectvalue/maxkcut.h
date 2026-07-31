// src/objectvalue/maxkcut.h
#ifndef KSUB_OBJECTVALUE_MAXKCUT_H
#define KSUB_OBJECTVALUE_MAXKCUT_H

#include "objectvalue_common.h"

namespace ksub {

using mygraph::node_id;
using mygraph::edge_id;

// ======================================================
// MAX-K-CUT  (biên dịch với -DKFUNC_MAXKCUT)
//
// f(x) = sum_{(u,v) : x[u]!=0, x[v]!=0, x[u]!=x[v]} w_{uv}[0]
// Giả sử mỗi cạnh chỉ xuất hiện đúng 1 lần trong g.edges.
// ======================================================

double kfunc_evaluate(const mygraph::tinyGraph &g,
                      const Assignment &x)
{
    const std::size_t m = g.m;
    double value = 0.0;

    for (edge_id eid = 0; eid < m; ++eid) {
        const auto &E = g.edges[eid];
        node_id u = E.u;
        node_id v = E.v;

        Label lu = (u < x.size() ? x[u] : 0);
        Label lv = (v < x.size() ? x[v] : 0);

        if (lu == 0 || lv == 0) continue;
        if (lu == lv) continue;

        // dùng trọng số topic 0
        value += E.weights[0];
    }
    return value;
}

double kfunc_marginal(const mygraph::tinyGraph &g,
                      node_id u,
                      Label new_label,
                      const Assignment &x)
{
    const Label old_label = (u < x.size() ? x[u] : 0);
    if (old_label == new_label) return 0.0;

    double delta = 0.0;

    // Dùng danh sách cạnh kề u để chỉ duyệt deg(u) cạnh
    if (u >= g.incident.size()) {
        // graph không có index, fallback: duyệt toàn bộ cạnh
        for (edge_id eid = 0; eid < g.m; ++eid) {
            const auto &E = g.edges[eid];
            node_id a = E.u;
            node_id b = E.v;

            if (a != u && b != u) continue;
            node_id v = (a == u ? b : a);
            Label lv = (v < x.size() ? x[v] : 0);
            double w = E.weights[0];

            double before = 0.0;
            if (old_label != 0 && lv != 0 && old_label != lv) {
                before = w;
            }

            double after = 0.0;
            if (new_label != 0 && lv != 0 && new_label != lv) {
                after = w;
            }

            delta += (after - before);
        }
        return delta;
    }

    // Trường hợp chuẩn: dùng incident[u]
    for (edge_id eid : g.incident[u]) {
        const auto &E = g.edges[eid];

        node_id v;
        if (g.undirected) {
            v = (E.u == u ? E.v : E.u);
        } else {
            // với đồ thị có hướng, incident[u] hiện đang chứa các cạnh đi ra từ u
            v = (E.u == u ? E.v : E.u);
        }

        Label lv = (v < x.size() ? x[v] : 0);
        double w = E.weights[0];

        double before = 0.0;
        if (old_label != 0 && lv != 0 && old_label != lv) {
            before = w;
        }

        double after = 0.0;
        if (new_label != 0 && lv != 0 && new_label != lv) {
            after = w;
        }

        delta += (after - before);
    }

    return delta;
}

// Overload: nhận f_x để tránh phải tính lại f(x) trong một số thuật toán.
// Với Max-K-Cut, marginal đã tính trực tiếp theo local edges, nên f_x không cần dùng.
double kfunc_marginal(const mygraph::tinyGraph &g,
                      node_id u,
                      Label new_label,
                      const Assignment &x,
                      double f_x)
{
    return kfunc_marginal(g, u, new_label, x);
}

} // namespace ksub

#endif // KSUB_OBJECTVALUE_MAXKCUT_H

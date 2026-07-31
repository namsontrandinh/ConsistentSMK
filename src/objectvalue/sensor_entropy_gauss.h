// src/objectvalue/sensor_entropy_gauss.h
#ifndef KSUB_OBJECTVALUE_SENSOR_ENTROPY_GAUSS_H
#define KSUB_OBJECTVALUE_SENSOR_ENTROPY_GAUSS_H

#include "objectvalue_common.h"

namespace ksub {

using mygraph::node_id;
using mygraph::edge_id;

// ======================================================
// SENSOR ENTROPY (Gaussian log-det) (compile with -DKFUNC_SENSOR_ENTROPY_GAUSS)
//
// f(x) = sum_{i=1..K} 0.5 * ( |S_i| * log(2*pi*e) + log det( Sigma^(i-1)[S_i,S_i] ) )
//   S_i = {u : x[u] = i}.
//
// Encoding trong tinyGraph (tạo bởi preIntelLab.py):
// - g.undirected == true
// - complete graph, m = n*(n-1)/2
// - edges lưu theo thứ tự: for u=0..n-1, for v=u+1..n-1 push (u,v)
// - Sigma^(t) (t=0..K-1):
//     diag: g.nodes[u].weights[t] = Sigma^(t)[u,u]
//     off : g.edges[eid(u,v)].weights[t] = Sigma^(t)[u,v], u<v
// ======================================================

static inline std::size_t complete_edge_id(std::size_t n, node_id a, node_id b)
{
    if (a == b) throw std::invalid_argument("complete_edge_id: a==b");
    node_id u = std::min(a, b);
    node_id v = std::max(a, b);

    // base = sum_{r=0}^{u-1} (n-1-r) = u*(n-1) - u*(u-1)/2
    std::size_t uu = static_cast<std::size_t>(u);
    std::size_t vv = static_cast<std::size_t>(v);
    std::size_t base = uu * (n - 1) - (uu * (uu - 1)) / 2;
    std::size_t off  = vv - uu - 1;
    return base + off;
}

// In-place Cholesky on SPD matrix A (row-major s*s).
// Stores L in lower triangle. Returns logdet(A) = 2 * sum log(L_ii).
static inline double logdet_spd_cholesky(std::vector<double> &A, std::size_t s)
{
    double sum_log_diag = 0.0;

    for (std::size_t i = 0; i < s; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double val = A[i*s + j];
            for (std::size_t k = 0; k < j; ++k) {
                val -= A[i*s + k] * A[j*s + k];
            }

            if (i == j) {
                if (val <= 0.0 || !std::isfinite(val)) {
                    throw std::runtime_error("Cholesky failed: Sigma is not SPD (increase ridge/spd_eps in preprocessor).");
                }
                double lij = std::sqrt(val);
                A[i*s + j] = lij;
                sum_log_diag += std::log(lij);
            } else {
                double ljj = A[j*s + j];
                A[i*s + j] = val / ljj;
            }
        }

        // optional: clear upper triangle
        for (std::size_t j = i + 1; j < s; ++j) A[i*s + j] = 0.0;
    }

    return 2.0 * sum_log_diag;
}

static inline double entropy_gauss_set_topic(const mygraph::tinyGraph &g,
                                             const std::vector<node_id> &S,
                                             std::size_t topic)
{
    const std::size_t s = S.size();
    if (s == 0) return 0.0;

    const double pi = std::acos(-1.0);
    const double LOG_2PIE = std::log(2.0 * pi) + 1.0; // log(2*pi*e)

    const std::size_t n = g.n;
    const std::size_t expected_m = n * (n - 1) / 2;

    if (!g.undirected || g.m != expected_m) {
        throw std::runtime_error("KFUNC_SENSOR_ENTROPY_GAUSS requires undirected complete graph produced by preIntelLab.py.");
    }
    if (topic >= g.K) {
        throw std::runtime_error("KFUNC_SENSOR_ENTROPY_GAUSS: topic out of range.");
    }

    // Build Sigma[S,S] into dense matrix A (row-major)
    std::vector<double> A(s * s, 0.0);

    for (std::size_t i = 0; i < s; ++i) {
        node_id ui = S[i];

        double di = g.nodes[ui].weights[topic];
        if (!(di > 0.0) || !std::isfinite(di)) di = 1e-12;
        A[i*s + i] = di;

        for (std::size_t j = i + 1; j < s; ++j) {
            node_id uj = S[j];
            std::size_t eid = complete_edge_id(n, ui, uj);
            if (eid >= g.edges.size()) {
                throw std::runtime_error("Graph edge order mismatch: eid out of range (not a complete graph in required order).");
            }
            double c = g.edges[eid].weights[topic];
            if (!std::isfinite(c)) c = 0.0;
            A[i*s + j] = c;
            A[j*s + i] = c;
        }
    }

    double logdet = 0.0;
    if (s == 1) {
        logdet = std::log(A[0]);
    } else {
        logdet = logdet_spd_cholesky(A, s);
    }

    return 0.5 * (static_cast<double>(s) * LOG_2PIE + logdet);
}

double kfunc_evaluate(const mygraph::tinyGraph &g,
                      const Assignment &x)
{
    const std::size_t n = g.n;
    const std::size_t K = g.K;
    if (n == 0 || K == 0) return 0.0;

    double value = 0.0;

    // separable theo nhãn: topic index = (label-1)
    for (std::size_t lab = 1; lab <= K; ++lab) {
        std::vector<node_id> S;
        S.reserve(64);

        for (node_id u = 0; u < static_cast<node_id>(n); ++u) {
            Label lu = (u < x.size() ? x[u] : static_cast<Label>(0));
            if (lu == static_cast<Label>(lab)) S.push_back(u);
        }

        value += entropy_gauss_set_topic(g, S, lab - 1);
    }

    return value;
}

static inline void gather_label_set(const mygraph::tinyGraph &g,
                                    const Assignment &x,
                                    Label lab,
                                    node_id skip,
                                    std::vector<node_id> &out)
{
    out.clear();
    if (lab == 0) return;
    const std::size_t n = g.n;
    for (node_id v = 0; v < static_cast<node_id>(n); ++v) {
        if (v == skip) continue;
        Label lv = (v < x.size() ? x[v] : static_cast<Label>(0));
        if (lv == lab) out.push_back(v);
    }
}

double kfunc_marginal(const mygraph::tinyGraph &g,
                      node_id u,
                      Label new_label,
                      const Assignment &x)
{
    const std::size_t K = g.K;

    const Label old_label = (u < x.size() ? x[u] : static_cast<Label>(0));
    if (old_label == new_label) return 0.0;

    if (new_label > static_cast<Label>(K) || old_label > static_cast<Label>(K)) {
        throw std::runtime_error("KFUNC_SENSOR_ENTROPY_GAUSS: label out of range.");
    }

    double delta = 0.0;
    std::vector<node_id> S_before, S_after;
    S_before.reserve(64);
    S_after.reserve(64);

    // (1) remove from old_label, if any
    if (old_label != 0) {
        const std::size_t topic_old = static_cast<std::size_t>(old_label - 1);
        gather_label_set(g, x, old_label, u, S_before); // S_old without u
        // before includes u
        S_after = S_before; // after is without u
        S_before.push_back(u);

        const double before = entropy_gauss_set_topic(g, S_before, topic_old);
        const double after  = entropy_gauss_set_topic(g, S_after,  topic_old);
        delta += (after - before);
    }

    // (2) add to new_label, if any
    if (new_label != 0) {
        const std::size_t topic_new = static_cast<std::size_t>(new_label - 1);
        gather_label_set(g, x, new_label, u, S_before); // S_new without u
        S_after = S_before;
        S_after.push_back(u);

        const double before = entropy_gauss_set_topic(g, S_before, topic_new);
        const double after  = entropy_gauss_set_topic(g, S_after,  topic_new);
        delta += (after - before);
    }

    return delta;
}

// Overload: với objective này, delta đã tính đúng cục bộ theo (old_label,new_label),
// nên không cần dùng f_x (vẫn giữ signature để đồng bộ interface).
double kfunc_marginal(const mygraph::tinyGraph &g,
                      node_id u,
                      Label new_label,
                      const Assignment &x,
                      double /*f_x*/)
{
    return kfunc_marginal(g, u, new_label, x);
}

} // namespace ksub

#endif // KSUB_OBJECTVALUE_SENSOR_ENTROPY_GAUSS_H

// mygraph.h
// Supports THREE binary formats:
//
// LEGACY SCALAR FORMAT:
//   [u64 n][u64 m][u8 undirected]
//   nodes: one double weight + one double alpha per node
//   edges: [u32 u][u32 v] + one double weight per edge
//
// OLD VECTOR FORMAT:
//   [u64 n][u64 m][u64 K][u8 undirected]
//   nodes: K doubles weights + K doubles alpha for each node
//   edges: [u32 u][u32 v] + K doubles weights for each edge
//
// NEW FORMAT:
//   [u64 n][u64 m][u64 K][u8 undirected]
//   nodes: K doubles weights + K doubles alpha for each node
//   part_id: n x [u32]
//   edges: [u32 u][u32 v] + K doubles weights for each edge
//
// write_binary() always writes the NEW FORMAT.
// read_binary() automatically detects and reads OLD or NEW FORMAT.

#ifndef MYGRAPH_H
#define MYGRAPH_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace mygraph {

using node_id = std::uint32_t;
using edge_id = std::size_t;

struct Edge {
    node_id u;
    node_id v;
    std::vector<double> weights; // size = K
};

struct Node {
    std::vector<double> weights; // size = K
    std::vector<double> alpha;   // size = K
};

struct tinyGraph {
    std::size_t n = 0;
    std::size_t m = 0;
    std::size_t K = 0;
    bool undirected = true;

    std::vector<Node> nodes;
    std::vector<Edge> edges;

    // Partition id for matroid:
    // part_id[u] in {0,...,p-1}.
    // For OLD binary files, read_binary() initializes all values to 0.
    std::vector<std::uint32_t> part_id;

    // indexes
    std::vector<std::vector<edge_id>> incident;
    std::vector<std::vector<edge_id>> incoming;
    std::vector<std::vector<node_id>> neighbors;

    inline void clear();

    inline void init(std::size_t n_nodes,
                     std::size_t k_topics,
                     bool undirected_flag);

    inline void build_index();

    inline edge_id add_edge(node_id u,
                            node_id v,
                            const std::vector<double>& w);

    inline double edge_weight(edge_id e,
                              std::size_t topic) const;

    template <class F>
    inline void for_each_edge(F&& f) const {
        for (edge_id e = 0; e < edges.size(); ++e) {
            const auto& E = edges[e];
            f(E.u, E.v, e);
        }
    }

    template <class F>
    inline void for_each_undirected_edge(F&& f) const {
        if (!undirected) {
            for_each_edge(std::forward<F>(f));
            return;
        }
        for_each_edge(std::forward<F>(f));
    }

    inline bool write_binary(const std::string& path) const;

    // Automatically detects OLD or NEW binary format.
    inline bool read_binary(const std::string& path);
};

// ======================================================
// Preprocess (txt -> bin) with optional partition file.
//
// Edge-list txt:
//   each line: u v [w]
//
// Optional partition txt:
//   each line: orig_id part
//
// Missing partition or missing node -> part_id = 0.
//
// Output is always NEW FORMAT.
// ======================================================
inline bool preprocess_edge_list_to_binary(
    const std::string& txt_path,
    const std::string& bin_path,
    std::size_t k_topics,
    bool undirected = true,
    bool randomize_node = false,
    unsigned seed = 42,
    const std::string& part_txt_path = "");

// =========================
// inline implementations
// =========================

inline void tinyGraph::clear() {
    n = 0;
    m = 0;
    K = 0;
    undirected = true;

    nodes.clear();
    edges.clear();
    part_id.clear();

    incident.clear();
    incoming.clear();
    neighbors.clear();
}

inline void tinyGraph::init(std::size_t n_nodes,
                            std::size_t k_topics,
                            bool undirected_flag)
{
    clear();

    n = n_nodes;
    K = k_topics;
    undirected = undirected_flag;

    nodes.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        nodes[i].weights.assign(K, 1.0);
        nodes[i].alpha.assign(K, 1.0);
    }

    part_id.assign(n, 0);

    incident.assign(n, std::vector<edge_id>());
    incoming.assign(n, std::vector<edge_id>());
    neighbors.assign(n, std::vector<node_id>());
}

inline void tinyGraph::build_index() {
    incident.assign(n, std::vector<edge_id>());
    incoming.assign(n, std::vector<edge_id>());
    neighbors.assign(n, std::vector<node_id>());

    for (edge_id e = 0; e < edges.size(); ++e) {
        const auto& E = edges[e];

        const node_id u = E.u;
        const node_id v = E.v;

        if (u >= n || v >= n) {
            continue;
        }

        if (undirected) {
            incident[u].push_back(e);
            incident[v].push_back(e);

            incoming[u].push_back(e);
            incoming[v].push_back(e);

            neighbors[u].push_back(v);
            neighbors[v].push_back(u);
        } else {
            // directed u -> v
            incident[u].push_back(e);
            incoming[v].push_back(e);
            neighbors[u].push_back(v);
        }
    }
}

inline edge_id tinyGraph::add_edge(
    node_id u,
    node_id v,
    const std::vector<double>& w)
{
    if (w.size() != K) {
        throw std::invalid_argument(
            "add_edge: weight vector size != K");
    }

    if (u >= n || v >= n) {
        throw std::out_of_range(
            "add_edge: node id out of range");
    }

    Edge e;
    e.u = u;
    e.v = v;
    e.weights = w;

    edges.push_back(std::move(e));

    const edge_id id =
        edges.size() - 1;

    m = edges.size();

    if (incident.size() != n ||
        incoming.size() != n ||
        neighbors.size() != n)
    {
        incident.assign(n, std::vector<edge_id>());
        incoming.assign(n, std::vector<edge_id>());
        neighbors.assign(n, std::vector<node_id>());
        build_index();
        return id;
    }

    if (undirected) {
        incident[u].push_back(id);
        incident[v].push_back(id);

        incoming[u].push_back(id);
        incoming[v].push_back(id);

        neighbors[u].push_back(v);
        neighbors[v].push_back(u);
    } else {
        incident[u].push_back(id);
        incoming[v].push_back(id);
        neighbors[u].push_back(v);
    }

    return id;
}

inline double tinyGraph::edge_weight(
    edge_id e,
    std::size_t topic) const
{
    if (e >= edges.size()) {
        throw std::out_of_range(
            "edge_weight: edge id out of range");
    }

    if (topic >= K) {
        throw std::out_of_range(
            "edge_weight: topic index out of range");
    }

    return edges[e].weights[topic];
}

// ======================================================
// NEW Binary format writer.
//
// [u64 n][u64 m][u64 K][u8 undirected]
//
// nodes:
//   for i=0..n-1:
//     K doubles weights
//     K doubles alpha
//
// part_id:
//   for i=0..n-1:
//     [u32 part_id[i]]
//
// edges:
//   for e=0..m-1:
//     [u32 u][u32 v]
//     K doubles weights
// ======================================================
inline bool tinyGraph::write_binary(
    const std::string& path) const
{
    if (part_id.size() != n) {
        throw std::invalid_argument(
            "write_binary: part_id.size() != n");
    }

    if (nodes.size() != n) {
        throw std::invalid_argument(
            "write_binary: nodes.size() != n");
    }

    if (edges.size() != m) {
        throw std::invalid_argument(
            "write_binary: edges.size() != m");
    }

    std::ofstream out(
        path,
        std::ios::binary);

    if (!out) {
        return false;
    }

    const std::uint64_t nn =
        static_cast<std::uint64_t>(n);

    const std::uint64_t mm =
        static_cast<std::uint64_t>(m);

    const std::uint64_t KK =
        static_cast<std::uint64_t>(K);

    const std::uint8_t und =
        undirected ? 1 : 0;

    out.write(
        reinterpret_cast<const char*>(&nn),
        sizeof(nn));

    out.write(
        reinterpret_cast<const char*>(&mm),
        sizeof(mm));

    out.write(
        reinterpret_cast<const char*>(&KK),
        sizeof(KK));

    out.write(
        reinterpret_cast<const char*>(&und),
        sizeof(und));

    // nodes
    for (std::size_t i = 0; i < n; ++i) {
        if (nodes[i].weights.size() != K ||
            nodes[i].alpha.size() != K)
        {
            return false;
        }

        for (std::size_t t = 0; t < K; ++t) {
            const double w =
                nodes[i].weights[t];

            out.write(
                reinterpret_cast<const char*>(&w),
                sizeof(double));
        }

        for (std::size_t t = 0; t < K; ++t) {
            const double a =
                nodes[i].alpha[t];

            out.write(
                reinterpret_cast<const char*>(&a),
                sizeof(double));
        }
    }

    // NEW FORMAT: serialized part_id
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint32_t pid =
            part_id[i];

        out.write(
            reinterpret_cast<const char*>(&pid),
            sizeof(pid));
    }

    // edges
    for (const auto& E : edges) {
        if (E.weights.size() != K) {
            return false;
        }

        out.write(
            reinterpret_cast<const char*>(&E.u),
            sizeof(node_id));

        out.write(
            reinterpret_cast<const char*>(&E.v),
            sizeof(node_id));

        for (std::size_t t = 0; t < K; ++t) {
            const double w =
                E.weights[t];

            out.write(
                reinterpret_cast<const char*>(&w),
                sizeof(double));
        }
    }

    return static_cast<bool>(out);
}

// ======================================================
// Binary reader supporting BOTH OLD and NEW FORMAT.
//
// OLD:
//   header + nodes + edges
//
// NEW:
//   header + nodes + part_id + edges
//
// Detection is based on exact file size computed from
// n, m and K stored in the header.
//
// OLD files receive:
//   part_id.assign(n, 0)
//
// This keeps old EDL/IC data usable while new files retain
// serialized partition information.
// ======================================================
inline bool tinyGraph::read_binary(
    const std::string& path)
{
    clear();

    std::ifstream in(path, std::ios::binary);
    if (!in) {
        return false;
    }

    // Get actual file size.
    in.seekg(0, std::ios::end);
    const std::streamoff end_pos = in.tellg();
    if (end_pos < 0) {
        clear();
        return false;
    }

    const std::uint64_t file_size =
        static_cast<std::uint64_t>(end_pos);

    in.seekg(0, std::ios::beg);

    // --------------------------------------------------
    // First read the common prefix [u64 n][u64 m].
    // --------------------------------------------------
    std::uint64_t nn = 0;
    std::uint64_t mm = 0;

    in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
    in.read(reinterpret_cast<char*>(&mm), sizeof(mm));

    if (!in || nn == 0) {
        clear();
        return false;
    }

    // Checked arithmetic helpers.
    const std::uint64_t UMAX =
        std::numeric_limits<std::uint64_t>::max();

    auto checked_mul =
        [UMAX](std::uint64_t a,
               std::uint64_t b,
               std::uint64_t& out) -> bool
    {
        if (a != 0 && b > UMAX / a) {
            return false;
        }
        out = a * b;
        return true;
    };

    auto checked_add =
        [UMAX](std::uint64_t a,
               std::uint64_t b,
               std::uint64_t& out) -> bool
    {
        if (b > UMAX - a) {
            return false;
        }
        out = a + b;
        return true;
    };

    // ==================================================
    // FORMAT A: LEGACY SCALAR
    //
    // [u64 n][u64 m][u8 undirected]
    // nodes:
    //   n * ([double weight][double alpha])
    // edges:
    //   m * ([u32 u][u32 v][double weight])
    //
    // This exactly matches the old fb.bin/mail.bin layout.
    // Convert it in memory to the current vector representation
    // by setting K=1.
    // ==================================================
    std::uint64_t legacy_node_bytes = 0;
    std::uint64_t legacy_edge_bytes = 0;
    std::uint64_t legacy_expected_size = 0;
    std::uint64_t legacy_tmp = 0;

    const std::uint64_t legacy_header_bytes =
        2ULL * static_cast<std::uint64_t>(sizeof(std::uint64_t)) +
        static_cast<std::uint64_t>(sizeof(std::uint8_t));

    if (checked_mul(
            nn,
            2ULL * static_cast<std::uint64_t>(sizeof(double)),
            legacy_node_bytes) &&
        checked_mul(
            mm,
            2ULL * static_cast<std::uint64_t>(sizeof(node_id)) +
                static_cast<std::uint64_t>(sizeof(double)),
            legacy_edge_bytes) &&
        checked_add(
            legacy_header_bytes,
            legacy_node_bytes,
            legacy_tmp) &&
        checked_add(
            legacy_tmp,
            legacy_edge_bytes,
            legacy_expected_size) &&
        file_size == legacy_expected_size)
    {
        // Rewind and parse the legacy scalar format.
        in.clear();
        in.seekg(0, std::ios::beg);

        std::uint8_t und = 0;

        in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
        in.read(reinterpret_cast<char*>(&mm), sizeof(mm));
        in.read(reinterpret_cast<char*>(&und), sizeof(und));

        if (!in) {
            clear();
            return false;
        }

        if (nn >
                static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max()) ||
            mm >
                static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max()))
        {
            clear();
            return false;
        }

        n = static_cast<std::size_t>(nn);
        m = static_cast<std::size_t>(mm);
        K = 1;
        undirected = (und != 0);

        nodes.resize(n);

        for (std::size_t i = 0; i < n; ++i) {
            double w = 0.0;
            double a = 0.0;

            in.read(reinterpret_cast<char*>(&w), sizeof(double));
            in.read(reinterpret_cast<char*>(&a), sizeof(double));

            if (!in) {
                clear();
                return false;
            }

            nodes[i].weights.assign(1, w);
            nodes[i].alpha.assign(1, a);
        }

        // Legacy files do not contain partition information.
        part_id.assign(n, 0);

        edges.resize(m);

        for (std::size_t e = 0; e < m; ++e) {
            node_id u = 0;
            node_id v = 0;
            double w = 0.0;

            in.read(reinterpret_cast<char*>(&u), sizeof(node_id));
            in.read(reinterpret_cast<char*>(&v), sizeof(node_id));
            in.read(reinterpret_cast<char*>(&w), sizeof(double));

            if (!in) {
                clear();
                return false;
            }

            if (u >= n || v >= n) {
                clear();
                return false;
            }

            edges[e].u = u;
            edges[e].v = v;
            edges[e].weights.assign(1, w);
        }

        build_index();
        return true;
    }

    // ==================================================
    // FORMAT B/C: VECTOR FORMAT
    //
    // Common header:
    // [u64 n][u64 m][u64 K][u8 undirected]
    //
    // OLD VECTOR:
    //   header + nodes + edges
    //
    // NEW VECTOR:
    //   header + nodes + part_id + edges
    // ==================================================
    in.clear();
    in.seekg(0, std::ios::beg);

    std::uint64_t KK = 0;
    std::uint8_t und = 0;

    in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
    in.read(reinterpret_cast<char*>(&mm), sizeof(mm));
    in.read(reinterpret_cast<char*>(&KK), sizeof(KK));
    in.read(reinterpret_cast<char*>(&und), sizeof(und));

    if (!in || nn == 0 || KK == 0) {
        clear();
        return false;
    }

    const std::uint64_t header_bytes =
        3ULL * static_cast<std::uint64_t>(sizeof(std::uint64_t)) +
        static_cast<std::uint64_t>(sizeof(std::uint8_t));

    // nodes = n * K * (weight + alpha) doubles
    std::uint64_t node_bytes = 0;
    std::uint64_t tmp1 = 0;
    std::uint64_t tmp2 = 0;

    if (!checked_mul(nn, KK, tmp1) ||
        !checked_mul(tmp1, 2ULL, tmp2) ||
        !checked_mul(
            tmp2,
            static_cast<std::uint64_t>(sizeof(double)),
            node_bytes))
    {
        clear();
        return false;
    }

    // each edge = u32 + v32 + K doubles
    std::uint64_t edge_weight_bytes = 0;

    if (!checked_mul(
            KK,
            static_cast<std::uint64_t>(sizeof(double)),
            edge_weight_bytes))
    {
        clear();
        return false;
    }

    std::uint64_t one_edge_bytes = 0;

    if (!checked_add(
            2ULL * static_cast<std::uint64_t>(sizeof(node_id)),
            edge_weight_bytes,
            one_edge_bytes))
    {
        clear();
        return false;
    }

    std::uint64_t edge_bytes = 0;

    if (!checked_mul(mm, one_edge_bytes, edge_bytes)) {
        clear();
        return false;
    }

    std::uint64_t header_plus_nodes = 0;
    std::uint64_t old_expected_size = 0;

    if (!checked_add(
            header_bytes,
            node_bytes,
            header_plus_nodes) ||
        !checked_add(
            header_plus_nodes,
            edge_bytes,
            old_expected_size))
    {
        clear();
        return false;
    }

    std::uint64_t part_bytes = 0;

    if (!checked_mul(
            nn,
            static_cast<std::uint64_t>(sizeof(std::uint32_t)),
            part_bytes))
    {
        clear();
        return false;
    }

    std::uint64_t new_expected_size = 0;

    if (!checked_add(
            old_expected_size,
            part_bytes,
            new_expected_size))
    {
        clear();
        return false;
    }

    bool has_part_id = false;

    if (file_size == new_expected_size) {
        has_part_id = true;
    }
    else if (file_size == old_expected_size) {
        has_part_id = false;
    }
    else {
        clear();
        return false;
    }

    if (nn >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::size_t>::max()) ||
        mm >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::size_t>::max()) ||
        KK >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::size_t>::max()))
    {
        clear();
        return false;
    }

    n = static_cast<std::size_t>(nn);
    m = static_cast<std::size_t>(mm);
    K = static_cast<std::size_t>(KK);
    undirected = (und != 0);

    nodes.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        nodes[i].weights.resize(K);
        nodes[i].alpha.resize(K);

        for (std::size_t t = 0; t < K; ++t) {
            double w = 0.0;
            in.read(reinterpret_cast<char*>(&w), sizeof(double));

            if (!in) {
                clear();
                return false;
            }

            nodes[i].weights[t] = w;
        }

        for (std::size_t t = 0; t < K; ++t) {
            double a = 0.0;
            in.read(reinterpret_cast<char*>(&a), sizeof(double));

            if (!in) {
                clear();
                return false;
            }

            nodes[i].alpha[t] = a;
        }
    }

    part_id.assign(n, 0);

    if (has_part_id) {
        for (std::size_t i = 0; i < n; ++i) {
            std::uint32_t pid = 0;

            in.read(reinterpret_cast<char*>(&pid), sizeof(pid));

            if (!in) {
                clear();
                return false;
            }

            part_id[i] = pid;
        }
    }

    edges.resize(m);

    for (std::size_t e = 0; e < m; ++e) {
        node_id u = 0;
        node_id v = 0;

        in.read(reinterpret_cast<char*>(&u), sizeof(node_id));
        in.read(reinterpret_cast<char*>(&v), sizeof(node_id));

        if (!in) {
            clear();
            return false;
        }

        if (u >= n || v >= n) {
            clear();
            return false;
        }

        edges[e].u = u;
        edges[e].v = v;
        edges[e].weights.resize(K);

        for (std::size_t t = 0; t < K; ++t) {
            double w = 0.0;

            in.read(reinterpret_cast<char*>(&w), sizeof(double));

            if (!in) {
                clear();
                return false;
            }

            edges[e].weights[t] = w;
        }
    }

    build_index();
    return true;
}

// ------------------------------------------------------
// helper: read orig_id -> part mapping
// ------------------------------------------------------
inline bool load_part_map_txt_(
    const std::string& part_txt_path,
    std::unordered_map<
        std::uint64_t,
        std::uint32_t>& partmap)
{
    std::ifstream pin(
        part_txt_path);

    if (!pin) {
        return false;
    }

    std::string line;

    while (std::getline(
        pin,
        line))
    {
        if (line.empty() ||
            line[0] == '#')
        {
            continue;
        }

        std::uint64_t u0 = 0;
        std::uint64_t p0 = 0;

        std::stringstream ss(
            line);

        ss >> u0 >> p0;

        if (!ss) {
            continue;
        }

        partmap[u0] =
            static_cast<std::uint32_t>(
                p0);
    }

    return true;
}

inline bool preprocess_edge_list_to_binary(
    const std::string& txt_path,
    const std::string& bin_path,
    std::size_t k_topics,
    bool undirected,
    bool randomize_node,
    unsigned seed,
    const std::string& part_txt_path)
{
    if (k_topics == 0) {
        return false;
    }

    std::ifstream in(
        txt_path);

    if (!in) {
        return false;
    }

    // optional partition map on original ids
    std::unordered_map<
        std::uint64_t,
        std::uint32_t> partmap;

    if (!part_txt_path.empty()) {
        if (!load_part_map_txt_(
                part_txt_path,
                partmap))
        {
            return false;
        }
    }

    std::unordered_map<
        std::uint64_t,
        node_id> idmap;

    idmap.reserve(
        1 << 16);

    struct RawEdge {
        std::uint64_t u_orig;
        std::uint64_t v_orig;
        double w;
    };

    std::vector<RawEdge> rawEdges;
    rawEdges.reserve(
        1 << 20);

    std::string line;

    while (std::getline(
        in,
        line))
    {
        if (line.empty() ||
            line[0] == '#')
        {
            continue;
        }

        std::uint64_t u0 = 0;
        std::uint64_t v0 = 0;
        double w = 1.0;

        std::stringstream ss(
            line);

        ss >> u0 >> v0;

        if (!ss) {
            continue;
        }

        if (ss >> w) {
            if (w <= 0.0) {
                w = 1e-6;
            }
        }

        rawEdges.push_back(
            {u0, v0, w});
    }

    in.close();

    node_id next_id = 0;

    std::vector<std::uint32_t>
        part_new;

    part_new.reserve(
        1 << 16);

    auto get_new_id =
        [&](std::uint64_t orig)
            -> node_id
    {
        const auto it =
            idmap.find(orig);

        if (it != idmap.end()) {
            return it->second;
        }

        const node_id nid =
            next_id++;

        idmap.emplace(
            orig,
            nid);

        std::uint32_t pid = 0;

        if (!partmap.empty()) {
            const auto itp =
                partmap.find(orig);

            if (itp != partmap.end()) {
                pid = itp->second;
            }
        }

        part_new.push_back(pid);

        return nid;
    };

    std::vector<Edge> edges;
    edges.reserve(
        rawEdges.size());

    for (const auto& re : rawEdges) {
        const node_id u =
            get_new_id(
                re.u_orig);

        const node_id v =
            get_new_id(
                re.v_orig);

        Edge E;
        E.u = u;
        E.v = v;
        E.weights.assign(
            k_topics,
            re.w);

        edges.push_back(
            std::move(E));
    }

    tinyGraph g;

    g.n =
        static_cast<std::size_t>(
            next_id);

    g.m =
        edges.size();

    g.K =
        k_topics;

    g.undirected =
        undirected;

    // nodes
    g.nodes.resize(
        g.n);

    for (std::size_t i = 0;
         i < g.n;
         ++i)
    {
        g.nodes[i].weights.assign(
            k_topics,
            1.0);

        g.nodes[i].alpha.assign(
            k_topics,
            1.0);
    }

    if (randomize_node) {
        std::mt19937 gen(
            seed);

        std::uniform_real_distribution<double>
            dist(
                1e-6,
                1.0);

        for (std::size_t i = 0;
             i < g.n;
             ++i)
        {
            for (std::size_t t = 0;
                 t < k_topics;
                 ++t)
            {
                g.nodes[i].weights[t] =
                    dist(gen);

                g.nodes[i].alpha[t] =
                    dist(gen);
            }
        }
    }

    // part_id
    g.part_id.assign(
        g.n,
        0);

    if (!part_new.empty()) {
        if (part_new.size() !=
            g.n)
        {
            return false;
        }

        for (std::size_t i = 0;
             i < g.n;
             ++i)
        {
            g.part_id[i] =
                part_new[i];
        }
    }

    // edges
    g.edges =
        std::move(edges);

    return g.write_binary(
        bin_path);
}

} // namespace mygraph

#endif // MYGRAPH_H

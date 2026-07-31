// src/data/preprocess_ic.cpp  (SUBMODULAR-ONLY for IC, ALWAYS normalize incoming)
// - Output ALWAYS directed (undirected=false)
// - If input_undirected=1: write both (u,v) and (v,u)
// - Dedup THÔ: keep-first, drop-later for duplicate (u,v)
// - Clamp negative weights to 0
// - ALWAYS normalize incoming weights: sum_in[v]>0 => w(u->v)/=sum_in[v]
// - Renumber node_id by first appearance in edges.txt (deterministic given file order)

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <random>
#include <cstdint>
#include <algorithm>

#include "mygraph.h"

// Pack (u,v) into uint64 key (u and v are node_id <= 2^32-1)
static inline std::uint64_t pack_uv(std::uint32_t u, std::uint32_t v) {
    return (static_cast<std::uint64_t>(u) << 32) | static_cast<std::uint64_t>(v);
}
static inline std::uint32_t unpack_u(std::uint64_t key) {
    return static_cast<std::uint32_t>(key >> 32);
}
static inline std::uint32_t unpack_v(std::uint64_t key) {
    return static_cast<std::uint32_t>(key & 0xFFFFFFFFULL);
}

static inline double clamp_nonneg(double x) {
    return (x >= 0.0 ? x : 0.0);
}
static inline double clamp01(double p) {
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return 1.0;
    return p;
}

static void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog
        << " edges.txt output.bin [input_undirected(0/1)=1] [randomize_node(0/1)=0] [seed=42]\n\n"
        << "edges.txt format:\n"
        << "  - Moi dong: u v [w]\n"
        << "    + u, v: so nguyen khong am (ID goc), co the khong lien tiep\n"
        << "    + w   : optional, mac dinh 1.0\n"
        << "    + w am: clamp ve 0\n"
        << "  - Dong rong hoac bat dau bang '#' se bi bo qua.\n\n"
        << "Renumber rule (node id):\n"
        << "  - ID goc trong edges.txt duoc map sang node_id 0..n-1 theo THU TU XUAT HIEN LAN DAU.\n\n"
        << "Output:\n"
        << "  - tinyGraph luu ALWAYS directed (undirected=false)\n"
        << "  - Neu input_undirected=1: ghi ca (u,v) va (v,u)\n"
        << "  - Dedup (u,v) THO: neu gap lan 2 -> BO QUA\n"
        << "  - ALWAYS normalize incoming: voi moi v, neu sum_in[v]>0 thi w(u->v)/=sum_in[v] (tong incoming = 1)\n";
}

int main(int argc, char** argv) {
    using namespace std;
    using namespace mygraph;

    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    const string edges_txt  = argv[1];
    const string output_bin = argv[2];

    bool input_undirected = true;
    bool randomize_node = false;
    unsigned seed = 42;

    if (argc >= 4) input_undirected = (std::atoi(argv[3]) != 0);
    if (argc >= 5) randomize_node = (std::atoi(argv[4]) != 0);
    if (argc >= 6) seed = static_cast<unsigned>(std::strtoul(argv[5], nullptr, 10));

    cout << "Preprocess IC (submodular-only, ALWAYS normalize incoming):\n";
    cout << "  edges.txt         = " << edges_txt  << "\n";
    cout << "  output.bin        = " << output_bin << "\n";
    cout << "  input_undirected  = " << (input_undirected ? "true" : "false") << "\n";
    cout << "  output directed   = true (undirected=false)\n";
    cout << "  randomize_node    = " << (randomize_node ? "true" : "false") << "\n";
    cout << "  seed              = " << seed << "\n";

    // --------- 1) Read edges and build renumber map ----------
    ifstream fin(edges_txt);
    if (!fin) {
        cerr << "Error: cannot open edges file: " << edges_txt << "\n";
        return 1;
    }

    // Map original ids to compact [0..n-1] (by first appearance in edges)
    unordered_map<std::uint64_t, node_id> idmap;
    idmap.reserve(1 << 20);

    node_id next_id = 0;
    auto get_new_id = [&](std::uint64_t orig) -> node_id {
        auto it = idmap.find(orig);
        if (it != idmap.end()) return it->second;
        node_id nid = next_id++;
        idmap.emplace(orig, nid);
        return nid;
    };

    // Dedup THÔ (keep-first): key=(u,v) -> weight_raw (>=0)
    unordered_map<std::uint64_t, double> keep_first;
    keep_first.reserve(1 << 20);

    size_t edge_line_no = 0;
    size_t num_lines = 0;
    size_t num_kept_directed = 0;
    size_t num_skipped_directed = 0;

    auto try_insert_keep_first = [&](node_id a, node_id b, double w_raw) {
        const std::uint64_t key = pack_uv(static_cast<std::uint32_t>(a), static_cast<std::uint32_t>(b));
        auto it = keep_first.find(key);
        if (it != keep_first.end()) {
            ++num_skipped_directed;
            return;
        }
        keep_first.emplace(key, w_raw);
        ++num_kept_directed;
    };

    string line;
    while (std::getline(fin, line)) {
        ++edge_line_no;
        if (line.empty() || line[0] == '#') continue;

        ++num_lines;

        stringstream ss(line);
        std::uint64_t u0, v0;
        if (!(ss >> u0 >> v0)) continue;

        node_id u = get_new_id(u0);
        node_id v = get_new_id(v0);

        double w = 1.0;
        if (ss >> w) w = clamp_nonneg(w);
        else w = 1.0;

        // Output ALWAYS directed, dedup THÔ
        try_insert_keep_first(u, v, w);
        if (input_undirected) {
            try_insert_keep_first(v, u, w);
        }
    }
    fin.close();

    const size_t n = static_cast<size_t>(next_id);
    if (n == 0) {
        cerr << "Error: empty graph (no nodes found in edges.txt).\n";
        return 1;
    }

    cout << "  Renumbered nodes n  = " << n << "\n";
    cout << "  lines_read          = " << num_lines << "\n";
    cout << "  kept_directed       = " << num_kept_directed << "\n";
    cout << "  skipped_duplicates  = " << num_skipped_directed << "\n";

    // --------- 2) Materialize edges in deterministic order (sort by (u,v)) ----------
    std::vector<std::uint64_t> keys;
    keys.reserve(keep_first.size());
    for (const auto &kv : keep_first) keys.push_back(kv.first);

    std::sort(keys.begin(), keys.end(), [](std::uint64_t a, std::uint64_t b) {
        const std::uint32_t au = unpack_u(a), av = unpack_v(a);
        const std::uint32_t bu = unpack_u(b), bv = unpack_v(b);
        if (au != bu) return au < bu;
        return av < bv;
    });

    // --------- 3) Build directed tinyGraph (submodular-only format) ----------
    tinyGraph g;
    g.init(n, /*undirected_flag=*/false); // ALWAYS directed
    g.edges.clear();
    g.edges.reserve(keys.size());

    for (std::uint64_t key : keys) {
        const node_id u = static_cast<node_id>(unpack_u(key));
        const node_id v = static_cast<node_id>(unpack_v(key));

        auto it = keep_first.find(key);
        double w_raw = (it != keep_first.end() ? it->second : 0.0);
        if (w_raw < 0.0) w_raw = 0.0; // safety

        // tạm lưu w_raw, sẽ normalize sau
        g.edges.push_back(Edge{u, v, w_raw});
    }
    g.m = g.edges.size();

    // --------- 4) ALWAYS normalize incoming (scalar) ----------
    std::vector<double> sum_in(g.n, 0.0);

    for (edge_id eid = 0; eid < g.edges.size(); ++eid) {
        const auto &E = g.edges[eid];
        const node_id v = E.v;
        if (v >= g.n) continue;
        if (E.weight > 0.0) sum_in[static_cast<size_t>(v)] += E.weight;
    }

    for (edge_id eid = 0; eid < g.edges.size(); ++eid) {
        auto &E = g.edges[eid];
        const node_id v = E.v;
        if (v >= g.n) continue;

        const double s = sum_in[static_cast<size_t>(v)];
        if (s > 0.0) E.weight = E.weight / s;
        // safety: IC cần xác suất
        E.weight = clamp01(E.weight);
    }

    // --------- 5) Node weights/alpha (optional randomize) ----------
    if (randomize_node) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(1e-6, 1.0);
        for (size_t i = 0; i < g.n; ++i) {
            g.nodes[i].weight = dist(gen);
            g.nodes[i].alpha  = dist(gen);
        }
    }
    // else: defaults Node{weight=1, alpha=1}

    // --------- 6) Write binary ----------
    if (!g.write_binary(output_bin)) {
        cerr << "Error: write_binary failed: " << output_bin << "\n";
        return 1;
    }

    cout << "Done. Saved binary graph to: " << output_bin << "\n";
    cout << "Summary:\n";
    cout << "  n = " << g.n << ", m = " << g.m
         << ", undirected(stored) = " << (g.undirected ? "true" : "false") << "\n";
    cout << "  normalize_incoming = true\n";
    cout << "  dedup              = keep-first, dropped " << num_skipped_directed << " duplicates\n";
    return 0;
}

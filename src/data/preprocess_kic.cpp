// preprocess_kic.cpp  (Partition-aware, FORMAT A ONLY)
// - Output ALWAYS directed (undirected=false)
// - If input_undirected=1: write both (u,v) and (v,u)
// - Dedup THÔ: keep-first, drop-later for duplicate (u,v)
// - Clamp negative weights to 0
// - Normalize incoming weights per topic (for directed graph)
// - Renumber node_id by first appearance in edges.txt (deterministic given file order)
// - Compress part_id to 0..p-1 by ascending raw part_id (deterministic)

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

static void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog
        << " edges.txt parts.txt output.bin K [input_undirected(0/1)=1] [randomize_node(0/1)=0] [seed=42]\n\n"
        << "edges.txt format:\n"
        << "  - Moi dong: u v [w]    OR    u v w1 w2 ... wK\n"
        << "    + u, v: so nguyen khong am (ID goc), co the khong lien tiep\n"
        << "    + w   : neu chi 1 so -> replicate cho K topic\n"
        << "    + w1..wK: neu >=K so -> lay K so dau; neu <K -> pad 0\n"
        << "    + trong so am: clamp ve 0\n"
        << "  - Dong rong hoac bat dau bang '#' se bi bo qua.\n\n"
        << "parts.txt format (BAT BUOC - Format A):\n"
        << "  - Moi dong: orig_node_id  part_id\n"
        << "    + orig_node_id: ID goc (dung he nhu trong edges.txt)\n"
        << "    + part_id     : so nguyen >= 0 (co the bat ky; se duoc renumber ve 0..p-1)\n"
        << "  - Moi orig_node_id xuat hien TOI DA 1 lan (neu trung -> loi)\n"
        << "  - Bat buoc DUNG 2 cot; sai dinh dang -> loi\n\n"
        << "Renumber rule (node id):\n"
        << "  - Cac ID goc trong edges.txt duoc map sang node_id 0..n-1 theo THU TU XUAT HIEN LAN DAU.\n"
        << "  - parts.txt dung ID GOC nen KHONG bi anh huong boi renumber node.\n\n"
        << "Renumber rule (part id):\n"
        << "  - part_id doc tu parts.txt se duoc 'compress' ve 0,1,2,... theo THU TU TANG DAN\n"
        << "    cua part_id goc (deterministic).\n\n"
        << "Strict checks:\n"
        << "  - Moi node (xuat hien trong edges.txt) PHAI co part_id trong parts.txt; thieu -> loi.\n"
        << "  - parts.txt co node khong xuat hien trong edges.txt -> loi.\n\n"
        << "Output:\n"
        << "  - tinyGraph luu ALWAYS directed (undirected=false)\n"
        << "  - Neu input_undirected=1: ghi ca (u,v) va (v,u)\n"
        << "  - Dedup (u,v) THO: neu gap lan 2 -> BO QUA\n"
        << "  - Normalize incoming theo tung topic: sum_in[v,t]>0 => w(u->v,t)/=sum_in[v,t]\n";
}

int main(int argc, char** argv) {
    using namespace std;
    using namespace mygraph;

    if (argc < 5) {
        print_usage(argv[0]);
        return 1;
    }

    const string edges_txt  = argv[1];
    const string parts_txt  = argv[2];
    const string output_bin = argv[3];
    const size_t K = static_cast<size_t>(std::strtoull(argv[4], nullptr, 10));

    if (K == 0) {
        cerr << "Error: K must be > 0.\n";
        return 1;
    }

    bool input_undirected = true;
    bool randomize_node = false;
    unsigned seed = 42;

    if (argc >= 6) input_undirected = (std::atoi(argv[5]) != 0);
    if (argc >= 7) randomize_node = (std::atoi(argv[6]) != 0);
    if (argc >= 8) seed = static_cast<unsigned>(std::strtoul(argv[7], nullptr, 10));

    cout << "Preprocess KIC (partition-aware, FORMAT A ONLY):\n";
    cout << "  edges.txt         = " << edges_txt  << "\n";
    cout << "  parts.txt         = " << parts_txt  << "\n";
    cout << "  output.bin        = " << output_bin << "\n";
    cout << "  K                 = " << K << "\n";
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
    unordered_map<uint64_t, node_id> idmap;
    idmap.reserve(1 << 20);

    node_id next_id = 0;
    auto get_new_id = [&](uint64_t orig) -> node_id {
        auto it = idmap.find(orig);
        if (it != idmap.end()) return it->second;
        node_id nid = next_id++;
        idmap.emplace(orig, nid);
        return nid;
    };

    // Dedup THO (keep-first): key=(u,v) -> weights[K]
    unordered_map<uint64_t, std::vector<double>> keep_first;
    keep_first.reserve(1 << 20);

    size_t edge_line_no = 0;
    size_t num_lines = 0;
    size_t num_kept_directed = 0;
    size_t num_skipped_directed = 0;

    auto try_insert_keep_first = [&](node_id a, node_id b, const std::vector<double>& ww) {
        const uint64_t key = pack_uv(static_cast<uint32_t>(a), static_cast<uint32_t>(b));
        auto it = keep_first.find(key);
        if (it != keep_first.end()) {
            ++num_skipped_directed;
            return;
        }
        keep_first.emplace(key, ww);
        ++num_kept_directed;
    };

    string line;
    while (std::getline(fin, line)) {
        ++edge_line_no;
        if (line.empty() || line[0] == '#') continue;

        ++num_lines;

        stringstream ss(line);
        uint64_t u0, v0;
        if (!(ss >> u0 >> v0)) continue;

        node_id u = get_new_id(u0);
        node_id v = get_new_id(v0);

        // parse weights
        vector<double> w(K, 1.0);   // default 1.0 if no weights provided
        vector<double> tmp;
        tmp.reserve(K);

        double x;
        while (ss >> x) tmp.push_back(x);

        if (!tmp.empty()) {
            if (tmp.size() == 1) {
                double val = clamp_nonneg(tmp[0]); // negative -> 0
                std::fill(w.begin(), w.end(), val);
            } else if (tmp.size() >= K) {
                for (size_t t = 0; t < K; ++t) w[t] = clamp_nonneg(tmp[t]);
            } else {
                // thiếu cột -> pad 0
                for (size_t t = 0; t < tmp.size(); ++t) w[t] = clamp_nonneg(tmp[t]);
                for (size_t t = tmp.size(); t < K; ++t) w[t] = 0.0;
            }
        } else {
            // keep default 1.0
            for (size_t t = 0; t < K; ++t) w[t] = 1.0;
        }

        // Output ALWAYS directed, dedup THO
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

    cout << "  Renumbered nodes n = " << n << "\n";
    cout << "  lines_read         = " << num_lines << "\n";
    cout << "  kept_directed      = " << num_kept_directed << "\n";
    cout << "  skipped_duplicates = " << num_skipped_directed << "\n";

    // --------- 2) Read parts.txt (FORMAT A ONLY) ----------
    ifstream fp(parts_txt);
    if (!fp) {
        cerr << "Error: cannot open parts file: " << parts_txt << "\n";
        return 1;
    }

    unordered_map<uint64_t, uint32_t> part_map_orig; // orig -> raw part_id
    part_map_orig.reserve(n * 2);

    vector<uint32_t> raw_part_ids;
    raw_part_ids.reserve(n);

    size_t parts_line_no = 0;
    while (std::getline(fp, line)) {
        ++parts_line_no;
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        unsigned long long a, b;
        if (!(ss >> a >> b)) {
            cerr << "Error: parts.txt invalid format at line " << parts_line_no
                 << " (expected: orig_node_id part_id)\n";
            return 1;
        }
        std::string extra;
        if (ss >> extra) {
            cerr << "Error: parts.txt has extra token at line " << parts_line_no
                 << " (expected exactly 2 columns)\n";
            return 1;
        }

        const uint64_t orig = static_cast<uint64_t>(a);
        const uint32_t pid_raw = static_cast<uint32_t>(b);

        if (part_map_orig.find(orig) != part_map_orig.end()) {
            cerr << "Error: parts.txt duplicate orig_node_id = " << orig
                 << " at line " << parts_line_no << "\n";
            return 1;
        }

        part_map_orig.emplace(orig, pid_raw);
        raw_part_ids.push_back(pid_raw);
    }
    fp.close();

    if (part_map_orig.empty()) {
        cerr << "Error: parts.txt has no valid data lines.\n";
        return 1;
    }

    // --------- 3) Strict consistency checks ----------
    // (i) Every orig id in parts must appear in edges
    for (const auto& kv : part_map_orig) {
        const uint64_t orig = kv.first;
        if (idmap.find(orig) == idmap.end()) {
            cerr << "Error: parts.txt contains orig_node_id = " << orig
                 << " but this node never appears in edges.txt\n";
            return 1;
        }
    }

    // (ii) Every node appearing in edges must have a part_id
    if (part_map_orig.size() != idmap.size()) {
        for (const auto& kv : idmap) {
            const uint64_t orig = kv.first;
            if (part_map_orig.find(orig) == part_map_orig.end()) {
                cerr << "Error: missing part_id for orig_node_id = " << orig
                     << " (node appears in edges.txt but not in parts.txt)\n";
                return 1;
            }
        }
        cerr << "Error: mismatch between edges nodes and parts nodes.\n";
        return 1;
    }

    cout << "  parts assigned     = " << part_map_orig.size() << " / " << n << "\n";

    // --------- 4) Renumber part_id to 0..p-1 (compress) ----------
    std::sort(raw_part_ids.begin(), raw_part_ids.end());
    raw_part_ids.erase(std::unique(raw_part_ids.begin(), raw_part_ids.end()), raw_part_ids.end());

    const size_t p = raw_part_ids.size();
    unordered_map<uint32_t, uint32_t> part_renum; // raw -> new
    part_renum.reserve(p * 2);

    for (size_t i = 0; i < p; ++i) {
        part_renum.emplace(raw_part_ids[i], static_cast<uint32_t>(i));
    }

    cout << "  unique parts (raw) = " << p << "\n";
    cout << "  part_id renumber   : raw_min=" << raw_part_ids.front()
         << ", raw_max=" << raw_part_ids.back()
         << "  -> new range [0.." << (p ? (p - 1) : 0) << "]\n";

    // Apply renumber in part_map_orig (orig -> new pid)
    for (auto &kv : part_map_orig) {
        uint32_t raw = kv.second;
        auto it = part_renum.find(raw);
        if (it == part_renum.end()) {
            cerr << "Internal error: missing renumber mapping for raw part_id=" << raw << "\n";
            return 1;
        }
        kv.second = it->second;
    }

    // --------- 5) Build directed tinyGraph (NEW format) ----------
    tinyGraph g;
    g.init(n, K, /*undirected_flag=*/false); // ALWAYS directed

    // part_id (already renumbered)
    g.part_id.assign(n, 0);
    for (const auto& kv : part_map_orig) {
        const uint64_t orig = kv.first;
        const uint32_t pid  = kv.second; // renumbered
        const node_id u = idmap[orig];
        g.part_id[static_cast<size_t>(u)] = pid;
    }

    // node weights/alpha
    if (randomize_node) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(1e-6, 1.0);
        for (size_t i = 0; i < g.n; ++i) {
            for (size_t t = 0; t < K; ++t) {
                g.nodes[i].weights[t] = dist(gen);
                g.nodes[i].alpha[t]   = dist(gen);
            }
        }
    } // else already initialized by g.init()

    // materialize edges in deterministic order (sort by (u,v))
    std::vector<std::uint64_t> keys;
    keys.reserve(keep_first.size());
    for (const auto &kv : keep_first) keys.push_back(kv.first);

    std::sort(keys.begin(), keys.end(), [](std::uint64_t a, std::uint64_t b) {
        const std::uint32_t au = unpack_u(a), av = unpack_v(a);
        const std::uint32_t bu = unpack_u(b), bv = unpack_v(b);
        if (au != bu) return au < bu;
        return av < bv;
    });

    g.edges.clear();
    g.edges.reserve(keys.size());
    for (std::uint64_t key : keys) {
        const uint32_t u = unpack_u(key);
        const uint32_t v = unpack_v(key);

        Edge E;
        E.u = static_cast<node_id>(u);
        E.v = static_cast<node_id>(v);

        auto it = keep_first.find(key);
        if (it != keep_first.end()) {
            E.weights = it->second;
        } else {
            E.weights.assign(K, 0.0);
        }
        if (E.weights.size() != K) E.weights.resize(K, 0.0);

        g.edges.push_back(std::move(E));
    }
    g.m = g.edges.size();

    // --------- 6) Normalize incoming per topic ----------
    // sum_in[v,t] = sum_{(u->v)} w_t(u,v) ; if sum_in>0 then divide to make sum=1
    std::vector<double> sum_in(g.n * g.K, 0.0);

    for (edge_id eid = 0; eid < g.m; ++eid) {
        const auto &E = g.edges[eid];
        const node_id v = E.v;
        if (v >= g.n) continue;
        for (size_t t = 0; t < g.K; ++t) {
            double wt = (t < E.weights.size() ? E.weights[t] : 0.0);
            if (wt > 0.0) sum_in[static_cast<size_t>(v) * g.K + t] += wt;
        }
    }

    for (edge_id eid = 0; eid < g.m; ++eid) {
        auto &E = g.edges[eid];
        const node_id v = E.v;
        if (v >= g.n) continue;
        for (size_t t = 0; t < g.K; ++t) {
            const double s = sum_in[static_cast<size_t>(v) * g.K + t];
            if (s > 0.0) {
                E.weights[t] = E.weights[t] / s;
                if (E.weights[t] < 0.0) E.weights[t] = 0.0; // safety
            }
        }
    }

    // --------- 7) Write binary ----------
    if (!g.write_binary(output_bin)) {
        cerr << "Error: write_binary failed: " << output_bin << "\n";
        return 1;
    }

    cout << "Done. Saved binary graph to: " << output_bin << "\n";
    cout << "Summary:\n";
    cout << "  n = " << g.n << ", m = " << g.m << ", K = " << g.K
         << ", undirected(stored) = " << (g.undirected ? "true" : "false") << "\n";
    cout << "  p = " << p << " (number of parts after renumber)\n";
    cout << "  dedup = keep-first, dropped " << num_skipped_directed << " duplicates\n";
    return 0;
}

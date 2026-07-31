// preprocess.cpp  (Partition-aware, FORMAT A ONLY: orig_node_id part_id)
// + Renumber (compress) part_id to 0,1,2,... by ascending original part_id
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

static inline double sanitize_w(double w) {
    return (w > 0.0) ? w : 1e-6;
}

static void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog
        << " edges.txt parts.txt output.bin K [undirected(0/1)=1] [randomize_node(0/1)=0] [seed=42]\n\n"
        << "edges.txt format:\n"
        << "  - Moi dong: u v [w]\n"
        << "    + u, v: so nguyen khong am (ID goc), co the khong lien tiep\n"
        << "    + w   : tuy chon; neu thieu -> 1.0; neu w<=0 -> sanitize 1e-6\n"
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
        << "  - parts.txt co node khong xuat hien trong edges.txt -> loi.\n";
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

    bool undirected = true;
    bool randomize_node = false;
    unsigned seed = 42;

    if (argc >= 6) undirected = (std::atoi(argv[5]) != 0);
    if (argc >= 7) randomize_node = (std::atoi(argv[6]) != 0);
    if (argc >= 8) seed = static_cast<unsigned>(std::strtoul(argv[7], nullptr, 10));

    cout << "Preprocess (FORMAT A ONLY, part_id will be renumbered):\n";
    cout << "  edges.txt      = " << edges_txt << "\n";
    cout << "  parts.txt      = " << parts_txt << "\n";
    cout << "  output.bin     = " << output_bin << "\n";
    cout << "  K              = " << K << "\n";
    cout << "  undirected     = " << (undirected ? "true" : "false") << "\n";
    cout << "  randomize_node = " << (randomize_node ? "true" : "false") << "\n";
    cout << "  seed           = " << seed << "\n";

    // --------- 1) Read edges and build renumber map ----------
    ifstream fin(edges_txt);
    if (!fin) {
        cerr << "Error: cannot open edges file: " << edges_txt << "\n";
        return 1;
    }

    struct RawEdge { uint64_t u0, v0; double w; };
    vector<RawEdge> rawEdges;
    rawEdges.reserve(1 << 20);

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

    string line;
    size_t edge_line_no = 0;
    while (std::getline(fin, line)) {
        ++edge_line_no;
        if (line.empty() || line[0] == '#') continue;

        uint64_t u0, v0;
        double w = 1.0;
        std::stringstream ss(line);
        ss >> u0 >> v0;
        if (!ss) continue;

        if (ss >> w) w = sanitize_w(w);

        rawEdges.push_back({u0, v0, w});
        (void)get_new_id(u0);
        (void)get_new_id(v0);
    }
    fin.close();

    const size_t n = static_cast<size_t>(next_id);
    if (n == 0 || K == 0) {
        cerr << "Error: empty graph or K=0.\n";
        return 1;
    }

    cout << "  Renumbered nodes n = " << n << "\n";
    cout << "  Raw edges read     = " << rawEdges.size() << "\n";

    // --------- 2) Read parts.txt (FORMAT A ONLY) ----------
    ifstream fp(parts_txt);
    if (!fp) {
        cerr << "Error: cannot open parts file: " << parts_txt << "\n";
        return 1;
    }

    unordered_map<uint64_t, uint32_t> part_map_orig; // orig -> part (raw)
    part_map_orig.reserve(n * 2);

    vector<uint32_t> raw_part_ids; // collect for renumber
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
    // Deterministic: sort unique raw part ids ascending
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
            // should never happen
            cerr << "Internal error: missing renumber mapping for raw part_id=" << raw << "\n";
            return 1;
        }
        kv.second = it->second;
    }

    // --------- 5) Build tinyGraph (NEW) ----------
    tinyGraph g;
    g.init(n, K, undirected);

    // part_id (now already renumbered)
    g.part_id.assign(n, 0);
    for (const auto& kv : part_map_orig) {
        const uint64_t orig = kv.first;
        const uint32_t pid  = kv.second; // renumbered
        const node_id u = idmap[orig];
        g.part_id[static_cast<size_t>(u)] = pid;
    }

    // nodes weights/alpha
    if (randomize_node) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(1e-6, 1.0);
        for (size_t i = 0; i < n; ++i) {
            for (size_t t = 0; t < K; ++t) {
                g.nodes[i].weights[t] = dist(gen);
                g.nodes[i].alpha[t]   = dist(gen);
            }
        }
    } // else already 1.0 by init()

    // edges
    g.edges.clear();
    g.edges.reserve(rawEdges.size());
    for (const auto& re : rawEdges) {
        const node_id u = idmap[re.u0];
        const node_id v = idmap[re.v0];
        Edge E;
        E.u = u;
        E.v = v;
        E.weights.assign(K, re.w);
        g.edges.push_back(std::move(E));
    }
    g.m = g.edges.size();

    // write .bin (NEW format)
    if (!g.write_binary(output_bin)) {
        cerr << "Error: write_binary failed: " << output_bin << "\n";
        return 1;
    }

    cout << "Done. Saved binary graph to: " << output_bin << "\n";
    cout << "Summary:\n";
    cout << "  n = " << g.n << ", m = " << g.m << ", K = " << g.K
         << ", undirected = " << (g.undirected ? "true" : "false") << "\n";
    cout << "  p = " << p << " (number of parts after renumber)\n";

    return 0;
}

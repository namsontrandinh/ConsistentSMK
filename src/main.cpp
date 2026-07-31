// src/main.cpp
// EDL-only executable for scalar submodular maximization
// under a knapsack constraint.
//
// Example:
// ./ic --graph fb.bin --B_factor 0.01 --alg edl --eps 0.1 --consistency 1 --seed 42 --mc 100 --csv fb_edl.csv

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <sys/resource.h>
#include <sys/stat.h>

#include "mygraph.h"
#include "algs/edl.h"

static long getPeakRSS_KB() {
    struct rusage r;

    if (getrusage(
            RUSAGE_SELF,
            &r) == 0)
    {
        return r.ru_maxrss;
    }

    return 0;
}

static bool file_is_empty_or_missing(
    const std::string& path)
{
    struct stat st;

    if (stat(
            path.c_str(),
            &st) != 0)
    {
        return true;
    }

    return st.st_size == 0;
}

static std::string csv_escape(
    const std::string& s)
{
    bool need_quote = false;

    for (char c : s) {
        if (c == ',' ||
            c == '"' ||
            c == '\n' ||
            c == '\r')
        {
            need_quote = true;
            break;
        }
    }

    if (!need_quote) {
        return s;
    }

    std::string out;
    out.reserve(
        s.size() + 2);

    out.push_back('"');

    for (char c : s) {
        if (c == '"') {
            out.push_back('"');
        }

        out.push_back(c);
    }

    out.push_back('"');

    return out;
}

static void append_csv_row(
    const std::string& csv_path,
    const std::string& header,
    const std::string& row)
{
    if (csv_path.empty()) {
        return;
    }

    const bool write_header =
        file_is_empty_or_missing(
            csv_path);

    std::ofstream fout(
        csv_path,
        std::ios::out |
        std::ios::app);

    if (!fout.good()) {
        std::cerr
            << "[WARN] Cannot open CSV: "
            << csv_path
            << "\n";
        return;
    }

    if (write_header) {
        fout
            << header
            << "\n";
    }

    fout
        << row
        << "\n";
}

static void print_usage(
    const char* prog)
{
    std::cerr
        << "Usage:\n"
        << "  "
        << prog
        << " --graph <graph.bin>"
        << " --B_factor <double>"
        << " --alg edl"
        << " [--eps <double>]"
        << " [--consistency <0|1>]"
        << " [--seed <uint64>]"
        << " [--mc <size_t>]"
        << " [--lambda <double>]"
        << " [--csv <output.csv>]\n";
}

int main(
    int argc,
    char** argv)
{
    using namespace std;
    using namespace mygraph;

    string graph_file;
    string algo = "edl";
    string csv_path;

    double B_factor = -1.0;
    double eps = 0.1;

    bool compute_consistency = true;

    string ic_seed = "42";
    string ic_mc = "100";
    string ic_lambda = "1.0";

    for (int i = 1;
         i < argc;
         ++i)
    {
        const string tok =
            argv[i];

        auto need_next =
            [&](const char* opt)
                -> const char*
        {
            if (i + 1 >= argc) {
                cerr
                    << "Error: missing value for "
                    << opt
                    << "\n";

                print_usage(
                    argv[0]);

                std::exit(1);
            }

            return argv[++i];
        };

        if (tok == "--graph") {
            graph_file =
                need_next(
                    "--graph");
        }
        else if (tok == "--B_factor") {
            B_factor =
                std::strtod(
                    need_next(
                        "--B_factor"),
                    nullptr);
        }
        else if (tok == "--alg") {
            algo =
                need_next(
                    "--alg");
        }
        else if (tok == "--eps") {
            eps =
                std::strtod(
                    need_next(
                        "--eps"),
                    nullptr);
        }
        else if (tok == "--consistency") {
            compute_consistency =
                std::atoi(
                    need_next(
                        "--consistency")) != 0;
        }
        else if (tok == "--seed") {
            ic_seed =
                need_next(
                    "--seed");
        }
        else if (tok == "--mc") {
            ic_mc =
                need_next(
                    "--mc");
        }
        else if (tok == "--lambda") {
            ic_lambda =
                need_next(
                    "--lambda");
        }
        else if (tok == "--csv") {
            csv_path =
                need_next(
                    "--csv");
        }
        else if (tok == "--help" ||
                 tok == "-h")
        {
            print_usage(
                argv[0]);

            return 0;
        }
        else {
            cerr
                << "Error: unknown argument: "
                << tok
                << "\n";

            print_usage(
                argv[0]);

            return 1;
        }
    }

    if (graph_file.empty()) {
        cerr
            << "Error: --graph is required.\n";
        return 1;
    }

    if (!(B_factor > 0.0)) {
        cerr
            << "Error: --B_factor must be > 0.\n";
        return 1;
    }

    if (!(eps > 0.0 &&
          eps < 1.0))
    {
        cerr
            << "Error: --eps must be in (0,1).\n";
        return 1;
    }

    if (algo != "edl" &&
        algo != "EDL")
    {
        cerr
            << "Error: this executable only supports EDL.\n";
        return 1;
    }

    // Fix the IC Monte-Carlo environment for reproducibility.
    if (::setenv(
            "KIC_SEED",
            ic_seed.c_str(),
            1) != 0)
    {
        cerr
            << "Error: cannot set KIC_SEED.\n";
        return 1;
    }

    if (::setenv(
            "KIC_MC",
            ic_mc.c_str(),
            1) != 0)
    {
        cerr
            << "Error: cannot set KIC_MC.\n";
        return 1;
    }

    if (::setenv(
            "KIC_LAMBDA",
            ic_lambda.c_str(),
            1) != 0)
    {
        cerr
            << "Error: cannot set KIC_LAMBDA.\n";
        return 1;
    }

    tinyGraph g;

    if (!g.read_binary(
            graph_file))
    {
        cerr
            << "Error: cannot read graph: "
            << graph_file
            << "\n";

        return 1;
    }

    if (g.K == 0) {
        cerr
            << "Error: graph has K=0.\n";
        return 1;
    }

    double total_cost = 0.0;

    for (std::size_t u = 0;
         u < g.n;
         ++u)
    {
        const double c =
            subm_obj_common::node_cost(
                g,
                static_cast<node_id>(u));

        if (c > 0.0) {
            total_cost += c;
        }
    }

    if (!(total_cost > 0.0)) {
        cerr
            << "Error: total positive node cost is zero.\n";
        return 1;
    }

    const double B =
        B_factor *
        total_cost;

    cout
        << fixed
        << setprecision(6);

    cout
        << "========================================\n"
        << "EDL experiment\n"
        << "========================================\n"
        << "graph                  = "
        << graph_file
        << "\n"
        << "n                      = "
        << g.n
        << "\n"
        << "m                      = "
        << g.m
        << "\n"
        << "K                      = "
        << g.K
        << "\n"
        << "B_factor               = "
        << B_factor
        << "\n"
        << "total node cost        = "
        << total_cost
        << "\n"
        << "B                      = "
        << B
        << "\n"
        << "eps                    = "
        << eps
        << "\n"
        << "KIC_SEED               = "
        << ic_seed
        << "\n"
        << "KIC_MC                 = "
        << ic_mc
        << "\n"
        << "KIC_LAMBDA             = "
        << ic_lambda
        << "\n"
        << "compute consistency    = "
        << (compute_consistency
                ? "true"
                : "false")
        << "\n";

    const long peak_before_kb =
        getPeakRSS_KB();

    algs::Result res;

    std::uint64_t C = 0;
    std::uint64_t cumulative_consistency = 0;
    std::uint64_t initial_changes = 0;
    std::uint64_t consistency_queries = 0;

    double consistency_time_sec = 0.0;

    if (compute_consistency) {
        const algs::EDLConsistencyStats stats =
            algs::run_EDL_consistency(
                g,
                B,
                eps);

        res =
            stats.final_result;

        C =
            stats.C;

        cumulative_consistency =
            stats.cumulative_consistency;

        initial_changes =
            stats.initial_changes;

        consistency_queries =
            stats.total_queries;

        consistency_time_sec =
            stats.total_time_sec;
    }
    else {
        res =
            algs::run_EDL(
                g,
                B,
                eps);
    }

    const long peak_after_kb =
        getPeakRSS_KB();

    if (peak_after_kb >
        peak_before_kb)
    {
        res.mem_mb =
            static_cast<double>(
                peak_after_kb -
                peak_before_kb) /
            1024.0;
    }
    else {
        res.mem_mb = 0.0;
    }

    if (res.inS.size() != g.n) {
        cerr
            << "Error: result solution size "
            << res.inS.size()
            << " differs from graph size "
            << g.n
            << ".\n";

        return 1;
    }

    double solution_cost = 0.0;
    std::size_t solution_size = 0;

    for (std::size_t u = 0;
         u < g.n;
         ++u)
    {
        if (res.inS[u]) {
            solution_cost +=
                subm_obj_common::node_cost(
                    g,
                    static_cast<node_id>(u));

            ++solution_size;
        }
    }

    cout
        << "----------------------------------------\n"
        << res.algo
        << " finished.\n"
        << "constraint              = "
        << res.constraint
        << "\n"
        << "f(S)                    = "
        << res.f_value
        << "\n"
        << "c(S)                    = "
        << solution_cost
        << "\n"
        << "|S|                     = "
        << solution_size
        << "\n"
        << "#calls final prefix     = "
        << res.queries
        << "\n"
        << "final-prefix time (s)   = "
        << res.time_sec
        << "\n"
        << "memory used (MB)        = "
        << res.mem_mb
        << "\n";

    if (compute_consistency) {
        cout
            << "C-consistency           = "
            << C
            << "\n"
            << "cumulative consistency  = "
            << cumulative_consistency
            << "\n"
            << "initial changes         = "
            << initial_changes
            << "\n"
            << "all-prefix queries      = "
            << consistency_queries
            << "\n"
            << "all-prefix time (s)     = "
            << consistency_time_sec
            << "\n";
    }

    if (!csv_path.empty()) {
        const string header =
            "algo,constraint,graph,n,m,K,"
            "B_factor,total_cost,B,eps,"
            "ic_seed,ic_mc,ic_lambda,"
            "f_value,solution_cost,solution_size,"
            "queries,time_sec,mem_mb,"
            "C,cumulative_consistency,initial_changes,"
            "consistency_queries,consistency_time_sec";

        ostringstream row;

        row
            << csv_escape(
                   res.algo)
            << ","
            << csv_escape(
                   res.constraint)
            << ","
            << csv_escape(
                   graph_file)
            << ","
            << g.n
            << ","
            << g.m
            << ","
            << g.K
            << ","
            << setprecision(17)
            << B_factor
            << ","
            << total_cost
            << ","
            << B
            << ","
            << eps
            << ","
            << csv_escape(
                   ic_seed)
            << ","
            << csv_escape(
                   ic_mc)
            << ","
            << csv_escape(
                   ic_lambda)
            << ","
            << res.f_value
            << ","
            << solution_cost
            << ","
            << solution_size
            << ","
            << res.queries
            << ","
            << setprecision(10)
            << res.time_sec
            << ","
            << setprecision(6)
            << res.mem_mb
            << ","
            << C
            << ","
            << cumulative_consistency
            << ","
            << initial_changes
            << ","
            << consistency_queries
            << ","
            << setprecision(10)
            << consistency_time_sec;

        append_csv_row(
            csv_path,
            header,
            row.str());

        cout
            << "CSV appended to         = "
            << csv_path
            << "\n";
    }

    return 0;
}

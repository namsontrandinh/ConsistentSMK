// src/main.cpp
// Executable for scalar submodular maximization under a knapsack
// constraint. Supports --alg edl and --alg atkc.
//
// Examples:
// ./ic --graph fb.bin --B_factor 0.01 --alg edl --eps 0.1 --consistency 1 --seed 42 --mc 100 --csv fb_edl.csv
// ./ic --graph fb.bin --B_factor 0.01 --alg atkc --delta 0.1 --atkc_eps 0.5 --consistency 1 --seed 42 --mc 100 --csv fb_atkc.csv
// ./ic --graph fb.bin --B_abs 80 --alg rkcl --delta 0.1 --consistency 1 --seed 42 --mc 100 --csv fb_rkcl.csv

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>

#include <sys/resource.h>
#include <sys/stat.h>

#include "mygraph.h"
#include "algs/edl.h"
#include "algs/atkc.h"
#include "algs/rkcl_im.h"

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
        << " (--B_factor <double> | --B_abs <double>)"
        << " --alg <edl|atkc|rkcl>"
        << " [--eps <double>]           (EDL only)\n"
        << "  [--delta <double>]         (ATKC/RKCL, cost-granularity; ATKC default 0.1, RKCL default 0 = auto-detect)\n"
        << "  [--atkc_eps <double>]      (ATKC only, precision, default 0.5)\n"
        << "  [--chase_cap <uint64>]     (ATKC/RKCL, 0 = theoretical N)\n"
        << "  [--max_scan <uint64>]      (ATKC only, 0 = unlimited)\n"
        << "  [--active_n <size_t>]      (ATKC/RKCL, 0 = full graph; truncate stream for quick testing)\n"
        << "  [--consistency <0|1>]"
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
    double eps = 0.1;          // used by EDL

    double atkc_delta = 0.1;   // used by ATKC (Algorithm 4, Line 1)
    double atkc_eps   = 0.5;   // used by ATKC (Algorithm 4, Line 1)
    std::uint64_t atkc_chase_cap     = 0;  // 0 = theoretical N
    std::uint64_t atkc_max_scan_iter = 0;  // 0 = unlimited candidate scan

    double rkcl_delta = 0.0;   // used by RKCL (Algorithm 2, Line 1). 0 = auto-detect.
    std::uint64_t rkcl_chase_cap = 0;      // 0 = theoretical N
    bool rkcl_strict_granularity = true;

    double B_abs = -1.0;       // absolute budget; overrides B_factor*total_cost when > 0.
    std::size_t active_n_override = 0; // 0 = use full graph (g.n)

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
        else if (tok == "--delta") {
            atkc_delta =
                std::strtod(
                    need_next(
                        "--delta"),
                    nullptr);
            rkcl_delta = atkc_delta;
        }
        else if (tok == "--atkc_eps") {
            atkc_eps =
                std::strtod(
                    need_next(
                        "--atkc_eps"),
                    nullptr);
        }
        else if (tok == "--chase_cap") {
            atkc_chase_cap =
                std::strtoull(
                    need_next(
                        "--chase_cap"),
                    nullptr,
                    10);
            rkcl_chase_cap = atkc_chase_cap;
        }
        else if (tok == "--B_abs") {
            B_abs =
                std::strtod(
                    need_next(
                        "--B_abs"),
                    nullptr);
        }
        else if (tok == "--active_n") {
            active_n_override =
                static_cast<std::size_t>(
                    std::strtoull(
                        need_next(
                            "--active_n"),
                        nullptr,
                        10));
        }
        else if (tok == "--max_scan") {
            atkc_max_scan_iter =
                std::strtoull(
                    need_next(
                        "--max_scan"),
                    nullptr,
                    10);
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

    if (!(B_factor > 0.0) && !(B_abs > 0.0)) {
        cerr
            << "Error: either --B_factor or --B_abs must be > 0.\n";
        return 1;
    }

    if (!(eps > 0.0 &&
          eps < 1.0))
    {
        cerr
            << "Error: --eps must be in (0,1).\n";
        return 1;
    }

    // Normalize algo name to lowercase for comparison.
    string algo_lc = algo;
    std::transform(
        algo_lc.begin(),
        algo_lc.end(),
        algo_lc.begin(),
        [](unsigned char c) { return std::tolower(c); });

    const bool is_edl  = (algo_lc == "edl");
    const bool is_atkc = (algo_lc == "atkc");
    const bool is_rkcl = (algo_lc == "rkcl");

    if (!is_edl && !is_atkc && !is_rkcl) {
        cerr
            << "Error: unsupported --alg '"
            << algo
            << "'. Supported: edl, atkc, rkcl.\n";
        return 1;
    }

    if (is_atkc) {
        if (!(atkc_delta > 0.0 && atkc_delta <= 1.0)) {
            cerr << "Error: --delta must be in (0,1].\n";
            return 1;
        }
        if (!(atkc_eps > 0.0 && atkc_eps < 1.0)) {
            cerr << "Error: --atkc_eps must be in (0,1).\n";
            return 1;
        }
    }

    if (is_rkcl) {
        if (rkcl_delta != 0.0 &&
            !(rkcl_delta > 0.0 && rkcl_delta <= 1.0))
        {
            cerr << "Error: --delta must be in (0,1] (or 0 for auto-detect).\n";
            return 1;
        }
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
        (B_abs > 0.0)
            ? B_abs
            : (B_factor * total_cost);

    cout
        << fixed
        << setprecision(6);

    cout
        << "========================================\n"
        << (is_atkc ? "ATKC experiment\n" : (is_rkcl ? "RKCL experiment\n" : "EDL experiment\n"))
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
        << "\n";

    if (B_abs > 0.0) {
        cout
            << "B_abs                  = "
            << B_abs
            << "\n";
    } else {
        cout
            << "B_factor               = "
            << B_factor
            << "\n";
    }

    cout
        << "total node cost        = "
        << total_cost
        << "\n"
        << "B                      = "
        << B
        << "\n";

    if (is_edl) {
        cout
            << "eps                    = "
            << eps
            << "\n";
    }
    if (is_atkc) {
        cout
            << "delta                  = "
            << atkc_delta
            << "\n"
            << "atkc_eps               = "
            << atkc_eps
            << "\n"
            << "chase_cap              = "
            << atkc_chase_cap
            << "\n"
            << "max_scan               = "
            << atkc_max_scan_iter
            << "\n";
    }
    if (is_rkcl) {
        cout
            << "delta (0=auto)         = "
            << rkcl_delta
            << "\n"
            << "chase_cap              = "
            << rkcl_chase_cap
            << "\n";
    }

    cout
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

    // ATKC-specific diagnostics (left at 0 for other branches so the
    // shared CSV schema stays numeric/parseable regardless of --alg).
    std::uint64_t atkc_q_param = 0;
    std::uint64_t atkc_J_param = 0;
    std::uint64_t atkc_N_param = 0;
    std::uint64_t atkc_total_lane_restarts = 0;
    std::uint64_t atkc_total_successful_exchanges = 0;

    // RKCL-specific diagnostics.
    double rkcl_delta_used = 0.0;
    double rkcl_mu_used = 0.0;
    double rkcl_lambda_used = 0.0;
    double rkcl_gamma_used = 0.0;
    std::uint64_t rkcl_q_param = 0;
    std::uint64_t rkcl_N_param = 0;
    std::uint64_t rkcl_total_restarts = 0;
    std::uint64_t rkcl_total_successful_exchanges = 0;

    const std::size_t active_n =
        (active_n_override > 0)
            ? active_n_override
            : g.n;

    if (is_edl) {
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
    }
    else if (is_atkc) {
        algs::ATKCParams params;
        params.delta = atkc_delta;
        params.eps   = atkc_eps;
        params.chase_cap = atkc_chase_cap;
        params.max_x_scan_per_iter = atkc_max_scan_iter;

        if (compute_consistency) {
            const algs::ATKCConsistencyStats stats =
                algs::run_ATKC_consistency_prefix(
                    g,
                    B,
                    params,
                    active_n);

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

            atkc_q_param = stats.params.q;
            atkc_J_param = stats.params.J;
            atkc_N_param = stats.params.used_N;
            atkc_total_lane_restarts = stats.total_lane_restarts;
            atkc_total_successful_exchanges = stats.total_successful_exchanges;
        }
        else {
            res =
                algs::run_ATKC_prefix(
                    g,
                    B,
                    params,
                    active_n);
        }
    }
    else if (is_rkcl) {
        algs::RKCLParams params;
        params.delta = rkcl_delta;
        params.chase_cap = rkcl_chase_cap;
        params.strict_granularity = rkcl_strict_granularity;

        if (compute_consistency) {
            const algs::RKCLConsistencyStats stats =
                algs::run_RKCL_consistency_prefix(
                    g,
                    B,
                    params,
                    active_n);

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

            rkcl_delta_used = stats.params.delta;
            rkcl_mu_used = stats.params.mu;
            rkcl_lambda_used = stats.params.lambda;
            rkcl_gamma_used = stats.params.gamma;
            rkcl_q_param = stats.params.q;
            rkcl_N_param = stats.params.used_N;
            rkcl_total_restarts = stats.total_restarts;
            rkcl_total_successful_exchanges = stats.total_successful_exchanges;
        }
        else {
            res =
                algs::run_RKCL_prefix(
                    g,
                    B,
                    params,
                    active_n);
        }
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

        if (is_atkc) {
            cout
                << "q (ATKC)                = "
                << atkc_q_param
                << "\n"
                << "J (ATKC, lanes=J+1)     = "
                << atkc_J_param
                << "\n"
                << "N used (ATKC)           = "
                << atkc_N_param
                << "\n"
                << "lane restarts           = "
                << atkc_total_lane_restarts
                << "\n"
                << "successful exchanges    = "
                << atkc_total_successful_exchanges
                << "\n";
        }
        if (is_rkcl) {
            cout
                << "delta used (RKCL)       = "
                << rkcl_delta_used
                << "\n"
                << "mu                      = "
                << rkcl_mu_used
                << "\n"
                << "lambda                  = "
                << rkcl_lambda_used
                << "\n"
                << "gamma (1+delta)         = "
                << rkcl_gamma_used
                << "\n"
                << "q (RKCL)                = "
                << rkcl_q_param
                << "\n"
                << "N used (RKCL)           = "
                << rkcl_N_param
                << "\n"
                << "restarts                = "
                << rkcl_total_restarts
                << "\n"
                << "successful exchanges    = "
                << rkcl_total_successful_exchanges
                << "\n";
        }
    }

    if (!csv_path.empty()) {
        const string header =
            "algo,constraint,graph,n,m,K,"
            "B_factor,B_abs,total_cost,B,eps,"
            "ic_seed,ic_mc,ic_lambda,"
            "f_value,solution_cost,solution_size,"
            "queries,time_sec,mem_mb,"
            "C,cumulative_consistency,initial_changes,"
            "consistency_queries,consistency_time_sec,"
            "atkc_delta,atkc_eps,atkc_chase_cap,atkc_max_scan,"
            "atkc_q,atkc_J,atkc_N,"
            "atkc_lane_restarts,atkc_successful_exchanges,"
            "rkcl_delta_used,rkcl_mu,rkcl_lambda,rkcl_gamma,"
            "rkcl_q,rkcl_N,rkcl_chase_cap,"
            "rkcl_restarts,rkcl_successful_exchanges";

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
            << B_abs
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
            << consistency_time_sec
            << ","
            << setprecision(17)
            << atkc_delta
            << ","
            << atkc_eps
            << ","
            << atkc_chase_cap
            << ","
            << atkc_max_scan_iter
            << ","
            << atkc_q_param
            << ","
            << atkc_J_param
            << ","
            << atkc_N_param
            << ","
            << atkc_total_lane_restarts
            << ","
            << atkc_total_successful_exchanges
            << ","
            << setprecision(17)
            << rkcl_delta_used
            << ","
            << rkcl_mu_used
            << ","
            << rkcl_lambda_used
            << ","
            << rkcl_gamma_used
            << ","
            << rkcl_q_param
            << ","
            << rkcl_N_param
            << ","
            << rkcl_chase_cap
            << ","
            << rkcl_total_restarts
            << ","
            << rkcl_total_successful_exchanges;

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
#include <CLI/CLI.hpp>
#include <fstream>
#include <iomanip>
#include "cli.hpp"

#if defined(_WIN32) // Windows API
#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <Windows.h>
#define SYMBA_OSAPI_WINDOWS
#elif defined(__unix__) or defined(__APPLE__) // POSIX
#include <sys/ioctl.h>
#define SYMBA_OSAPI_POSIX
#else // Something else?
#error CLI module could not determine OS interface.
#endif

namespace cli {

void init(int argc_, char** argv_) {
    argc = argc_;
    argv = argv_;
    
    using namespace CLI;
    App app("A market simulator built in C++.", "symba");

    using namespace args;
    app.add_option("-o,--output-dir", output_dir, "Output directory for simulation results")
        ->required()
        ->check(ExistingDirectory);
    app.add_option("-I,--n-agents", n_agents, "Number of agents in the simulation")
        ->default_val(500)
        ->check(Range(10, 10000));
    app.add_option("-J,--n-stocks", n_stocks, "Number of stocks in the simulation")
        ->default_val(1)
        ->check(Range(1, 100));
    app.add_option("-T,--n-steps", n_steps, "Number of simulation time steps")
        ->default_val(3875)
        ->check(Range(281, 10000));
    app.add_option("-S,--n-rounds", n_rounds, "Number of simulations")
        ->default_val(1)
        ->check(Range(1, 10000));
    app.add_option("--rate", rate, "Risk-free rate")
        ->default_val(0.01)
        ->check(Range(0.0, 2.0));
    app.add_option("--plot", plot, "Output plots")
        ->default_val(false);
    app.add_option("--type-neb", type_neb, "Agents' cognitive traits")
        ->default_val("Classic")
        ->check(IsMember({
            "Classic", "Algorithmic", "Human", "LossAversion", "Positivity", "Negativity", "DelayDiscounting", "Fear",
            "Greed", "LearningRate"
        }));
    app.add_option("--hp-gesture", hp_gesture, "Transaction gesture scalar")
        ->default_val(1.0)
        ->check(Range(1.0, 10.0));
    app.add_option("--liquidation-floor", liquidation_floor, "YTD drawdown")
        ->default_val(50)
        ->check(Range(0, 100));
    app.add_option("--leader-type", leader_type, "Type of agent leader other agents emulate")
        ->default_val("NoCluster")
        ->check(IsMember({"Worst", "Best", "Static", "Noise", "NoCluster"}));
    app.add_option("--cluster-limit", cluster_limit, "Percentage of agents imitating agent leader")
        ->default_val(1)
        ->check(Range(0, 100));
    app.add_flag("-v", verbosity, "Level of output verbosity");
    
    try {
        app.parse(argc, argv);
    }
    catch (const ParseError& e) {
        std::exit(app.exit(e));
    }

    // Log streams
    #if defined(SYMBA_OSAPI_WINDOWS)
    static ofstream cnull_("nul");
    #elif defined(SYMBA_OSAPI_POSIX)
    static ofstream cnull_("/dev/null");
    #else
    #error Unknown OS API.
    #endif

    cnull.rdbuf(cnull_.rdbuf());

    log0.rdbuf(cout.rdbuf());

    if (verbosity >= 1) log1.rdbuf(cout.rdbuf());
    else log1.rdbuf(cnull.rdbuf());
}

namespace terminal {

pair<int, int> size() {
    // I hate this
    #if defined(SYMBA_OSAPI_WINDOWS)
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    bool ok = GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);

    if (!ok) return pair(0, 0); // Probably no terminal

    int width = csbi.srWindow.Right - csbi.srWindow.Left + 1;
    int height = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;

    if (width == 1 and height == 1) return pair(0, 0); // No terminal
    else return pair(width, height);

    #elif defined(SYMBA_OSAPI_POSIX)
    struct winsize w;
    ioctl(fileno(stdout), TIOCGWINSZ, &w);
    int width = w.ws_col;
    int height = w.ws_row;

    return pair(width, height);

    #else
    #error Unknown OS API.
    #endif
}

int width() {
    return size().first;
}

int height() {
    return size().second;
}

} // namespace terminal

namespace args {

void dump(fs::path filename, bool skip_non_model_args) {
    json j = to_json();

    if (skip_non_model_args) {
        // Exclude non-model parameters
        j.erase("output-dir");
    }

    ofstream ofs(filename);
    ofs << setw(4) << j;
}

json to_json() {
    json j;

    j["output-dir"] = output_dir.string();
    j["n-agents"] = n_agents;
    j["n-stocks"] = n_stocks;
    j["n-steps"] = n_steps;
    j["n-rounds"] = n_rounds;
    j["rate"] = rate;
    j["plot"] = plot;
    j["type-neb"] = type_neb;
    j["hp-gesture"] = hp_gesture;
    j["liquidation-floor"] = liquidation_floor;
    j["leader-type"] = leader_type;
    j["cluster-limit"] = cluster_limit;

    return j;
}

} // namespace args
} // namespace cli
#include <CLI/CLI.hpp>
#include <fstream>
#include <iomanip>
#include "cli.hpp"

CliArgs::CliArgs(int argc, char** argv) {
    using namespace CLI;
    App app("A market simulator built in C++.", "symba");

    app.add_option("-o,--output-dir", output_dir, "Output directory for simulation results")
        ->required()
        ->check(ExistingDirectory);
    app.add_option("-I,--n-agents", n_agents, "Number of agents in the simulation")
        ->default_val(500)
        ->check(Range(0, 9999));
    app.add_option("-J,--n-stocks", n_stocks, "Number of stocks in the simulation")
        ->default_val(1)
        ->check(Range(0, 99));
    app.add_option("-T,--n-steps", n_steps, "Number of simulation time steps")
        ->default_val(3875)
        ->check(Range(282, 9999));
    app.add_option("-S,--n-rounds", n_rounds, "Number of simulations")
        ->default_val(1)
        ->check(Range(1, 9999));
    app.add_option("--rate", rate, "Risk-free rate")
        ->default_val(0.01)
        ->check(Range(0.01, 1.99));
    app.add_option("--plot", plot, "Output plots")
        ->default_val(false);
    app.add_option("--type-neb", type_neb, "Agents' cognitive traits")
        ->default_val("Classic")
        ->check(IsMember({
            "Classic", "Algorithmic", "Human", "LossAversion", "Positivity", "Negativity", "DelayDiscounting", "Fear",
            "Greed", "LearningRate"
        }));
    app.add_option("--hp-gesture", hp_gesture, "Transaction gesture scalar")
        ->default_val(0)
        ->check(Range(0, 9));
    app.add_option("--liquidation-floor", liquidation_floor, "YTD drawdown")
        ->default_val(50)
        ->check(Range(1, 99));
    app.add_option("--leader-type", leader_type, "Type of agent leader other agents emulate")
        ->default_val("NoCluster")
        ->check(IsMember({"Worst", "Best", "Static", "Noise", "NoCluster"}));
    app.add_option("--cluster-limit", cluster_limit, "Percentage of agents imitating agent leader")
        ->default_val(1)
        ->check(Range(1, 99));

    try {
        app.parse(argc, argv);
    }
    catch (const ParseError& e) {
        std::exit(app.exit(e));
    }
}

void CliArgs::dump(fs::path filename, bool skip_non_model_args) {
    json j = (*this);

    if (skip_non_model_args) {
        // Exclude non-model parameters
        j.erase("output-dir");
    }

    ofstream ofs(filename);
    ofs << setw(4) << j;
}

void to_json(json& j, const CliArgs& args) {
    j["output-dir"] = args.output_dir.string();
    j["n-agents"] = args.n_agents;
    j["n-stocks"] = args.n_stocks;
    j["n-steps"] = args.n_steps;
    j["n-rounds"] = args.n_rounds;
    j["rate"] = args.rate;
    j["plot"] = args.plot;
    j["type-neb"] = args.type_neb;
    j["hp-gesture"] = args.hp_gesture;
    j["liquidation-floor"] = args.liquidation_floor;
    j["leader-type"] = args.leader_type;
    j["cluster-limit"] = args.cluster_limit;
}
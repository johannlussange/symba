#include <CLI/CLI.hpp>
#include "cli.hpp"

CliArgs::CliArgs(int argc, char** argv) {
    using namespace CLI;
    App app("A market simulator built in C++.", "symba");

    app.add_option("-o,--output-dir", output_dir, "Output directory for simulation results")
        ->required()
        ->check(ExistingDirectory);
    app.add_option("-i,--n-agents", n_agents, "Number of agents in the simulations")
        ->default_val(500)
        ->check(NonNegativeNumber);
    app.add_option("-j,--n-stocks", n_stocks, "Number of stocks in the simulation")
        ->default_val(1)
        ->check(NonNegativeNumber);
    app.add_option("--time", time, "Amount of time steps")
        ->default_val(3875)
        ->check(NonNegativeNumber);
    app.add_option("--rate", rate, "No description")
        ->default_val(0.01)
        ->check(NonNegativeNumber);
    app.add_option("--plot", plot, "No description")
        ->default_val("Off")
        ->check(IsMember({"Off", "On"}));
    app.add_option("--pd-condition", pd_condition, "No description")
        ->default_val("PDOff")
        ->check(IsMember({"PDOff", "PDOn"}));
    app.add_option("--type-neb", type_neb, "No description")
        ->default_val("Algorithmic")
        ->check(IsMember({
            "Classic", "Algorithmic", "Human", "LossAversion", "Positivity", "Negativity", "DelayDiscounting", "Fear",
            "Greed", "LearningRate"
        }));
    app.add_option("--hp-gesture", hp_gesture, "No description")
        ->default_val(0)
        ->check(Range(0, 3));
    app.add_option("--hp-true-mu", hp_true_mu, "No description")
        ->default_val(0.1);
    app.add_option("--hp-accuracy", hp_accuracy, "No description")
        ->default_val(10);
    app.add_option("--liquidation-floor", liquidation_floor, "No description")
        ->default_val(50);
    app.add_option("--leader-type", leader_type, "No description")
        ->default_val("NoCluster")
        ->check(IsMember({"Worst", "Best", "Static", "Noise", "NoCluster"}));
    app.add_option("--cluster-limit", cluster_limit, "No description")
        ->default_val(1);
    app.add_option("-s", s, "No Description")
        ->default_val(0);

    try {
        app.parse(argc, argv);
    }
    catch (const ParseError& e) {
        std::exit(app.exit(e));
    }
}
#pragma once
#include <filesystem>
#include <string>

using namespace std;
namespace fs = std::filesystem;

/// Process command-line arguments.
struct CliArgs {
    fs::path output_dir;
    int n_agents;
    int n_stocks;
    int n_steps;
    int n_rounds;
    double rate;
    bool plot;
    string type_neb;
    int hp_gesture;
    int liquidation_floor;
    string leader_type;
    int cluster_limit;

    CliArgs(int argc, char** argv);
};
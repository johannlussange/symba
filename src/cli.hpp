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
    int time;
    double rate;
    string plot;
    string pd_condition;
    string type_neb;
    int hp_gesture;
    double hp_true_mu;
    int hp_accuracy;
    int liquidation_floor;
    string leader_type;
    int cluster_limit;
    int s;

    CliArgs(int argc, char** argv);
};
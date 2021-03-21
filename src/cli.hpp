#pragma once
#include <filesystem>
#include <string>
#include <string_view>
#include <nlohmann/json.hpp>

using namespace std;
namespace fs = std::filesystem;
using json = nlohmann::json;

/// Processed command-line arguments.
struct CliArgs {
    fs::path output_dir;
    int n_agents;
    int n_stocks;
    int n_steps;
    int n_rounds;
    double rate;
    bool plot;
    string type_neb;
    double hp_gesture;
    int liquidation_floor;
    string leader_type;
    int cluster_limit;

    CliArgs(int argc, char** argv);

    /**
     * @brief Write the contents of CliArgs to a file.
     * @param filename: path to the target file.
     * @param skip_non_model_args: skip arguments, unrelated to the SYMBA model (like output-dir).
    */
    void dump(fs::path filename, bool skip_non_model_args = true);
};

void to_json(json& j, const CliArgs& args);
#pragma once
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <nlohmann/json.hpp>

using namespace std;
namespace fs = std::filesystem;
using json = nlohmann::json;

namespace cli {

/// Initialize command-line interface. Call this in main() as soon as possible!
void init(int argc, char** argv);

inline int argc;
inline char** argv;

/// Terminal parameters
namespace terminal {
    /// Terminal size - columns x rows (0x0 for no terminal).
    pair<int, int> size();

    /// Number of columns in the terminal window (0 for no terminal).
    int width();

    /// Number of rows in the terminal window (0 for no terminal).
    int height();
};

/// Processed command-line arguments.
namespace args {
    inline fs::path output_dir;
    inline int n_jobs;
    inline int n_agents;
    inline int n_stocks;
    inline int n_steps;
    inline int n_rounds;
    inline double rate;
    inline bool plot;
    inline string type_neb;
    inline double hp_gesture;
    inline int liquidation_floor;
    inline string leader_type;
    inline int cluster_limit;
    inline int verbosity;

    /**
     * @brief Write the contents of args to a file.
     * @param filename: path to the target file.
     * @param skip_non_model_args: skip arguments, unrelated to the SYMBA model (like output-dir).
    */
    void dump(fs::path filename, bool skip_non_model_args = true);
    json to_json();
}

/// Output stream that discards contents. Rebind log streams to this using `log0.rdbuf(cnull.rdbuf())` to discard logs.
inline std::ostream cnull(nullptr);

/// Output stream for messages with a verbosity of 0 (i.e. always outputs). Equivalent to cout.
inline std::ostream log0(nullptr);

/// Output stream for messages with a verbosity of 1 (flag -v).
inline std::ostream log1(nullptr);

} // namespace cli
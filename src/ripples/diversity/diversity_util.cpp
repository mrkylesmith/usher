#include "diversity_util.hpp"
#include <array>
#include <iostream>
#include <tbb/task_scheduler_init.h>

std::optional<var_map> parse_diversity_flags(int argc, char *argv[]) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    uint32_t num_threads;

    std::string num_threads_message = "Number of threads to use when possible "
                                      "[DEFAULT uses all available cores, " +
                                      std::to_string(num_cores) +
                                      " detected on this machine]";

    po::options_description diversity_options("diversity options");
    diversity_options.add_options()(
        "input-mat,i", po::value<std::string>()->required(),
        "Input mutation-annotated tree file. [REQUIRED]")(
        "chronumental,C", po::value<std::string>()->required(),
        "Chronumental dates file corresponding to the given MAT. [REQUIRED]")(
        "month,m", po::value<std::string>()->required(),
        "Given month to calculate standing genetic diversity for. [REQUIRED]")(
        "output-directory,o", po::value<std::string>(),
        "Output directory to write all output files to [Optional]. ")(
        "threads,T",
        po::value<uint32_t>(&num_threads)->default_value(num_cores),
        num_threads_message.c_str())("help,h", "Print help message.");

    po::options_description all_options;
    all_options.add(diversity_options);
    po::positional_options_description p;
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(all_options)
                      .positional(p)
                      .run(),
                  vm);
        po::notify(vm);
    } catch (std::exception &e) {
        std::cerr << diversity_options << std::endl;
        // Either returning with error or help message
        return std::nullopt;
    }
    return vm;
}

std::unordered_map<std::string, std::vector<std::string>> init_month_dict() {
    std::array<std::string, 4> years = {"2020", "2021", "2022", "2023"};
    std::array<std::string, 12> months = {"01", "02", "03", "04", "05", "06",
                                          "07", "08", "09", "10", "11", "12"};
    using Samples = std::vector<std::string>;
    using Bins = std::unordered_map<std::string, Samples>;
    Bins bins;
    bins.reserve(years.size() * months.size());

    // Parse into format like, 2022-11
    for (const auto &y : years) {
        for (const auto &m : months) {
            bins.insert(std::make_pair(y + "-" + m, Samples()));
        }
    }
    return bins;
}

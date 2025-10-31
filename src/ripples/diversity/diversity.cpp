#include "diversity.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>

genetic_diversity::genetic_diversity(var_map &vm)
    : mat_filename_(vm["input-mat"].as<string>()),
      chronumental_filename_(vm["chronumental"].as<string>()),
      month_(vm["month"].as<string>()),
      num_threads_(vm["threads"].as<uint32_t>()), init_(num_threads_),
      tree_(load_mat()), parser_(chronumental_filename_) {}

MAT::Tree genetic_diversity::load_mat() {
    fprintf(stdout, "Loading input MAT file: %s.\n", mat_filename_.c_str());

    if (mat_filename_.find(".pb\0") == std::string::npos) {
        fprintf(stderr,
                "ERROR: Input file ending not recognized. Must be .pb\n");
        exit(1);
    }

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(mat_filename_);
    T.uncondense_leaves();
    return T;
}

genetic_diversity::Bins genetic_diversity::bin_samples_by_month() {
    // Initialize the month bins
    auto bins = init_month_dict();

    std::vector<std::string> samples;
    bool header = true;
    bool exclude_internal_nodes = true;
    int sample_column_num = 0;
    int date_column_num = 1;

    // If file has header, skip over first header line
    // Assuming header is just a single first line in TSV file
    if (header) {
        parser_.next_line();
    }

    for (; !parser_.done(); parser_.next_line()) {
        auto sample = parser_.get_value(sample_column_num);
        // Skip internal nodes
        if (exclude_internal_nodes) {
            if (sample.find("node_") != std::string_view::npos) {
                continue;
            }
        }
        auto date = parser_.get_value(date_column_num);
        if (date.empty() || date.size() < 10) {
            continue;
        }
        auto month = string{date.substr(0, 7)};
        // Check if sample data out of range
        if (bins.find(month) == bins.end()) {
            continue;
        }
        bins[month].emplace_back(sample);
    }

    // Check if month is valid
    if (bins.find(month_) == bins.end()) {
        const std::string message = "Erorr: Month: " + month_ +
                                    " not valid, please enter month "
                                    "between 2020-2023 (eg. '2020-08')";
        throw std::runtime_error(message);
    }
    return bins;
}

float genetic_diversity::compute_subtree_diversity(MAT::Tree &subtree) {
    // Phylogenetic entropy index for a given subtree
    float diversity = 0.0;

    subtree.uncondense_leaves();
    size_t subtree_num_leaves = subtree.get_num_leaves();
    for (const auto n : subtree.depth_first_expansion()) {
        float branch_length = static_cast<float>(n->branch_length);
        size_t num_descendants = subtree.get_num_leaves(n);
        if (branch_length <= 0 || num_descendants == 0) {
            continue;
        }

        float proportion = static_cast<float>(num_descendants) /
                           static_cast<float>(subtree_num_leaves);
        diversity += (branch_length * proportion * log(proportion));
    }
    return -diversity;
}

float genetic_diversity::calculate(const Bins &bins) {
    // Get samples for the given month
    const auto &samples = bins.at(month_);
    if (samples.size() == 0) {
        const std::string message =
            "Erorr: Zero samples found in given month subtree: " + month_;
        throw std::runtime_error(message);
    }

    std::cout << "Computing standing genetic diversity for month: " << month_
              << ", with " << samples.size() << " samples.\n";

    // Build a subtree containing just the samples for the given month
    MAT::Tree subtree = get_subtree(tree_, samples, true);

    float diversity_score = compute_subtree_diversity(subtree);
    return diversity_score;
}


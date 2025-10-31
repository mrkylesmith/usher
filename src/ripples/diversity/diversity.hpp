#pragma once

#include "diversity_util.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "src/ripples/util/text_parser.hpp"
#include "src/usher_graph.hpp"

namespace MAT = Mutation_Annotated_Tree;

class genetic_diversity {
  public:
    using string = std::string;
    using Samples = std::vector<string>;
    using Bins = std::unordered_map<string, Samples>;

    genetic_diversity(var_map &vm);
    genetic_diversity(const genetic_diversity &) = delete;
    genetic_diversity &operator=(const genetic_diversity &) = delete;

    std::string_view month() const noexcept { return month_; }
    Bins bin_samples_by_month();
    float calculate(const Bins &bins);

  private:
    MAT::Tree load_mat();
    void fill_bins();
    float compute_subtree_diversity(MAT::Tree &subtree);

    string mat_filename_;
    string chronumental_filename_;
    string month_;
    uint32_t num_threads_;
    tbb::task_scheduler_init init_;
    MAT::Tree tree_;
    text_parser parser_;
};

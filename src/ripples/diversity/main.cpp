#include "diversity.hpp"
#include "diversity_util.hpp"
#include <iostream>

// Program to calculate the standing genetic diversity of a given month subtree.
int main(int argc, char *argv[]) {
    auto vm = parse_diversity_flags(argc, argv);
    if (!vm) {
        return 1;
    }

    // Calculate the standing genetic diversity for the given month subtree
    genetic_diversity diversity(vm.value());
    auto bins = diversity.bin_samples_by_month();
    float score = diversity.calculate(bins);
    std::cout << "Month: " << diversity.month() << ", Diversity: " << score
              << '\n';
    return 0;
}

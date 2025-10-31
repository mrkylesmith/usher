#pragma once

#include <array>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <optional>

namespace po = boost::program_options;
using var_map = boost::program_options::variables_map;
using parsed_options = boost::program_options::parsed_options;

std::optional<var_map> parse_diversity_flags(int argc, char *argv[]);

std::unordered_map<std::string, std::vector<std::string>> init_month_dict();

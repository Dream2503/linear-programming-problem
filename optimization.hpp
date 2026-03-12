#pragma once
#include <iomanip>
#include <map>
#include "linear-algebra/linalg.hpp"

namespace optimization {
    inline struct FormatSettings {
        bool verbose = false;
        std::ostream* out = &std::cout;
    } GLOBAL_FORMATTING;

    enum class Optimization : bool { MINIMIZE, MAXIMIZE };
    enum class Solution : uint8_t { UNOPTIMIZED, OPTIMIZED, INFEASIBLE, UNBOUNDED, ALTERNATE };
    class LPP;
    class ComputationalTable;

    std::vector<std::map<algebra::Variable, algebra::Fraction>> basic_feasible_solutions(const std::vector<algebra::Equation>&);
} // namespace lpp

#include "src/lpp.hpp"
#include "src/computation_table.hpp"

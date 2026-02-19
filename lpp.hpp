#pragma once
#include <iomanip>
#include <map>
#include <utility>
#include "algebra/algebra.hpp"

namespace lpp {
    using namespace algebra;

    enum class Optimization { MAXIMIZE, MINIMIZE };
    enum class Solution { UNOPTIMIZED, OPTIMIZED, INFEASIBLE, UNBOUNDED, ALTERNATE };

    class LPP;
    class ComputationalTable;

    std::vector<std::map<Variable, Fraction>> basic_feasible_solutions(const std::vector<Equation>&);

    namespace detail {
        std::vector<std::vector<int>> generate_combination(int, int);
    } // namespace detail
} // namespace lpp

#include "detail.hpp"
#include "src/lpp.hpp"

#pragma once
#include <map>
#include "algebra/algebra.hpp"

namespace lpp {
    using namespace algebra;
    class LPP;

    std::vector<std::map<Variable, Fraction>> basic_feasible_solutions(const std::vector<Equation>&);

    namespace detail {
        std::vector<std::vector<int>> generate_combination(int, int);
    }
} // namespace lpp

#include "detail.hpp"
#include "src/lpp.hpp"

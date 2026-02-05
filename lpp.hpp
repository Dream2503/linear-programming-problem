#pragma once
#include <assert.h>
#include <cmath>
#include <iostream>
#include <ranges>
#include <vector>

namespace lpp {
    enum class OptimizeType { MAXIMIZE, MINIMIZE };

    class Variable;
    class Polynomial;
    class Inequation;
    class Equation;
    class LPP;

    double optimize(LPP);
} // namespace lpp

#include "types.hpp"
#include "simplex.hpp"

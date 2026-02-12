#pragma once
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <numeric>
#include <ranges>
#include <vector>

namespace algebra {
    class Fraction;
    class Variable;
    class Polynomial;
    class Inequation;
    class Equation;

    template <typename>
    class Matrix;
} // namespace algebra

#include "src/fraction.hpp"
#include "src/variable.hpp"
#include "src/polynomial.hpp"
#include "src/inequation.hpp"
#include "src/matrix.hpp"

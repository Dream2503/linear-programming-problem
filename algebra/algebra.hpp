#pragma once
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <numeric>
#include <ranges>
#include <vector>

namespace algebra {
    enum class RelationalOperator { LT, LE, GT, GE, EQ };

    class Fraction;
    class Variable;
    class Polynomial;
    class Inequation;
    class Equation;
    class Interval;

    template <typename>
    class Matrix;

    namespace detail {
        RelationalOperator invert_relational_operator(RelationalOperator);
    }

    inline std::ostream& operator<<(std::ostream& out, const RelationalOperator opr) {
        switch (opr) {
        case RelationalOperator::LT:
            out << '<';
            break;

        case RelationalOperator::LE:
            out << "<=";
            break;

        case RelationalOperator::GT:
            out << '>';
            break;

        case RelationalOperator::GE:
            out << ">=";
            break;

        case RelationalOperator::EQ:
            out << '=';
        }
        return out;
    }
} // namespace algebra

#include "src/detail.hpp"
#include "src/fraction.hpp"
#include "src/variable.hpp"
#include "src/polynomial.hpp"
#include "src/inequation.hpp"
#include "src/interval.hpp"
#include "src/matrix.hpp"

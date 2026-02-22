#pragma once
#include "fraction.hpp"

namespace algebra::detail {
    inline RelationalOperator invert_relational_operator(const RelationalOperator opr) {
        switch (opr) {
        case RelationalOperator::LT:
            return RelationalOperator::GT;

        case RelationalOperator::LE:
            return RelationalOperator::GE;

        case RelationalOperator::GT:
            return RelationalOperator::LT;

        case RelationalOperator::GE:
            return RelationalOperator::LE;

        default:
            return RelationalOperator::EQ;
        }
    }

    inline bool evaluate_relational_operator(const Fraction& lhs, const RelationalOperator opr, const Fraction& rhs) {
        switch (opr) {
        case RelationalOperator::LT:
            return lhs < rhs;

        case RelationalOperator::LE:
            return lhs <= rhs;

        case RelationalOperator::GT:
            return lhs > rhs;

        case RelationalOperator::GE:
            return lhs >= rhs;

        default:
            return lhs == rhs;
        }
    }
} // namespace algebra::detail

#pragma once

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
} // namespace algebra::detail

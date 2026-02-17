#pragma once

class algebra::Interval {
public:
    Inequation lhs;
    RelationalOperator opr;
    Polynomial rhs;

    constexpr Interval() = default;

    Interval(const Inequation& lhs, const RelationalOperator opr, const Polynomial& rhs) : lhs(lhs), opr(opr), rhs(rhs) {}

    friend std::ostream& operator<<(std::ostream& out, const Interval& interval) {
        return out << interval.lhs << ' ' << interval.opr << interval.rhs;
    }
};

inline algebra::Interval operator<(const algebra::Inequation& lhs, const algebra::Fraction& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Interval operator<=(const algebra::Inequation& lhs, const algebra::Fraction& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Interval operator>(const algebra::Inequation& lhs, const algebra::Fraction& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Interval operator>=(const algebra::Inequation& lhs, const algebra::Fraction& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Interval operator==(const algebra::Inequation& lhs, const algebra::Fraction& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::EQ, rhs);
}

inline algebra::Interval operator<(const algebra::Fraction& lhs, const algebra::Inequation& rhs) { return rhs > lhs; }

inline algebra::Interval operator<=(const algebra::Fraction& lhs, const algebra::Inequation& rhs) { return rhs >= lhs; }

inline algebra::Interval operator>(const algebra::Fraction& lhs, const algebra::Inequation& rhs) { return rhs < lhs; }

inline algebra::Interval operator>=(const algebra::Fraction& lhs, const algebra::Inequation& rhs) { return rhs <= lhs; }

inline algebra::Interval operator==(const algebra::Fraction& lhs, const algebra::Inequation& rhs) { return rhs == lhs; }

inline algebra::Interval operator<(const algebra::Inequation& lhs, const algebra::Variable& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Interval operator<=(const algebra::Inequation& lhs, const algebra::Variable& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Interval operator>(const algebra::Inequation& lhs, const algebra::Variable& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Interval operator>=(const algebra::Inequation& lhs, const algebra::Variable& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Interval operator==(const algebra::Inequation& lhs, const algebra::Variable& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::EQ, rhs);
}

inline algebra::Interval operator<(const algebra::Variable& lhs, const algebra::Inequation& rhs) { return rhs > lhs; }

inline algebra::Interval operator<=(const algebra::Variable& lhs, const algebra::Inequation& rhs) { return rhs >= lhs; }

inline algebra::Interval operator>(const algebra::Variable& lhs, const algebra::Inequation& rhs) { return rhs < lhs; }

inline algebra::Interval operator>=(const algebra::Variable& lhs, const algebra::Inequation& rhs) { return rhs <= lhs; }

inline algebra::Interval operator==(const algebra::Variable& lhs, const algebra::Inequation& rhs) { return rhs == lhs; }

inline algebra::Interval operator<(const algebra::Inequation& lhs, const algebra::Polynomial& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Interval operator<=(const algebra::Inequation& lhs, const algebra::Polynomial& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Interval operator>(const algebra::Inequation& lhs, const algebra::Polynomial& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Interval operator>=(const algebra::Inequation& lhs, const algebra::Polynomial& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Interval operator==(const algebra::Inequation& lhs, const algebra::Polynomial& rhs) {
    return algebra::Interval(lhs, algebra::RelationalOperator::EQ, rhs);
}

inline algebra::Interval operator<(const algebra::Polynomial& lhs, const algebra::Inequation& rhs) { return rhs > lhs; }

inline algebra::Interval operator<=(const algebra::Polynomial& lhs, const algebra::Inequation& rhs) { return rhs >= lhs; }

inline algebra::Interval operator>(const algebra::Polynomial& lhs, const algebra::Inequation& rhs) { return rhs < lhs; }

inline algebra::Interval operator>=(const algebra::Polynomial& lhs, const algebra::Inequation& rhs) { return rhs <= lhs; }

inline algebra::Interval operator==(const algebra::Polynomial& lhs, const algebra::Inequation& rhs) { return rhs == lhs; }

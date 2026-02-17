#pragma once

class algebra::Inequation {

public:
    RelationalOperator opr;
    Polynomial lhs, rhs;

    constexpr Inequation() = default;

    Inequation(const Polynomial& polynomial, const RelationalOperator opr, const Polynomial& rhs) : opr(opr), lhs(polynomial), rhs(rhs) {}

    Inequation& operator*=(const Fraction& value) {
        lhs *= value;
        rhs *= value;

        if (value < 0) {
            opr = detail::invert_relational_operator(opr);
        }
        return *this;
    }

    Inequation operator*(const Fraction& value) const {
        Inequation inequation = *this;
        inequation *= value;
        return inequation;
    }

    Inequation& operator*=(const Variable& value) {
        lhs *= value;
        rhs *= value;

        if (value.coefficient < 0) {
            opr = detail::invert_relational_operator(opr);
        }
        return *this;
    }

    Inequation operator*(const Variable& value) const {
        Inequation inequation = *this;
        inequation *= value;
        return inequation;
    }

    Inequation& operator*=(const Polynomial& value) {
        lhs *= value;
        rhs *= value;
        return *this;
    }

    Inequation operator*(const Polynomial& value) const {
        Inequation inequation = *this;
        inequation *= value;
        return inequation;
    }

    Inequation& operator/=(const Fraction& value) {
        lhs /= value;
        rhs /= value;

        if (value < 0) {
            opr = detail::invert_relational_operator(opr);
        }
        return *this;
    }

    Inequation operator/(const Fraction& value) const {
        Inequation inequation = *this;
        inequation /= value;
        return inequation;
    }

    Inequation& operator/=(const Variable& value) {
        lhs /= value;
        rhs /= value;

        if (value.coefficient < 0) {
            opr = detail::invert_relational_operator(opr);
        }
        return *this;
    }

    Inequation operator/(const Variable& value) const {
        Inequation inequation = *this;
        inequation /= value;
        return inequation;
    }

    friend std::ostream& operator<<(std::ostream& out, const Inequation& inequation) {
        return out << inequation.lhs << ' ' << inequation.opr << ' ' << inequation.rhs;
    }
};

class algebra::Equation : public Inequation {
public:
    Equation(const Polynomial& polynomial, const Polynomial& rhs) : Inequation(polynomial, RelationalOperator::EQ, rhs) {}
};

inline algebra::Inequation operator<(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Inequation operator<=(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Inequation operator>(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Inequation operator>=(const algebra::Variable& lhs, const algebra::Fraction& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Equation operator==(const algebra::Variable& lhs, const algebra::Fraction& rhs) { return algebra::Equation(lhs, rhs); }

inline algebra::Inequation operator<(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return rhs > lhs; }

inline algebra::Inequation operator<=(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return rhs >= lhs; }

inline algebra::Inequation operator>(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return rhs < lhs; }

inline algebra::Inequation operator>=(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return rhs <= lhs; }

inline algebra::Equation operator==(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return rhs == lhs; }

inline algebra::Inequation operator<(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Inequation operator<=(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Inequation operator>(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Inequation operator>=(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Equation operator==(const algebra::Variable& lhs, const algebra::Variable& rhs) { return algebra::Equation(lhs, rhs); }

inline algebra::Inequation operator<(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Inequation operator<=(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Inequation operator>(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Inequation operator>=(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Equation operator==(const algebra::Polynomial& lhs, const algebra::Fraction& rhs) { return algebra::Equation(lhs, rhs); }

inline algebra::Inequation operator<(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return rhs > lhs; }

inline algebra::Inequation operator<=(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return rhs >= lhs; }

inline algebra::Inequation operator>(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return rhs < lhs; }

inline algebra::Inequation operator>=(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return rhs <= lhs; }

inline algebra::Equation operator==(const algebra::Fraction& lhs, const algebra::Polynomial& rhs) { return rhs == lhs; }

inline algebra::Inequation operator<(const algebra::Polynomial& lhs, const algebra::Variable& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Inequation operator<=(const algebra::Polynomial& lhs, const algebra::Variable& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Inequation operator>(const algebra::Polynomial& lhs, const algebra::Variable& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Inequation operator>=(const algebra::Polynomial& lhs, const algebra::Variable& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Equation operator==(const algebra::Polynomial& lhs, const algebra::Variable& rhs) { return algebra::Equation(lhs, rhs); }

inline algebra::Inequation operator<(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return rhs > lhs; }

inline algebra::Inequation operator<=(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return rhs >= lhs; }

inline algebra::Inequation operator>(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return rhs < lhs; }

inline algebra::Inequation operator>=(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return rhs <= lhs; }

inline algebra::Equation operator==(const algebra::Variable& lhs, const algebra::Polynomial& rhs) { return rhs == lhs; }

inline algebra::Inequation operator<(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LT, rhs);
}

inline algebra::Inequation operator<=(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::LE, rhs);
}

inline algebra::Inequation operator>(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GT, rhs);
}

inline algebra::Inequation operator>=(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) {
    return algebra::Inequation(lhs, algebra::RelationalOperator::GE, rhs);
}

inline algebra::Equation operator==(const algebra::Polynomial& lhs, const algebra::Polynomial& rhs) { return algebra::Equation(lhs, rhs); }

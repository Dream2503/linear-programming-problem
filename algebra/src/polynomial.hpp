#pragma once

class algebra::Polynomial {
public:
    std::vector<Variable> expression;

    constexpr Polynomial() = default;

    Polynomial(const Variable& variable) { expression.push_back(variable); }
    Polynomial(const std::vector<Variable>& terms) : expression(terms) { std::ranges::sort(expression); }

    Polynomial operator-() const {
        Polynomial res = *this;

        for (Variable& variable : res.expression) {
            variable *= -1;
        }
        return res;
    }

    Polynomial& operator+=(const Polynomial& value) {
        for (const Variable& variable : value.expression) {
            *this += variable;
        }
        return *this;
    }

    Polynomial operator+(const Polynomial& value) const {
        Polynomial res = *this;
        res += value;
        return res;
    }

    Polynomial& operator+=(const Variable& value) {
        const auto itr = std::ranges::lower_bound(expression, value, [](const Variable& lhs, const Variable& rhs) -> bool {
            const bool lhs_const = lhs.variables.empty();
            const bool rhs_const = rhs.variables.empty();

            if (lhs_const != rhs_const) {
                return !lhs_const;
            }
            return lhs.variables < rhs.variables;
        });

        if (itr != expression.end() && itr->variables == value.variables) {
            itr->coefficient += value.coefficient;
        } else {
            expression.insert(itr, value);
        }
        return *this;
    }

    Polynomial operator+(const Variable& value) const {
        Polynomial polynomial = *this;
        polynomial += value;
        return polynomial;
    }

    Polynomial& operator+=(const Fraction& value) {
        *this += Variable(value);
        return *this;
    }

    Polynomial operator+(const Fraction& value) const { return *this + Variable(value); }

    Polynomial& operator-=(const Polynomial& value) {
        *this += -value;
        return *this;
    }

    Polynomial operator-(const Polynomial& value) const {
        Polynomial res = *this;
        res -= value;
        return res;
    }

    Polynomial& operator-=(const Variable& value) {
        *this += -value;
        return *this;
    }

    Polynomial operator-(const Variable& value) const { return *this + -value; }

    Polynomial operator-=(const Fraction& value) {
        *this += -value;
        return *this;
    }

    Polynomial operator-(const Fraction& value) const { return *this + -value; }

    Polynomial& operator*=(const Polynomial& value) {
        for (const Variable& variable : value.expression) {
            *this *= variable;
        }
        return *this;
    }

    Polynomial operator*(const Polynomial& value) const {
        Polynomial res = *this;
        res *= value;
        return res;
    }

    Polynomial& operator*=(const Variable& value) {
        for (Variable& variable : expression) {
            variable *= value;
        }
        return *this;
    }

    Polynomial operator*(const Variable& value) const {
        Polynomial res = *this;
        res *= value;
        return res;
    }

    Polynomial& operator*=(const Fraction& value) {
        for (Variable& variable : expression) {
            variable *= value;
        }
        return *this;
    }

    Polynomial operator*(const Fraction& value) const {
        Polynomial polynomial = *this;
        polynomial *= value;
        return polynomial;
    }

    Polynomial partial_substitution();

    Fraction complete_substitution();

    explicit operator Fraction() const {
        assert(expression.size() == 1);
        return static_cast<Fraction>(expression[0]);
    }

    friend std::ostream& operator<<(std::ostream& out, const Polynomial& polynomial) {
        if (!polynomial.expression.empty()) {
            out << polynomial.expression.front();

            for (const Variable& variable : polynomial.expression | std::views::drop(1)) {
                out << (variable.coefficient < Fraction() ? " - " : " + ") << std::abs(variable);
            }
        }
        return out;
    }
};

inline algebra::Polynomial operator+(const algebra::Variable& lhs, const algebra::Variable& rhs) { return algebra::Polynomial(lhs) + rhs; }

inline algebra::Polynomial operator+(const algebra::Variable& lhs, const algebra::Fraction& rhs) { return algebra::Polynomial(lhs) + rhs; }

inline algebra::Polynomial operator+(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return algebra::Polynomial(lhs) + rhs; }

inline algebra::Polynomial operator-(const algebra::Variable& lhs, const algebra::Variable& rhs) { return algebra::Polynomial(lhs) - rhs; }

inline algebra::Polynomial operator-(const algebra::Variable& lhs, const algebra::Fraction& rhs) { return algebra::Polynomial(lhs) - rhs; }

inline algebra::Polynomial operator-(const algebra::Fraction& lhs, const algebra::Variable& rhs) { return algebra::Polynomial(lhs) - rhs; }

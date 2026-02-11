#pragma once

class algebra::Polynomial {
public:
    std::vector<Variable> expression;

    constexpr Polynomial() = default;

    Polynomial(const std::vector<Variable>& terms) : expression(terms) { std::ranges::sort(expression); }

    Polynomial& operator+=(const Variable& value) {
        const auto itr = std::ranges::lower_bound(expression, value, [](const Variable& lhs, const Variable& rhs) -> bool {
            return std::tie(lhs.name, lhs.exponent) < std::tie(rhs.name, rhs.exponent);
        });

        if (itr != expression.end() && itr->name == value.name && itr->exponent == value.exponent) {
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

    Polynomial& operator-=(const Variable& value) {
        *this += -value;
        return *this;
    }

    Polynomial operator-(const Variable& value) const {
        Polynomial polynomial = *this;
        polynomial -= value;
        return polynomial;
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

    constexpr bool is_fraction() const { return expression.size() == 1 && expression[0].is_fraction(); }

    constexpr explicit operator Fraction() const {
        assert(is_fraction());
        return static_cast<Fraction>(expression[0]);
    }

    friend std::ostream& operator<<(std::ostream& out, const Polynomial& polynomial) {
        if (!polynomial.expression.empty()) {
            out << polynomial.expression.front();

            for (const Variable& variable : polynomial.expression | std::views::drop(1)) {
                out << (variable.coefficient < 0 ? " - " : " + ") << std::abs(variable);
            }
        }
        return out;
    }
};

inline algebra::Polynomial operator+(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Polynomial(std::vector{lhs, rhs});
}

inline algebra::Polynomial operator-(const algebra::Variable& lhs, const algebra::Variable& rhs) {
    return algebra::Polynomial(std::vector{lhs, -rhs});
}

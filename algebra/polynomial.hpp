#pragma once

class lpp::Polynomial {
public:
    std::vector<Variable> expression;

    Polynomial(const std::vector<Variable>& terms) : expression(terms) { std::ranges::sort(expression); }

    Polynomial& operator+=(const Variable& value) {
        const auto itr = std::ranges::lower_bound(expression, value);

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

inline lpp::Polynomial operator+(const lpp::Variable& lhs, const lpp::Variable& rhs) { return lpp::Polynomial(std::vector{lhs, rhs}); }

inline lpp::Polynomial operator-(const lpp::Variable& lhs, const lpp::Variable& rhs) { return lpp::Polynomial(std::vector{lhs, -rhs}); }

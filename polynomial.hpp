#pragma once

class Polynomial {
public:
    std::vector<Variable> terms;

    Polynomial(const std::vector<Variable>& terms) : terms(terms) {}

    friend std::ostream& operator<<(std::ostream& out, const Polynomial& polynomial) {
        if (!polynomial.terms.empty()) {
            out << polynomial.terms.front();

            for (const Variable& variable : polynomial.terms | std::views::drop(1)) {
                out << (variable.coefficient < 0 ? " - " : " + ") << std::abs(variable);
            }
        }
        return out;
    }
};

inline Polynomial operator+(const Variable& lhs, const Variable& rhs) {
    return Polynomial(std::vector{lhs, rhs});
}

inline Polynomial operator-(const Variable& lhs, const Variable& rhs) {
    return Polynomial(std::vector{lhs, -rhs});
}
#pragma once

class algebra::Variable {
    struct Var {
        std::string name;
        Fraction exponent;

        constexpr std::strong_ordering operator<=>(const Var&) const = default;
    };

    void print_variables(std::ostream& out) const {
        for (const auto& [name, exponent] : variables) {
            if (exponent != 1) {
                out << '(';
            }
            out << name;

            if (exponent != 1) {
                out << '^';

                if (exponent.denominator != 1) {
                    out << '(' << exponent << ')';
                } else {
                    out << exponent;
                }
            }
            if (exponent != 1) {
                out << ')';
            }
        }
    }

public:
    Fraction coefficient = 1;
    std::vector<Var> variables;

    Variable() = default;

    Variable(const std::string& name) : variables({{name, 1}}) {}

    Variable(const Fraction& coefficient) : coefficient(coefficient) {}

    Variable operator-() const {
        Variable res = *this;
        res.coefficient = -res.coefficient;
        return res;
    }

    Variable& operator*=(const Fraction& value) {
        this->coefficient *= value;
        return *this;
    }

    Variable operator*(const Fraction& value) const {
        Variable res = *this;
        res *= value;
        return res;
    }

    Variable& operator*=(const Variable& value) {
        coefficient *= value.coefficient;

        for (const Var& var : value.variables) {
            const auto itr = std::ranges::lower_bound(variables, var);

            if (itr != variables.end() && itr->name == var.name) {
                if ((itr->exponent += var.exponent) == 0) {
                    variables.erase(itr);
                }
            } else {
                variables.insert(itr, var);
            }
        }
        return *this;
    }

    Variable operator*(const Variable& value) const {
        Variable res = *this;
        res *= value;
        return res;
    }

    Variable& operator/=(const Fraction& value) {
        this->coefficient /= value;
        return *this;
    }

    Variable operator/(const Fraction& value) const {
        Variable res = *this;
        res /= value;
        return res;
    }

    Variable& operator/=(const Variable& value) {
        coefficient /= value.coefficient;

        for (const Var& var : value.variables) {
            const auto itr = std::ranges::lower_bound(variables, var);

            if (itr != variables.end() && itr->name == var.name) {
                if ((itr->exponent -= var.exponent) == 0) {
                    variables.erase(itr);
                }
            } else {
                variables.emplace(itr, var.name, -var.exponent);
            }
        }
        return *this;
    }

    Variable operator/(const Variable& value) const {
        Variable res = *this;
        res /= value;
        return res;
    }

    Variable& operator^=(const Fraction& value) {
        coefficient ^= value;

        for (auto& [_, exponent] : variables) {
            exponent *= value;
        }
        return *this;
    }

    Variable operator^(const Fraction& value) const {
        Variable variable = *this;
        variable ^= value;
        return variable;
    }

    constexpr std::strong_ordering operator<=>(const Variable& value) const {
        const bool is_const = this->variables.empty();
        const bool value_const = value.variables.empty();

        if (is_const != value_const) {
            return is_const ? std::strong_ordering::greater : std::strong_ordering::less;
        }
        return std::tie(variables, coefficient) <=> std::tie(value.variables, value.coefficient);
    }

    constexpr bool operator==(const Variable&) const = default;

    Variable partial_substitution();

    Fraction complete_substitution();

    Variable basis() const {
        Variable res = *this;
        res.coefficient = 1;

        for (auto& [_, exponent] : res.variables) {
            exponent = 1;
        }
        return res;
    }

    constexpr explicit operator Fraction() const {
        assert(variables.empty());
        return coefficient;
    }

    friend std::ostream& operator<<(std::ostream& out, const Variable& variable) {
        if (variable.variables.size() == 0) {
            out << variable.coefficient;
            return out;
        }
        if (variable.coefficient == 0) {
            out << '0';
        } else if (variable.coefficient == 1) {
            variable.print_variables(out);
        } else if (variable.coefficient == -1) {
            out << '-';
            variable.print_variables(out);
        } else if (variable.coefficient.denominator != 1) {
            if (variable.coefficient.numerator == 1) {
                variable.print_variables(out);
                out << '/' << variable.coefficient.denominator;
            } else {
                out << '(' << variable.coefficient.numerator;
                variable.print_variables(out);
                out << '/' << variable.coefficient.denominator << ')';
            }
        } else {
            out << variable.coefficient;
            variable.print_variables(out);
        }
        return out;
    }
};

namespace std {
    inline algebra::Variable abs(algebra::Variable variable) {
        variable.coefficient = std::abs(variable.coefficient);
        return variable;
    }
} // namespace std

constexpr algebra::Variable operator*(const algebra::Fraction& value, const algebra::Variable& variable) { return variable * value; }

constexpr algebra::Variable operator/(const algebra::Fraction& value, const algebra::Variable& variable) { return variable / value; }

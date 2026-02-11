#pragma once

class algebra::Variable {
public:
    static constexpr std::string CONSTANT = "CONSTANT";
    std::string name;
    Fraction coefficient = 1, exponent = 1;

    constexpr Variable() = default;

    constexpr Variable(const std::string& name) : name(name) {}

    constexpr Variable(const Fraction value) : name(CONSTANT), coefficient(value) {}

    constexpr Variable operator-() const {
        Variable variable = *this;
        variable.coefficient = -variable.coefficient;
        return variable;
    }

    constexpr Variable& operator*=(const Fraction& value) {
        coefficient *= value;
        return *this;
    }

    constexpr Variable operator*(const Fraction& value) const {
        Variable variable = *this;
        variable *= value;
        return variable;
    }

    constexpr Variable& operator^=(const Fraction& value) {
        coefficient ^= value;

        if (name != CONSTANT) {
            exponent *= value;
        }
        return *this;
    }

    constexpr Variable operator^(const Fraction& value) const {
        Variable variable = *this;
        variable ^= value;
        return variable;
    }

    constexpr Variable& operator*=(const Variable& variable) {
        if (name == CONSTANT) {
            name = variable.name;
            exponent = variable.exponent;
        } else {
            name += variable.name;
            exponent += variable.exponent;
        }
        coefficient *= variable.coefficient;
        return *this;
    }

    constexpr Variable operator*(const Variable& value) const {
        Variable variable = *this;
        variable *= value;
        return variable;
    }

    constexpr std::strong_ordering operator<=>(const Variable& other) const {
        return std::tie(name, exponent, coefficient) <=> std::tie(other.name, other.exponent, other.coefficient);
    }

    constexpr bool operator==(const Variable&) const = default;

    constexpr Fraction substitute(const Fraction value) const { return name != CONSTANT ? coefficient * (value ^ exponent) : Fraction(); }

    constexpr Variable basis() const { return Variable(name); }

    constexpr bool is_fraction() const { return name == CONSTANT; }

    constexpr explicit operator Fraction() const {
        assert(is_fraction());
        return coefficient;
    }

    constexpr friend Variable operator*(const Fraction& value, const Variable& variable) { return variable * value; }

    friend std::ostream& operator<<(std::ostream& out, const Variable& variable) {
        if (variable.exponent == 0 || variable.name == CONSTANT) {
            out << variable.coefficient;
            return out;
        }
        if (variable.coefficient == 0) {
            out << '0';
        } else if (variable.coefficient == 1) {
            out << variable.name;
        } else {
            if (variable.coefficient == -1) {
                out << '-' << variable.name;
            } else if (variable.coefficient.denominator != 1) {
                out << '(' << variable.coefficient.numerator << variable.name << '/' << variable.coefficient.denominator << ')';
            } else {
                out << variable.coefficient.numerator << variable.name;
            }
            if (variable.exponent != 1) {
                out << '^';

                if (variable.exponent.denominator != 1) {
                    out << '(' << variable.exponent << ')';
                } else {
                    out << variable.exponent;
                }
            }
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

#pragma once

class Variable {
public:
    std::string name;
    double coefficient, exponent;

    explicit Variable(const std::string& name, const double coefficient = 1, const double exponent = 1) :
        name(name), coefficient(coefficient), exponent(exponent) {}

    Variable operator-() const {
        return Variable(name, -coefficient, exponent);
    }

    Variable operator*(const double value) const { return Variable(name, coefficient * value, exponent); }

    Variable& operator*=(const double value) {
        coefficient *= value;
        return *this;
    }

    Variable operator^(const double value) const { return Variable(name, std::pow(coefficient, value), exponent * value); }

    Variable& operator^=(const double value) {
        coefficient = std::pow(coefficient, value);
        exponent *= value;
        return *this;
    }

    Variable operator*(const Variable& variable) const {
        return Variable(name + variable.name, coefficient * variable.coefficient, exponent + variable.exponent);
    }

    Variable& operator*=(const Variable& variable) {
        name += variable.name;
        coefficient *= variable.coefficient;
        exponent += variable.exponent;
        return *this;
    }

    friend Variable operator*(const double value, const Variable& variable) { return variable * value; }

    friend std::ostream& operator<<(std::ostream& out, const Variable& variable) {
        if (variable.exponent == 0) {
            out << variable.coefficient;
        } else {
            if (variable.coefficient == 0) {
                out << '0';
            } else if (variable.coefficient == 1) {
                out << variable.name;
            } else {
                if (variable.coefficient == -1) {
                    out << '-';
                } else {
                    out << variable.coefficient;
                }
                out << variable.name;

                if (variable.exponent != 1) {
                    out << '^' << variable.exponent;
                }
            }
        }
        return out;
    }
};

namespace std {
    inline Variable abs(const Variable& variable) { return Variable(variable.name, abs(variable.coefficient), variable.exponent); }
} // namespace std

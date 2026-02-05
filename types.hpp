#pragma once

class lpp::Variable {
public:
    static constexpr std::string CONSTANT = "CONSTANT";
    std::string name;
    double coefficient = 1, exponent = 1;

    Variable(const std::string& name) : name(name) {}
    Variable(const double value) : name(CONSTANT), coefficient(value) {}

    Variable operator-() const {
        Variable variable = *this;
        variable.coefficient = -variable.coefficient;
        return variable;
    }

    Variable& operator*=(const double value) {
        coefficient *= value;
        return *this;
    }

    Variable operator*(const double value) const {
        Variable variable = *this;
        variable *= value;
        return variable;
    }

    Variable& operator^=(const double value) {
        coefficient = std::pow(coefficient, value);

        if (name != CONSTANT) {
            exponent *= value;
        }
        return *this;
    }

    Variable operator^(const double value) const {
        Variable variable = *this;
        variable ^= value;
        return variable;
    }

    Variable& operator*=(const Variable& variable) {
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

    Variable operator*(const Variable& value) const {
        Variable variable = *this;
        variable *= value;
        return variable;
    }

    double substitute(const double value) const { return name != CONSTANT ? coefficient * std::pow(value, exponent) : 0; }

    explicit operator double() const {
        assert(name == CONSTANT);
        return coefficient;
    }

    friend Variable operator*(const double value, const Variable& variable) { return variable * value; }

    friend std::ostream& operator<<(std::ostream& out, const Variable& variable) {
        if (variable.exponent == 0 || variable.name == CONSTANT) {
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
    inline lpp::Variable abs(const lpp::Variable& variable) {
        lpp::Variable res = variable;
        res.coefficient = std::abs(res.coefficient);
        return res;
    }
} // namespace std

class lpp::Polynomial {
public:
    std::vector<Variable> terms;

    Polynomial(const std::vector<Variable>& terms) : terms(terms) {}

    Polynomial& operator+=(const Variable& value) {
        for (Variable& variable : terms) {
            if (variable.name == value.name && variable.exponent == value.exponent) {
                variable.coefficient += value.coefficient;
                return *this;
            }
        }
        terms.push_back(value);
        return *this;
    }

    Polynomial operator+(const Variable& value) const {
        Polynomial polynomial = *this;
        polynomial += value;
        return polynomial;
    }

    Polynomial& operator*=(const double value) {
        for (Variable& variable : terms) {
            variable *= value;
        }
        return *this;
    }

    Polynomial operator*(const double value) const {
        Polynomial polynomial = *this;
        polynomial *= value;
        return polynomial;
    }

    Polynomial partial_substitution();

    double complete_substitution();

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

class lpp::Inequation {
public:
    enum class Operator { LT, LE, GT, GE, EQ } opr;
    double rhs;
    Polynomial lhs;

    Inequation(const Polynomial& polynomial, const Operator opr, const double rhs) : opr(opr), rhs(rhs), lhs(polynomial) {}

    Inequation& operator*=(const double value) {
        lhs *= value;
        rhs *= value;

        if (value < 0) {
            switch (opr) {
            case Operator::LT:
                opr = Operator::GT;
                break;

            case Operator::LE:
                opr = Operator::GT;
                break;

            case Operator::GT:
                opr = Operator::LT;
                break;

            case Operator::GE:
                opr = Operator::LE;
                break;

            default:
                break;
            }
        }
        return *this;
    }

    Inequation operator*(const double value) const {
        Inequation inequation = *this;
        inequation *= value;
        return inequation;
    }

    friend std::ostream& operator<<(std::ostream& out, const Inequation& inequation) {
        out << inequation.lhs << ' ';

        switch (inequation.opr) {
        case Operator::LT:
            out << '<';
            break;

        case Operator::LE:
            out << "<=";
            break;

        case Operator::GT:
            out << '>';
            break;

        case Operator::GE:
            out << ">=";
            break;

        case Operator::EQ:
            out << '=';
        }
        return out << ' ' << inequation.rhs;
    }
};

class lpp::Equation : public Inequation {
public:
    Equation(const Polynomial& polynomial, const double rhs) : Inequation(polynomial, Operator::EQ, rhs) {}
};

class lpp::LPP {
public:
    OptimizeType type;
    Polynomial objective;
    std::vector<Inequation> constraints;

    LPP(const OptimizeType type, const Polynomial& objective, std::initializer_list<Inequation>&& constraints) :
        type(type), objective(objective), constraints(constraints) {}

    friend std::ostream& operator<<(std::ostream& out, const LPP& lpp) {
        out << (lpp.type == OptimizeType::MAXIMIZE ? "Maximize " : "Minimize ") << lpp.objective << std::endl;
        out << "subject to  " << lpp.constraints.front() << std::endl;

        for (const Inequation& constraint : lpp.constraints | std::views::drop(1)) {
            out << "\t    " << constraint << std::endl;
        }
        return out;
    }

    LPP& standardize() {
        if (type == OptimizeType::MINIMIZE) {
            objective *= -1;
        }
        int i = 1;

        for (Inequation& constraint : constraints) {
            if (constraint.rhs < 0) {
                constraint *= -1;
            }
            if (constraint.opr != Inequation::Operator::EQ) {
                constraint = Equation(constraint.lhs + Variable(std::string("s") + std::to_string(i++)), constraint.rhs);
            }
        }
        return *this;
    }
};


inline lpp::Polynomial operator+(const lpp::Variable& lhs, const lpp::Variable& rhs) { return lpp::Polynomial(std::vector{lhs, rhs}); }
inline lpp::Polynomial operator-(const lpp::Variable& lhs, const lpp::Variable& rhs) { return lpp::Polynomial(std::vector{lhs, -rhs}); }

inline lpp::Inequation operator<(const lpp::Polynomial& polynomial, const double value) {
    return lpp::Inequation(polynomial, lpp::Inequation::Operator::LT, value);
}

inline lpp::Inequation operator<=(const lpp::Polynomial& polynomial, const double value) {
    return lpp::Inequation(polynomial, lpp::Inequation::Operator::LE, value);
}

inline lpp::Inequation operator>(const lpp::Polynomial& polynomial, const double value) {
    return lpp::Inequation(polynomial, lpp::Inequation::Operator::GT, value);
}

inline lpp::Inequation operator>=(const lpp::Polynomial& polynomial, const double value) {
    return lpp::Inequation(polynomial, lpp::Inequation::Operator::GE, value);
}

inline lpp::Equation operator==(const lpp::Polynomial& polynomial, const double value) { return lpp::Equation(polynomial, value); }

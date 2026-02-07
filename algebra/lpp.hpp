#pragma once

namespace lpp {
    enum class Optimization { MAXIMIZE, MINIMIZE };
    using BV = std::vector<Variable>;
    using C = std::map<Variable, Fraction>;
    using B = std::vector<Fraction>;
    using CoefficientMatrix = std::map<Variable, std::vector<Fraction>>;
    using MR = std::map<Variable, Fraction>;
} // namespace lpp

class lpp::LPP {
public:
    Optimization type;
    Polynomial objective;
    std::vector<Inequation> constraints;

    LPP(const Optimization type, const Polynomial& objective, std::initializer_list<Inequation>&& constraints) :
        type(type), objective(objective), constraints(constraints) {}

    friend std::ostream& operator<<(std::ostream& out, const LPP& lpp) {
        out << (lpp.type == Optimization::MAXIMIZE ? "Maximize" : "Minimize") << "   " << lpp.objective << std::endl;
        out << "subject to " << lpp.constraints.front() << std::endl;

        for (const Inequation& constraint : lpp.constraints | std::views::drop(1)) {
            out << "\t   " << constraint << std::endl;
        }
        return out;
    }

    std::tuple<BV, C, B, CoefficientMatrix, MR> prepare_computational_table() {
        BV bv = standardize();
        C c;
        B b;
        CoefficientMatrix coefficient_matrix;

        for (const Variable& variable : objective.expression) {
            c[variable.basis()] = variable.coefficient;
        }
        for (const Inequation& constraint : constraints) {
            for (const Variable& variable : constraint.lhs.expression) {
                c.emplace(variable.basis(), 0);
            }
            b.push_back(constraint.rhs);
        }
        for (const Variable& variable : c | std::views::keys) {
            coefficient_matrix.emplace(variable, std::vector<Fraction>());
        }
        for (const Inequation& constraint : constraints) {
            for (std::vector<Fraction>& fractions : coefficient_matrix | std::views::values) {
                fractions.push_back(0);
            }
            for (const Variable& variable : constraint.lhs.expression) {
                coefficient_matrix[variable.basis()].back() = variable.coefficient;
            }
        }
        return {bv, c, b, coefficient_matrix, {}};
    }

    std::vector<Variable> standardize() {
        std::vector<Variable> res;
        int i = 1;

        if (type == Optimization::MINIMIZE) {
            objective *= -1;
        }
        for (Inequation& constraint : constraints) {
            if (constraint.rhs < 0) {
                constraint *= -1;
            }
            if (constraint.opr != Inequation::Operator::EQ) {
                res.push_back(Variable(std::string("s") + std::to_string(i++)));
                constraint = Equation(
                    constraint.lhs +
                        (constraint.opr == Inequation::Operator::LT || constraint.opr == Inequation::Operator::LE ? res.back() : -res.back()),
                    constraint.rhs);
            }
        }
        std::ranges::sort(res);
        return res;
    }
};

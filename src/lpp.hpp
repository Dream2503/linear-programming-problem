#pragma once

class lpp::LPP {
public:
    LPP() = default;

    LPP standardize(const bool dual = false) const {
        LPP lpp = *this;
        int i = 1;

        if (lpp.type == Optimization::MINIMIZE) {
            lpp.objective *= -1;
            lpp.type = Optimization::MAXIMIZE;
        }
        for (Inequation& constraint : lpp.constraints) {
            if (dual) {
                if (constraint.opr == RelationalOperator::GE) {
                    constraint = constraint.invert();
                }
            } else {
                if (static_cast<Fraction>(constraint.rhs) < 0) {
                    constraint = constraint.invert();
                }
            }
            if (constraint.opr != RelationalOperator::EQ) {
                Variable variable(std::string("s") + std::to_string(i++));
                constraint = Equation(constraint.lhs + (constraint.opr == RelationalOperator::LE ? variable : -variable), constraint.rhs);
            }
        }
        return lpp;
    }

    LPP canonicalize() const {
        LPP lpp = *this;

        for (const Inequation& constraint : constraints) {
            if (constraint.opr == RelationalOperator::EQ) {
                lpp.constraints.push_back(Inequation(constraint.lhs, RelationalOperator::LE, constraint.rhs));
                lpp.constraints.push_back(Inequation(constraint.lhs, RelationalOperator::GE, constraint.rhs));
            }
        }
        std::erase_if(lpp.constraints, [](const Inequation& inequation) -> bool { return inequation.opr == RelationalOperator::EQ; });

        for (Inequation& constraint : lpp.constraints) {
            if (type == Optimization::MAXIMIZE && constraint.opr == RelationalOperator::GE ||
                type == Optimization::MINIMIZE && constraint.opr == RelationalOperator::LE) {
                constraint = constraint.invert();
            }
        }
        return lpp;
    }

    inline static const Variable B{"@"}, M{"M"}, Z{"Z"}; // '@' for ascii arrangement
    Optimization type;
    Polynomial objective;
    std::vector<Inequation> constraints, restrictions;

    LPP(const Optimization type, const Polynomial& objective, const std::vector<Inequation>& constraints,
        const std::vector<Inequation>& restrictions) :
        type(type), objective(objective), constraints(constraints), restrictions(restrictions) {} // need to add variables integrity check

    static Inequation unrestrict(const Variable& variable) { return Inequation(variable, {}, inf); }

    ComputationalTable optimize(const std::string& = "simplex", bool = false) const;

    LPP dual(const std::string& = "w") const;

    friend std::ostream& operator<<(std::ostream& out, const LPP& lpp) {
        const int size = lpp.restrictions.size();
        out << (lpp.type == Optimization::MAXIMIZE ? "Maximize" : "Minimize") << "   " << lpp.objective << std::endl;
        out << "subject to " << lpp.constraints.front() << std::endl;

        for (const Inequation& constraint : lpp.constraints | std::views::drop(1)) {
            out << "\t   " << constraint << std::endl;
        }
        out << "\t   ";

        for (int i = 0; i < size; i++) {
            if (lpp.restrictions[i].lhs.is_fraction() && static_cast<Fraction>(lpp.restrictions[i].lhs) == inf ||
                lpp.restrictions[i].rhs.is_fraction() && static_cast<Fraction>(lpp.restrictions[i].rhs) == inf) {
                out << lpp.restrictions[i].lhs << " is unrestricted";
            } else {
                out << lpp.restrictions[i];
            }
            if (i < size - 1) {
                out << ", ";
            }
        }
        return out << std::endl;
    }
};

inline std::vector<std::map<algebra::Variable, algebra::Fraction>> lpp::basic_feasible_solutions(const std::vector<Equation>& equations) {
    std::map<Variable, std::vector<Fraction>> variables;

    for (const Equation& equation : equations) {
        for (const Variable& variable : equation.lhs.expression) {
            variables[variable.basis()].push_back(variable.coefficient);
        }
        variables[LPP::B].push_back(static_cast<Fraction>(equation.rhs));
    }
    const int col = std::ranges::max(variables | std::views::drop(1) | std::views::values, {}, &std::vector<Fraction>::size).size();
    const int row = variables.size() - 1, n = std::max(row, col), k = std::min(row, col);
    const std::vector<std::vector<int>> combinations = detail::generate_combination(n, k);
    std::vector<std::map<Variable, Fraction>> result;

    for (const std::vector<int>& combination : combinations) {
        Matrix<Fraction> B(k, k), C(k, 1);
        std::vector<Variable> X;
        X.reserve(k);

        for (int i = 0; i < k; i++) {
            const auto itr = std::next(variables.begin(), combination[i] + 1); // B

            for (int j = 0; j < k; j++) {
                B[j, i] = itr->second[j];
            }
            X.push_back(itr->first);
            C[i, 0] = variables[LPP::B][i];
        }
        Matrix<Fraction> res = B.inverse() * C;
        std::map<Variable, Fraction> element;

        for (int i = 0; i < k; i++) {
            element[X[i]] = res[i, 0];
        }
        result.push_back(element);
    }
    return result;
}

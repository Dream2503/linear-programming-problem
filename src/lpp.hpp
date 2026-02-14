#pragma once

namespace lpp {
    enum class Optimization { MAXIMIZE, MINIMIZE };
    using BV = std::vector<Variable>;
    using C = std::map<Variable, Variable>;
    using CoefficientMatrix = std::map<Variable, std::vector<Fraction>>;
    using Result = std::map<Variable, Fraction>;
} // namespace lpp

class lpp::LPP {
public:
    inline static const Variable B{"@"}, M{"M"}, Z{"Z"};
    Optimization type;
    Polynomial objective;
    std::vector<Inequation> constraints;

    LPP(const Optimization type, const Polynomial& objective, std::initializer_list<Inequation>&& constraints) :
        type(type), objective(objective), constraints(constraints) {} // need to add variables integrity check

    std::variant<Result, std::string> optimize() const { return compute(standardize().prepare_computational_table()); }

    friend std::ostream& operator<<(std::ostream& out, const LPP& lpp) {
        out << (lpp.type == Optimization::MAXIMIZE ? "Maximize" : "Minimize") << "   " << lpp.objective << std::endl;
        out << "subject to " << lpp.constraints.front() << std::endl;

        for (const Inequation& constraint : lpp.constraints | std::views::drop(1)) {
            out << "\t   " << constraint << std::endl;
        }
        return out;
    }

private:
    static Fraction extract_M_coefficient(const Polynomial& polynomial) {
        const auto itr = std::ranges::find(polynomial.expression, M, &Variable::basis);

        if (itr != polynomial.expression.end()) {
            return itr->coefficient;
        }
        return static_cast<Fraction>(polynomial);
    }

    LPP standardize() const {
        LPP lpp = *this;
        int i = 1;

        if (lpp.type == Optimization::MINIMIZE) {
            lpp.objective *= -1;
        }
        for (Inequation& constraint : lpp.constraints) {
            if (constraint.rhs < 0) {
                constraint *= -1;
            }
            if (constraint.opr != Inequation::Operator::EQ) {
                Variable variable(std::string("s") + std::to_string(i++));
                constraint = Equation(constraint.lhs + (constraint.opr == Inequation::Operator::LE ? variable : -variable), constraint.rhs);
            }
        }
        return lpp;
    }

    LPP canonicalize() const {
        LPP lpp = *this;

        for (const Inequation& constraint : constraints) {
            if (constraint.opr == Inequation::Operator::EQ) {
                lpp.constraints.push_back(Inequation(constraint.lhs, Inequation::Operator::LE, constraint.rhs));
                lpp.constraints.push_back(Inequation(constraint.lhs, Inequation::Operator::GE, constraint.rhs));
            }
        }
        std::erase_if(lpp.constraints, [](const Inequation& inequation) -> bool { return inequation.opr == Inequation::Operator::EQ; });

        for (Inequation& constraint : lpp.constraints) {
            if (type == Optimization::MAXIMIZE && constraint.opr == Inequation::Operator::GE ||
                type == Optimization::MINIMIZE && constraint.opr == Inequation::Operator::LE) {
                constraint *= -1;
            }
        }
        return lpp;
    }

    std::tuple<BV, C, CoefficientMatrix> prepare_computational_table() const {
        const int size = constraints.size();
        int j = 1;
        BV bv;
        C c;
        CoefficientMatrix coefficient_matrix;
        Matrix<Fraction> unit_matrix = Matrix<Fraction>::make_identity(size);

        for (const Variable& variable : objective.expression) {
            c.emplace(variable.basis(), variable.coefficient);
        }
        for (const Inequation& constraint : constraints) {
            for (const Variable& variable : constraint.lhs.expression) {
                c.emplace(variable.basis(), 0);
            }
            coefficient_matrix[B].push_back(constraint.rhs);
        }
        for (const Variable& variable : c | std::views::keys) {
            coefficient_matrix.emplace(variable, std::vector<Fraction>());
        }
        for (const Inequation& constraint : constraints) {
            for (std::vector<Fraction>& fractions : coefficient_matrix | std::views::drop(1) | std::views::values) { // B
                fractions.push_back(0);
            }
            for (const Variable& variable : constraint.lhs.expression) {
                coefficient_matrix[variable.basis()].back() = variable.coefficient;
            }
        }
        for (int i = 0; i < size; i++) {
            const auto itr = std::ranges::find_if(
                coefficient_matrix | std::views::drop(1),
                [&unit_matrix, i](const std::vector<Fraction>& fractions) -> bool { return fractions == unit_matrix[i]; },
                &std::map<Variable, std::vector<Fraction>>::value_type::second);

            if (itr != coefficient_matrix.end()) {
                bv.push_back(itr->first);
            } else {
                const Variable var("A" + std::to_string(j++));
                coefficient_matrix[var] = unit_matrix[i];
                c.emplace(var, -M);
                bv.push_back(var);
            }
        }
        return {bv, c, coefficient_matrix};
    }

    std::variant<Result, std::string> compute(const std::tuple<BV, C, CoefficientMatrix>& table) const {
        const int size = constraints.size();
        auto [bv, c, coefficient_matrix] = table;
        Result res;
        Matrix<Fraction> unit_matrix = Matrix<Fraction>::make_identity(size);

        while (true) {
            std::vector<Fraction> mr;
            std::vector<Polynomial> zj_cj;

            for (const Variable& variable : c | std::views::keys) {
                Polynomial polynomial;

                for (int i = 0; i < size; i++) {
                    polynomial += c[bv[i]] * coefficient_matrix[variable][i];
                }
                zj_cj.push_back(polynomial - c[variable]);
            }
            if (std::ranges::all_of(zj_cj, [](const Polynomial& polynomial) -> bool { return extract_M_coefficient(polynomial) >= 0; })) {
                if (std::ranges::contains(bv, 'A', [](const Variable& variable) -> char { return variable.variables[0].name[0]; })) {
                    return "No feasible solution";
                }
                for (const Variable& variable : objective.expression | std::views::transform(&Variable::basis)) {
                    const int idx = std::ranges::find(bv, variable) - bv.begin();
                    res[variable] = idx < size ? coefficient_matrix[B][idx] : 0;
                    res[Z] += static_cast<Fraction>(c[variable]) * res[variable];
                }
                res[Z] *= type == Optimization::MINIMIZE ? -1 : 1;
                return res;
            }
            const auto ev =
                std::next(coefficient_matrix.begin(), std::ranges::min_element(zj_cj, {}, extract_M_coefficient) - zj_cj.begin() + 1); // B

            for (int i = 0; i < size; i++) {
                mr.push_back(ev->second[i] < 0 ? INT32_MAX : coefficient_matrix[B][i] / ev->second[i]);
            }
            int lv = std::ranges::min_element(mr) - mr.begin();
            bool unbounded = true;

            if (!std::ranges::contains(bv, 'A', [](const Variable& variable) -> char { return variable.variables[0].name[0]; })) {
                for (int k = 0; k < size && mr[lv] != INT32_MAX; k++) {
                    std::vector<int> candidates;

                    for (int i = 0; i < size; i++) {
                        if (mr[i] == mr[lv]) {
                            candidates.push_back(i);
                        }
                    }
                    if (candidates.size() > 1) {
                        for (const int candidate : candidates) {
                            mr[candidate] = coefficient_matrix[bv[k]][candidate] / coefficient_matrix[ev->first][k];
                        }
                        lv = std::ranges::min_element(mr) - mr.begin();
                    } else {
                        unbounded = false;
                        break;
                    }
                }
            } else {
                unbounded = false;
            }
            if (unbounded || mr[lv] == INT32_MAX) {
                return "Unbounded Solution";
            }
            if (bv[lv].variables[0].name[0] == 'A') {
                coefficient_matrix.erase(bv[lv]);
                c.erase(bv[lv]);
            }
            CoefficientMatrix new_coefficient_matrix = coefficient_matrix;
            const std::pair pivot = {ev->first, lv};
            bv.erase(bv.begin() + lv);
            bv.insert(bv.begin() + lv, ev->first);

            for (int i = 0; i < size; i++) {
                new_coefficient_matrix[bv[i]] = unit_matrix[i];
            }
            for (const Variable& variable : coefficient_matrix | std::views::keys |
                     std::views::filter([&bv](const Variable& element) -> bool { return !std::ranges::contains(bv, element); })) {
                for (int i = 0; i < size; i++) {
                    if (i == pivot.second) {
                        new_coefficient_matrix[variable][i] /= coefficient_matrix[pivot.first][pivot.second];
                    } else {
                        new_coefficient_matrix[variable][i] = (coefficient_matrix[variable][i] * coefficient_matrix[pivot.first][pivot.second] -
                                                               coefficient_matrix[variable][pivot.second] * coefficient_matrix[pivot.first][i]) /
                            coefficient_matrix[pivot.first][pivot.second];
                    }
                }
            }
            coefficient_matrix = std::move(new_coefficient_matrix);
        }
    }
};

inline std::vector<lpp::Result> lpp::basic_feasible_solutions(const std::vector<Equation>& equations) {
    std::map<Variable, std::vector<Fraction>> variables;

    for (const Equation& equation : equations) {
        for (const Variable& variable : equation.lhs.expression) {
            variables[variable.basis()].push_back(variable.coefficient);
        }
        variables[LPP::B].push_back(equation.rhs);
    }
    const int col = std::ranges::max(variables | std::views::drop(1) | std::views::values, {}, &std::vector<Fraction>::size).size();
    const int row = variables.size() - 1, n = std::max(row, col), k = std::min(row, col);
    const std::vector<std::vector<int>> combinations = detail::generate_combination(n, k);
    std::vector<Result> result;

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
        Result element;

        for (int i = 0; i < k; i++) {
            element[X[i]] = res[i, 0];
        }
        result.push_back(element);
    }
    return result;
}

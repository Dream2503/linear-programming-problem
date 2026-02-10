#pragma once

namespace lpp {
    enum class Optimization { MAXIMIZE, MINIMIZE };
    using BV = std::vector<Variable>;
    using C = std::map<Variable, Variable>;
    using CoefficientMatrix = std::map<Variable, std::vector<Fraction>>;
} // namespace lpp

class lpp::LPP {
public:
    static constexpr auto B = Variable(""), M = Variable("M");
    Optimization type;
    Polynomial objective;
    std::vector<Inequation> constraints;

    LPP(const Optimization type, const Polynomial& objective, std::initializer_list<Inequation>&& constraints) :
        type(type), objective(objective), constraints(constraints) {}

    std::variant<std::map<Variable, Fraction>, std::string> optimize() const { return compute(standardize().prepare_computational_table()); }

    friend std::ostream& operator<<(std::ostream& out, const LPP& lpp) {
        out << (lpp.type == Optimization::MAXIMIZE ? "Maximize" : "Minimize") << "   " << lpp.objective << std::endl;
        out << "subject to " << lpp.constraints.front() << std::endl;

        for (const Inequation& constraint : lpp.constraints | std::views::drop(1)) {
            out << "\t   " << constraint << std::endl;
        }
        return out;
    }

private:
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
                constraint =
                    Equation(constraint.lhs +
                                 (constraint.opr == Inequation::Operator::LT || constraint.opr == Inequation::Operator::LE ? variable : -variable),
                             constraint.rhs);
            }
        }
        return lpp;
    }

    std::tuple<BV, C, CoefficientMatrix> prepare_computational_table() const {
        const int size = constraints.size();
        int j = 1;
        std::vector unit_matrix(size, std::vector<Fraction>(size));
        BV bv;
        C c;
        CoefficientMatrix coefficient_matrix;

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
            unit_matrix[i][i] = 1;
        }
        for (int i = 0; i < size; i++) {
            const auto itr = std::ranges::find_if(
                coefficient_matrix, [&unit_matrix, i](const std::vector<Fraction>& fractions) -> bool { return fractions == unit_matrix[i]; },
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

    std::variant<std::map<Variable, Fraction>, std::string> compute(const std::tuple<BV, C, CoefficientMatrix>& table) const {
        const int size = constraints.size();
        auto [bv, c, coefficient_matrix] = table;
        std::map<Variable, Fraction> res;
        std::vector unit_matrix(size, std::vector<Fraction>(size));

        const auto extract_M_coefficient = [](const Polynomial& polynomial) -> Fraction {
            const auto itr = std::ranges::find(polynomial.expression, M, [](const Variable& variable) -> Variable { return variable.basis(); });

            if (itr != polynomial.expression.end()) {
                return itr->coefficient;
            }
            return static_cast<Fraction>(polynomial);
        };

        for (int i = 0; i < size; i++) {
            unit_matrix[i][i] = 1;
        }
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
            if (std::ranges::all_of(
                    zj_cj, [extract_M_coefficient](const Polynomial& polynomial) -> bool { return extract_M_coefficient(polynomial) >= 0; })) {
                if (std::ranges::contains(bv, 'A', [](const Variable& variable) -> char { return variable.name[0]; })) {
                    return "No feasible solution";
                }
                for (const Variable& variable :
                     objective.expression | std::views::transform([](const Variable& var) -> Variable { return var.basis(); })) {
                    const int idx = std::ranges::find(bv, variable) - bv.begin();
                    res[variable] = idx < size ? coefficient_matrix[B][idx] : 0;
                }
                return res;
            }
            const auto ev = std::next(coefficient_matrix.begin(),
                                      std::ranges::min_element(zj_cj, {},
                                                               extract_M_coefficient) -
                                          zj_cj.begin() + 1); // B

            for (int i = 0; i < size; i++) {
                mr.push_back(ev->second[i] < 0 ? INT32_MAX : coefficient_matrix[B][i] / ev->second[i]);
            }
            const int lv = std::ranges::min_element(mr) - mr.begin();

            if (mr[lv] == INT32_MAX) {
                return "Unbounded Solution";
            }
            const std::pair pivot = {ev->first, lv};

            if (bv[lv].name[0] == 'A') {
                coefficient_matrix.erase(bv[lv]);
                c.erase(bv[lv]);
            }
            bv.erase(bv.begin() + lv);
            bv.insert(bv.begin() + lv, ev->first);
            CoefficientMatrix new_coefficient_matrix = coefficient_matrix;

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
            coefficient_matrix.swap(new_coefficient_matrix);
        }
    }
};

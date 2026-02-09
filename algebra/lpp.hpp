#pragma once

namespace lpp {
    enum class Optimization { MAXIMIZE, MINIMIZE };
    using BV = std::vector<Variable>;
    using C = std::map<Variable, Fraction>;
    using CoefficientMatrix = std::map<Variable, std::vector<Fraction>>;
} // namespace lpp

class lpp::LPP {
public:
    static constexpr auto B = Variable("B");
    Optimization type;
    Polynomial objective;
    std::vector<Inequation> constraints;

    LPP(const Optimization type, const Polynomial& objective, std::initializer_list<Inequation>&& constraints) :
        type(type), objective(objective), constraints(constraints) {}

    std::map<Variable, Fraction> optimize() { return compute(standardize().prepare_computational_table()); }

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
        BV bv;
        C c;
        CoefficientMatrix coefficient_matrix;
        const int size = constraints.size();
        std::vector unit_matrix(size, std::vector<Fraction>(size));

        for (const Variable& variable : objective.expression) {
            c[variable.basis()] = variable.coefficient;
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
            bv.push_back(std::ranges::find_if(
                             coefficient_matrix,
                             [&unit_matrix, i](const std::vector<Fraction>& fractions) -> bool { return fractions == unit_matrix[i]; },
                             [](const std::pair<Variable, std::vector<Fraction>>& element) -> std::vector<Fraction> { return element.second; })
                             ->first);
        }
        return {bv, c, coefficient_matrix};
    }

    std::map<Variable, Fraction> compute(const std::tuple<BV, C, CoefficientMatrix>& table) const {
        const int size = constraints.size();
        auto [bv, c, coefficient_matrix] = table;
        std::map<Variable, Fraction> res;
        std::vector unit_matrix(size, std::vector<Fraction>(size));

        for (int i = 0; i < size; i++) {
            unit_matrix[i][i] = 1;
        }
        while (true) {
            std::vector<Fraction> zj_cj, mr;

            for (const Variable& variable : c | std::views::keys) {
                Fraction fraction;

                for (int i = 0; i < size; i++) {
                    fraction += c[bv[i]] * coefficient_matrix[variable][i];
                }
                zj_cj.push_back(fraction - c[variable]);
            }
            if (std::ranges::all_of(zj_cj, [](const Fraction& fraction) -> bool { return fraction >= 0; })) {
                break;
            }
            const auto entering_variable = std::next(coefficient_matrix.begin(), std::ranges::min_element(zj_cj) - zj_cj.begin() + 1); // B

            for (int i = 0; i < size; i++) {
                mr.push_back(entering_variable->second[i] < 0 ? INT32_MAX : coefficient_matrix[B][i] / entering_variable->second[i]);
            }
            const int leaving_variable = std::ranges::min_element(mr) - mr.begin();
            const std::pair pivot = {entering_variable->first, leaving_variable};
            bv.erase(bv.begin() + leaving_variable);
            bv.insert(bv.begin() + leaving_variable, entering_variable->first);
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
        for (const Variable& variable : objective.expression | std::views::transform([](const Variable& var) -> Variable { return var.basis(); })) {
            const int idx = std::ranges::find(bv, variable) - bv.begin();
            res[variable] = idx < size ? coefficient_matrix[B][idx] : 0;
        }
        return res;
    }
};

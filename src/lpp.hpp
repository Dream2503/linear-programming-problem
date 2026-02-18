#pragma once

class lpp::LPP {
    LPP() = default;

    LPP standardize() const {
        LPP lpp = *this;
        int i = 1;

        if (lpp.type == Optimization::MINIMIZE) {
            lpp.objective *= -1;
        }
        for (Inequation& constraint : lpp.constraints) {
            if (static_cast<Fraction>(constraint.rhs) < 0) {
                constraint *= -1;
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
                constraint *= -1;
            }
        }
        return lpp;
    }

public:
    inline static const Variable B{"@"}, M{"M"}, Z{"Z"};
    Optimization type;
    Polynomial objective;
    std::vector<Inequation> constraints, restrictions;

    LPP(const Optimization type, const Polynomial& objective, const std::vector<Inequation>& constraints,
        const std::vector<Inequation>& restrictions) :
        type(type), objective(objective), constraints(constraints), restrictions(restrictions) {} // need to add variables integrity check

    static Inequation unrestrict(const Variable& variable) { return Inequation(variable, {}, inf); }

    std::variant<std::map<Variable, Fraction>, std::string> optimize(bool debug = false) const;

    LPP dual(const std::string& basis = "w") const;

    friend std::ostream& operator<<(std::ostream& out, const LPP& lpp) {
        const int size = lpp.restrictions.size();
        out << (lpp.type == Optimization::MAXIMIZE ? "Maximize" : "Minimize") << "   " << lpp.objective << std::endl;
        out << "subject to " << lpp.constraints.front() << std::endl;

        for (const Inequation& constraint : lpp.constraints | std::views::drop(1)) {
            out << "\t   " << constraint << std::endl;
        }
        out << "\t   ";

        for (int i = 0; i < size; i++) {
            if (static_cast<Fraction>(lpp.restrictions[i].rhs) == inf) {
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

class lpp::ComputationalTable {
    static Fraction extract_M_coefficient(const Polynomial& polynomial) {
        const auto itr = std::ranges::find(polynomial.expression, LPP::M, &Variable::basis);

        if (itr != polynomial.expression.end()) {
            return itr->coefficient;
        }
        return static_cast<Fraction>(polynomial);
    }

public:
    Solution solution;
    std::vector<Variable> basis_vector;
    std::map<Variable, Variable> cost;
    std::map<Variable, std::vector<Fraction>> coefficient_matrix;
    std::vector<Polynomial> zj_cj;
    std::vector<Fraction> mr;

    ComputationalTable() = default;

    explicit ComputationalTable(const LPP& lpp) : solution(Solution::UNOPTIMIZED) {
        const int size = lpp.constraints.size();
        int j = 1;
        Matrix<Fraction> unit_matrix = Matrix<Fraction>::make_identity(size);

        for (const Variable& variable : lpp.objective.expression) {
            cost.emplace(variable.basis(), variable.coefficient);
        }
        for (const Inequation& constraint : lpp.constraints) {
            for (const Variable& variable : constraint.lhs.expression) {
                cost.emplace(variable.basis(), 0);
            }
            coefficient_matrix[LPP::B].push_back(static_cast<Fraction>(constraint.rhs));
        }
        for (const Variable& variable : cost | std::views::keys) {
            coefficient_matrix.emplace(variable, std::vector<Fraction>());
        }
        for (const Inequation& constraint : lpp.constraints) {
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
                basis_vector.push_back(itr->first);
            } else {
                const Variable var("A" + std::to_string(j++));
                coefficient_matrix[var] = unit_matrix[i];
                cost.emplace(var, -LPP::M);
                basis_vector.push_back(var);
            }
        }
    }

    Solution optimize(const bool debug = false) {
        if (solution != Solution::UNOPTIMIZED) {
            return solution;
        }
        const int size = coefficient_matrix[LPP::B].size();
        Matrix<Fraction> unit_matrix = Matrix<Fraction>::make_identity(size);

        while (true) {
            zj_cj.clear();
            mr.clear();

            for (const Variable& variable : cost | std::views::keys) {
                Polynomial polynomial;

                for (int i = 0; i < size; i++) {
                    polynomial += cost[basis_vector[i]] * coefficient_matrix[variable][i];
                }
                zj_cj.push_back(polynomial - cost[variable]);
            }
            if (std::ranges::all_of(zj_cj, [](const Polynomial& polynomial) -> bool { return extract_M_coefficient(polynomial) >= 0; })) {
                return solution =
                           std::ranges::contains(basis_vector, 'A', [](const Variable& variable) -> char { return variable.variables[0].name[0]; })
                    ? Solution::INFEASIBLE
                    : Solution::OPTIMIZED;
            }
            const auto ev =
                std::next(coefficient_matrix.begin(), std::ranges::min_element(zj_cj, {}, extract_M_coefficient) - zj_cj.begin() + 1); // B

            for (int i = 0; i < size; i++) {
                mr.push_back(ev->second[i] < 0 ? inf : coefficient_matrix[LPP::B][i] / ev->second[i]);
            }
            if (debug) {
                std::cout << *this;
            }
            int lv = std::ranges::min_element(mr) - mr.begin();
            bool unbounded = true;

            if (!std::ranges::contains(basis_vector, 'A', [](const Variable& variable) -> char { return variable.variables[0].name[0]; })) {
                for (int k = 0; k < size && mr[lv] != inf; k++) {
                    std::vector<int> candidates;

                    for (int i = 0; i < size; i++) {
                        if (mr[i] == mr[lv]) {
                            candidates.push_back(i);
                        }
                    }
                    if (candidates.size() > 1) {
                        for (const int candidate : candidates) {
                            mr[candidate] = coefficient_matrix[basis_vector[k]][candidate] / coefficient_matrix[ev->first][k];
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
            if (unbounded || mr[lv] == inf) {
                return solution = Solution::UNBOUNDED;
            }
            if (basis_vector[lv].variables[0].name[0] == 'A') {
                coefficient_matrix.erase(basis_vector[lv]);
                cost.erase(basis_vector[lv]);
            }
            const std::pair pivot = {ev->first, lv};
            std::map<Variable, std::vector<Fraction>> new_coefficient_matrix = coefficient_matrix;
            basis_vector.erase(basis_vector.begin() + lv);
            basis_vector.insert(basis_vector.begin() + lv, ev->first);

            for (int i = 0; i < size; i++) {
                new_coefficient_matrix[basis_vector[i]] = unit_matrix[i];
            }
            for (const Variable& variable : coefficient_matrix | std::views::keys |
                     std::views::filter([this](const Variable& element) -> bool { return !std::ranges::contains(basis_vector, element); })) {
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

    friend std::ostream& operator<<(std::ostream& out, const ComputationalTable& computational_table) {
        static constexpr int TAB_SIZE = 11;
        const int size = computational_table.basis_vector.size(), columns = 3 + computational_table.coefficient_matrix.size();
        auto print_partition = [columns, &out]() -> void {
            out << std::right << std::setfill('-') << '+';

            for (int i = 0; i < columns; i++) {
                out << std::setw(TAB_SIZE) << "" << '+';
            }
            out << std::endl << std::left << std::setfill(' ');
        };
        out << std::left << ' ';

        for (int i = 0; i < 3; i++) {
            out << std::setw(TAB_SIZE) << "" << ' ';
        }
        for (const Variable& variable : computational_table.coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            out << std::setw(TAB_SIZE) << computational_table.cost.at(variable) << ' ';
        }
        out << std::endl;
        print_partition();
        out << '|';

        for (const char* string : {"BV", "C", "B"}) {
            out << std::setw(TAB_SIZE) << string << '|';
        }
        for (const Variable& variable : computational_table.coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            out << std::setw(TAB_SIZE) << variable.variables[0].name << '|';
        }
        out << std::setw(TAB_SIZE) << "MR" << '|' << std::endl;
        print_partition();

        for (int i = 0; i < size; i++) {
            out << '|' << std::setw(TAB_SIZE) << computational_table.basis_vector[i] << '|' << std::setw(TAB_SIZE)
                << computational_table.cost.at(computational_table.basis_vector[i]) << '|';

            for (const std::vector<Fraction>& fractions : computational_table.coefficient_matrix | std::views::values) {
                out << std::setw(TAB_SIZE) << fractions[i] << '|';
            }
            if (computational_table.solution == Solution::UNOPTIMIZED) {
                out << std::setw(TAB_SIZE) << (computational_table.mr[i] == inf ? "-ve" : std::to_string(computational_table.mr[i])) << '|';
            } else {
                out << std::setw(TAB_SIZE) << "" << '|';
            }
            out << std::endl;
        }
        print_partition();
        out << '|';

        for (int i = 0; i < 2; i++) {
            out << std::setw(TAB_SIZE) << "" << '|';
        }
        out << std::setw(TAB_SIZE) << "Zj-Cj" << '|';

        for (const Polynomial& polynomial : computational_table.zj_cj) {
            out << std::setw(TAB_SIZE) << polynomial << '|';
        }
        out << std::setw(TAB_SIZE) << "" << '|' << std::endl;
        print_partition();
        return out << std::endl;
    }
};

inline lpp::LPP lpp::LPP::dual(const std::string& basis) const {
    LPP canonical = *this;

    for (Inequation& constraint : canonical.constraints) {
        if (type == Optimization::MAXIMIZE && constraint.opr == RelationalOperator::GE ||
            type == Optimization::MINIMIZE && constraint.opr == RelationalOperator::LE) {
            constraint *= -1;
        }
    }
    const int objective_size = canonical.objective.expression.size(), constraints_size = canonical.constraints.size();
    LPP res;
    ComputationalTable computational_table(*this);
    res.type = canonical.type == Optimization::MAXIMIZE ? Optimization::MINIMIZE : Optimization::MAXIMIZE;
    res.constraints.resize(objective_size);
    res.restrictions.reserve(constraints_size);

    for (int i = 0; i < constraints_size; i++) {
        auto itr = std::next(computational_table.coefficient_matrix.begin()); // B
        Variable variable(basis + std::to_string(i + 1));
        res.objective += computational_table.coefficient_matrix[B][i] * variable;

        for (int j = 0; j < objective_size; j++) {
            res.constraints[j].lhs += itr->second[i] * variable;
            ++itr;
        }
        res.restrictions.push_back(canonical.constraints[i].opr == RelationalOperator::EQ ? unrestrict(variable) : variable >= 0);
    }
    auto itr = computational_table.cost.begin();

    for (int i = 0; i < objective_size; i++) {
        res.constraints[i].opr = static_cast<Fraction>(canonical.restrictions[i].rhs) == inf ? RelationalOperator::EQ
            : res.type == Optimization::MAXIMIZE                                             ? RelationalOperator::LE
                                                                                             : RelationalOperator::GE;
        res.constraints[i].rhs = itr->second;
        ++itr;
    }
    return res;
}

inline std::variant<std::map<algebra::Variable, algebra::Fraction>, std::string> lpp::LPP::optimize(const bool debug) const {
    const LPP lpp = standardize();

    if (debug) {
        cout << std::endl << "Standard Form: " << std::endl << lpp;
    }
    ComputationalTable computational_table(lpp);

    switch (computational_table.optimize(debug)) {
    case Solution::OPTIMIZED:
        {
            const int size = constraints.size();
            std::map<Variable, Fraction> res;

            for (const Variable& variable : objective.expression | std::views::transform(&Variable::basis)) {
                const int idx = std::ranges::find(computational_table.basis_vector, variable) - computational_table.basis_vector.begin();
                res[variable] = idx < size ? computational_table.coefficient_matrix[B][idx] : 0;
                res[Z] += static_cast<Fraction>(computational_table.cost[variable]) * res[variable];
            }
            res[Z] *= type == Optimization::MINIMIZE ? -1 : 1;
            if (debug) {
                std::cout << computational_table;
            }
            return res;
        }

    case Solution::INFEASIBLE:
        return "In Feasible Solution";

    case Solution::UNBOUNDED:
        return "Unbounded Solution";

    case Solution::UNOPTIMIZED:
        std::unreachable();
    }
    std::unreachable();
}

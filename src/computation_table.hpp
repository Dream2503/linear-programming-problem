#pragma once

class optimization::ComputationalTable {
    static algebra::Fraction extract_coefficient_M(const algebra::Polynomial& polynomial) {
        const auto itr = std::ranges::find(polynomial.expression, LPP::M, &algebra::Variable::basis);

        if (itr != polynomial.expression.end()) {
            return itr->coefficient;
        }
        return static_cast<algebra::Fraction>(polynomial);
    }

    void compute_zj_cj() {
        const int size = coefficient_matrix[LPP::B].size();
        zj_cj.clear();

        for (const algebra::Variable& variable : cost | std::views::keys) {
            algebra::Polynomial polynomial;

            for (int i = 0; i < size; i++) {
                polynomial += cost[basis_vector[i]] * coefficient_matrix[variable][i];
            }
            zj_cj.push_back(polynomial - cost[variable]);
        }
    }

public:
    LPP lpp;
    Solution solution;
    std::vector<algebra::Variable> basis_vector;
    std::map<algebra::Variable, algebra::Variable> cost;
    std::map<algebra::Variable, std::vector<algebra::Fraction>> coefficient_matrix;
    std::vector<algebra::Polynomial> zj_cj;
    std::vector<algebra::Fraction> mr;

    explicit ComputationalTable(const LPP& lpp) : lpp(lpp), solution(Solution::UNOPTIMIZED) {
        const int size = lpp.constraints.size();
        linalg::Matrix<algebra::Fraction> unit_matrix = linalg::Matrix<algebra::Fraction>::make_identity(size);

        for (const algebra::Variable& variable : lpp.objective.expression) {
            cost.emplace(variable.basis(), variable.coefficient);
        }
        for (const algebra::Inequation& constraint : lpp.constraints) {
            for (const algebra::Variable& variable : constraint.lhs.expression) {
                cost.emplace(variable.basis(), 0);
            }
            coefficient_matrix[LPP::B].push_back(static_cast<algebra::Fraction>(constraint.rhs));
        }
        for (const algebra::Variable& variable : cost | std::views::keys) {
            coefficient_matrix.emplace(variable, std::vector<algebra::Fraction>());
        }
        for (const algebra::Inequation& constraint : lpp.constraints) {
            for (std::vector<algebra::Fraction>& fractions : coefficient_matrix | std::views::drop(1) | std::views::values) { // B
                fractions.push_back(0);
            }
            for (const algebra::Variable& variable : constraint.lhs.expression) {
                coefficient_matrix[variable.basis()].back() = variable.coefficient;
            }
        }
        for (int i = 0, j = 1; i < size; i++) {
            const auto itr =
                std::ranges::find_if(coefficient_matrix | std::views::drop(1), // B
                                     [&unit_matrix, i](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& element) -> bool {
                                         return element.second == unit_matrix[i];
                                     });

            if (itr != coefficient_matrix.end()) {
                basis_vector.push_back(itr->first);
            } else {
                const algebra::Variable variable("A" + std::to_string(j++));
                coefficient_matrix[variable] = unit_matrix[i];
                cost.emplace(variable, -LPP::M);
                basis_vector.push_back(variable);
            }
        }
    }

    ComputationalTable(const std::map<algebra::Variable, algebra::Variable>& cost, const std::vector<algebra::Variable>& basis_vector,
                       const std::map<algebra::Variable, std::vector<algebra::Fraction>>& coefficient_matrix, const Solution solution,
                       const LPP& lpp = LPP()) :
        lpp(lpp), solution(solution), basis_vector(basis_vector), cost(cost), coefficient_matrix(coefficient_matrix) {
        compute_zj_cj();
    }

    std::variant<std::vector<std::map<algebra::Variable, algebra::Fraction>>, std::string> get_solutions(const std::string& method = "simplex") {
        auto add_solution = [this]() -> std::map<algebra::Variable, algebra::Fraction> {
            std::map<algebra::Variable, algebra::Fraction> res;
            const int size = basis_vector.size();

            for (const algebra::Variable& variable :
                 cost | std::views::keys | std::views::filter([](const algebra::Variable& var) -> bool { return var.variables[0].name[0] != 's'; })) {
                const int idx = std::ranges::find(basis_vector, variable) - basis_vector.begin();
                res[variable] = idx < size ? coefficient_matrix[LPP::B][idx] : 0;
                res[LPP::Z] += static_cast<algebra::Fraction>(cost[variable]) * res[variable];
            }
            res[LPP::Z] *= lpp.type == Optimization::MINIMIZE ? -1 : 1;

            if (GLOBAL_FORMATTING.verbose) {
                *GLOBAL_FORMATTING.out << *this;
            }
            return res;
        };
        std::vector<std::map<algebra::Variable, algebra::Fraction>> res;

        while (true) {
            solution = method == "simplex" ? optimize_simplex() : optimize_dual_simplex();

            switch (solution) {
            case Solution::OPTIMIZED:
                res.push_back(add_solution());
                return res;

            case Solution::INFEASIBLE:
                return "Infeasible Solution\n";

            case Solution::UNBOUNDED:
                return "Unbounded Solution\n";

            case Solution::ALTERNATE:
                {
                    std::map<algebra::Variable, algebra::Fraction> sol = add_solution();

                    if (!std::ranges::contains(res, sol)) {
                        res.push_back(sol);
                    } else {
                        return res;
                    }
                    break;
                }

            case Solution::UNOPTIMIZED:
                std::unreachable();
            }
        }
        return res;
        std::unreachable();
    }

    Solution optimize_simplex() {
        if (solution != Solution::UNOPTIMIZED && solution != Solution::ALTERNATE) {
            return solution;
        }
        const int size = coefficient_matrix[LPP::B].size();
        linalg::Matrix<algebra::Fraction> unit_matrix = linalg::Matrix<algebra::Fraction>::make_identity(size);

        while (true) {
            if (solution != Solution::ALTERNATE) {
                compute_zj_cj();

                if (std::ranges::all_of(zj_cj,
                                        [](const algebra::Polynomial& polynomial) -> bool { return extract_coefficient_M(polynomial) >= 0; })) {
                    solution = std::ranges::contains(basis_vector, 'A',
                                                     [](const algebra::Variable& variable) -> char { return variable.variables[0].name[0]; })
                        ? Solution::INFEASIBLE
                        : Solution::OPTIMIZED;
                    const int zj_cj_size = zj_cj.size();
                    auto itr = std::next(coefficient_matrix.begin()); // B

                    for (int i = 0; i < zj_cj_size; i++) {
                        if (zj_cj[i].is_fraction() && static_cast<algebra::Fraction>(zj_cj[i]) == 0 &&
                            !std::ranges::contains(basis_vector, itr->first)) {
                            solution = Solution::ALTERNATE;
                            break;
                        }
                        ++itr;
                    }
                    return solution;
                }
            }
            auto ev = std::next(coefficient_matrix.begin()); // B
            mr.clear();

            if (solution == Solution::ALTERNATE) {
                const int zj_cj_size = zj_cj.size();

                for (int i = 0; i < zj_cj_size; i++) {
                    if (static_cast<algebra::Fraction>(zj_cj[i]) == 0 && !std::ranges::contains(basis_vector, ev->first)) {
                        solution = Solution::UNOPTIMIZED;
                        break;
                    }
                    ++ev;
                }
            } else {
                ev = std::next(coefficient_matrix.begin(), std::ranges::min_element(zj_cj, {}, extract_coefficient_M) - zj_cj.begin() + 1); // B
            }
            for (int i = 0; i < size; i++) {
                mr.push_back(ev->second[i] <= 0 ? algebra::inf : coefficient_matrix[LPP::B][i] / ev->second[i]);
            }
            if (GLOBAL_FORMATTING.verbose) {
                *GLOBAL_FORMATTING.out << *this;
            }
            int lv = std::ranges::min_element(mr) - mr.begin();
            bool unbounded = true;

            if (!std::ranges::contains(basis_vector, 'A', [](const algebra::Variable& variable) -> char { return variable.variables[0].name[0]; })) {
                for (int k = 0; k < size && mr[lv] != algebra::inf; k++) {
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
            if (unbounded || mr[lv] == algebra::inf) {
                return solution = Solution::UNBOUNDED;
            }
            if (cost[basis_vector[lv]].variables == LPP::M.variables) {
                coefficient_matrix.erase(basis_vector[lv]);
                cost.erase(basis_vector[lv]);
                auto itr = std::ranges::find(lpp.objective.expression, basis_vector[lv], &algebra::Variable::basis);

                if (itr != lpp.objective.expression.end()) {
                    lpp.objective.expression.erase(itr);
                }
            }
            const std::pair pivot = {ev->first, lv};
            std::map<algebra::Variable, std::vector<algebra::Fraction>> new_coefficient_matrix = coefficient_matrix;
            basis_vector.erase(basis_vector.begin() + lv);
            basis_vector.insert(basis_vector.begin() + lv, ev->first);

            for (int i = 0; i < size; i++) {
                new_coefficient_matrix[basis_vector[i]] = unit_matrix[i];
            }
            for (const algebra::Variable& variable : coefficient_matrix | std::views::keys |
                     std::views::filter([this](const algebra::Variable& element) -> bool { return !std::ranges::contains(basis_vector, element); })) {
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

    Solution optimize_dual_simplex() {
        solution = Solution::UNBOUNDED;
        const int size = coefficient_matrix[LPP::B].size();
        linalg::Matrix<algebra::Fraction> unit_matrix = linalg::Matrix<algebra::Fraction>::make_identity(size);

        while (true) {
            compute_zj_cj();

            if (std::ranges::all_of(zj_cj,
                                    [](const algebra::Polynomial& polynomial) -> bool { return static_cast<algebra::Fraction>(polynomial) >= 0; }) &&
                std::ranges::all_of(coefficient_matrix[LPP::B], [](const algebra::Fraction& fraction) -> bool { return fraction >= 0; })) {
                return solution = Solution::OPTIMIZED;
            }
            if (GLOBAL_FORMATTING.verbose) {
                *GLOBAL_FORMATTING.out << *this;
            }
            const int lv = std::ranges::min_element(coefficient_matrix[LPP::B]) - coefficient_matrix[LPP::B].begin();
            const auto begin = std::next(coefficient_matrix.begin()); // B

            auto ev = std::ranges::max_element(coefficient_matrix | std::views::drop(1),
                                               [&](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& lhs,
                                                   const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& rhs) -> bool {
                                                   const algebra::Fraction &lhs_value = lhs.second[lv], rhs_value = rhs.second[lv];

                                                   if (lhs_value >= 0) {
                                                       return true;
                                                   }
                                                   if (rhs_value >= 0) {
                                                       return false;
                                                   }
                                                   const int lhs_idx = std::distance(begin, coefficient_matrix.find(lhs.first));
                                                   const int rhs_idx = std::distance(begin, coefficient_matrix.find(rhs.first));
                                                   return static_cast<algebra::Fraction>(zj_cj[lhs_idx]) / lhs_value <
                                                       static_cast<algebra::Fraction>(zj_cj[rhs_idx]) / rhs_value;
                                               });
            assert(ev != coefficient_matrix.end());

            const std::pair pivot = {ev->first, lv};
            std::map<algebra::Variable, std::vector<algebra::Fraction>> new_coefficient_matrix = coefficient_matrix;
            basis_vector.erase(basis_vector.begin() + lv);
            basis_vector.insert(basis_vector.begin() + lv, ev->first);

            for (int i = 0; i < size; i++) {
                new_coefficient_matrix[basis_vector[i]] = unit_matrix[i];
            }
            for (const algebra::Variable& variable : coefficient_matrix | std::views::keys |
                     std::views::filter([this](const algebra::Variable& element) -> bool { return !std::ranges::contains(basis_vector, element); })) {
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

    std::vector<algebra::Interval> cost_variation() {
        std::vector<algebra::Interval> res;

        if (solution == Solution::UNOPTIMIZED) {
            optimize_simplex();
        }
        for (const auto& [variable, value] : cost | std::views::filter([](const std::pair<algebra::Variable, algebra::Variable>& element) -> bool {
                                                 return element.first.variables[0].name[0] != 's';
                                             })) {
            const auto itr = std::ranges::find(basis_vector, variable);
            algebra::Variable var("C" + variable.variables[0].name);

            if (itr != basis_vector.end()) {
                const int idx = itr - basis_vector.begin();
                int i = 0;
                algebra::Fraction min = algebra::inf, max = -algebra::inf;

                for (const auto& [name, fractions] : coefficient_matrix | std::views::drop(1)) { // B
                    if (!std::ranges::contains(basis_vector, name)) {
                        if (fractions[idx] > 0) {
                            max = std::max(max, -static_cast<algebra::Fraction>(zj_cj[i] / fractions[idx]));
                        } else if (fractions[idx] < 0) {
                            min = std::min(min, -static_cast<algebra::Fraction>(zj_cj[i] / fractions[idx]));
                        }
                    }
                    i++;
                }
                res.push_back(max + value < var < min + value);
            } else {
                res.push_back(-algebra::inf < var <
                              zj_cj[std::distance(coefficient_matrix.begin(), coefficient_matrix.find(variable)) - 1] + value); // B
            }
        }
        return res;
    }

    std::vector<algebra::Interval> RHS_variation() {
        const int size = basis_vector.size();
        int i = 0;
        std::vector<algebra::Interval> res;

        for (const auto& [name, fractions] :
             coefficient_matrix | std::views::filter([](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& element) -> bool {
                 return element.first.variables[0].name[0] == 's';
             })) {
            algebra::Variable var("B" + std::to_string(i + 1));
            algebra::Fraction min = algebra::inf, max = -algebra::inf;

            for (int j = 0; j < size; j++) {
                if (fractions[j] > 0) {
                    max = std::max(max, -coefficient_matrix[LPP::B][j] / fractions[j]);
                } else if (fractions[j] < 0) {
                    min = std::min(min, -coefficient_matrix[LPP::B][j] / fractions[j]);
                }
            }
            res.push_back(max + lpp.constraints[i].rhs < var < min + lpp.constraints[i].rhs);
            i++;
        }
        return res;
    }

    void add_variable(const algebra::Variable& variable, const linalg::Matrix<algebra::Fraction>& coefficients) {
        int i = 0;
        const int size = basis_vector.size();
        linalg::Matrix<algebra::Fraction> res(size, size);
        lpp.objective += variable;

        for (const auto& [name, fractions] :
             coefficient_matrix | std::views::filter([](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& element) -> bool {
                 return element.first.variables[0].name[0] == 's';
             })) {
            for (int j = 0; j < size; j++) {
                res[i, j] = fractions[j];
            }
            i++;
        }
        res = coefficients * res;
        coefficient_matrix[variable.basis()] = res[0];
        cost[variable.basis()] = variable.coefficient;
        solution = Solution::UNOPTIMIZED;
    }

    void remove_variable(const algebra::Variable& variable) {
        if (std::ranges::contains(basis_vector, variable)) {
            cost[variable.basis()] = -LPP::M;
            solution = Solution::UNOPTIMIZED;
        } else {
            const int idx = std::ranges::distance(cost.begin(), cost.find(variable));
            lpp.objective -= cost[variable] * variable;
            zj_cj.erase(zj_cj.begin() + idx);
            cost.erase(variable);
            coefficient_matrix.erase(variable);
        }
    }

    void add_constraint(const algebra::Inequation& inequation) {
        std::vector<std::pair<std::string, algebra::Fraction>> substituent;
        const std::map<algebra::Variable, algebra::Fraction> temp =
            std::get<std::vector<std::map<algebra::Variable, algebra::Fraction>>>(get_solutions())[0];
        substituent.reserve(cost.size());

        for (const auto& [key, value] : temp) {
            substituent.emplace_back(key.variables[0].name, value);
        }
        if (static_cast<bool>(inequation.substitute(substituent))) {
            return;
        }

        auto range = coefficient_matrix | std::views::keys |
            std::views::filter([](const algebra::Variable& variable) -> bool { return variable.variables[0].name[0] == 's'; }) |
            std::views::transform([](const algebra::Variable& variable) -> int { return std::stoi(variable.variables[0].name.substr(1)); });
        const int size = basis_vector.size() + 1;
        const algebra::Variable slack("s" + std::to_string(*std::ranges::max_element(range) + 1));
        algebra::Polynomial transformation;
        substituent.resize(size);
        coefficient_matrix[LPP::B].push_back(static_cast<algebra::Fraction>(inequation.rhs));

        for (const algebra::Variable& variable : inequation.lhs.expression) {
            coefficient_matrix[variable.basis()].push_back(variable.coefficient);

            if (std::ranges::contains(basis_vector, variable.basis())) {
                transformation += variable;
            }
        }
        for (algebra::Variable& variable : transformation.expression) {
            const int idx = std::ranges::find(basis_vector, variable.basis()) - basis_vector.begin();
            const std::string name = "R" + std::to_string(idx + 1);
            variable = algebra::Variable(variable.coefficient, name);
            substituent[idx] = {name, 0};
        }
        for (std::vector<algebra::Fraction>& fractions : coefficient_matrix | std::views::values) {
            if (fractions.size() < size) {
                fractions.push_back(0);
            }
        }
        basis_vector.push_back(slack);
        cost.emplace(slack, 0);
        coefficient_matrix.emplace(slack, std::vector<algebra::Fraction>(size, 0));
        coefficient_matrix[slack].back() = 1;

        for (std::vector<algebra::Fraction>& fractions : coefficient_matrix | std::views::values) {
            for (int i = 0; i < size - 1; i++) {
                substituent[i].second = fractions[i];
            }
            fractions.back() -= static_cast<algebra::Fraction>(transformation.substitute(substituent));
        }
        solution = Solution::UNOPTIMIZED;
    }

    friend std::ostream& operator<<(std::ostream& out, const ComputationalTable& computational_table) {
        static constexpr int TAB_SIZE = 13;
        const int size = computational_table.basis_vector.size(), columns = 3 + computational_table.coefficient_matrix.size();
        auto print_partition = [columns, &out] -> void {
            out << std::right << std::setfill('-') << '+';

            for (int i = 0; i < columns; i++) {
                out << std::setw(TAB_SIZE) << "" << '+';
            }
            out << std::endl << std::left << std::setfill(' ');
        };
        auto format = [](const auto& value) -> std::string {
            const std::string res = std::to_string(value);
            const int remaining = TAB_SIZE - res.size();
            return std::string(remaining / 2, ' ') + res + std::string(remaining - remaining / 2, ' ');
        };
        out << std::left << ' ';

        for (int i = 0; i < 3; i++) {
            out << std::setw(TAB_SIZE) << "" << ' ';
        }
        for (const algebra::Variable& variable : computational_table.coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            out << format(computational_table.cost.at(variable)) << ' ';
        }
        out << std::endl;
        print_partition();
        out << '|';

        for (const std::string string : {"BV", "C", "B"}) {
            out << format(string) << '|';
        }
        for (const algebra::Variable& variable : computational_table.coefficient_matrix | std::views::drop(1) | std::views::keys) { // B
            out << format(variable.variables[0].name) << '|';
        }
        out << format(std::string("MR")) << '|' << std::endl;
        print_partition();

        for (int i = 0; i < size; i++) {
            out << '|' << format(computational_table.basis_vector[i]) << '|' << std::setw(TAB_SIZE)
                << format(computational_table.cost.at(computational_table.basis_vector[i])) << '|';

            for (const std::vector<algebra::Fraction>& fractions : computational_table.coefficient_matrix | std::views::values) {
                out << format(fractions[i]) << '|';
            }
            if (!computational_table.mr.empty()) {
                out << format(computational_table.mr[i] == algebra::inf ? "-ve" : std::to_string(computational_table.mr[i])) << '|';
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
        out << format(std::string("Zj-Cj")) << '|';

        for (const algebra::Polynomial& polynomial : computational_table.zj_cj) {
            out << format(polynomial) << '|';
        }
        out << std::setw(TAB_SIZE) << "" << '|' << std::endl;
        print_partition();
        return out << std::endl;
    }
};

inline optimization::LPP optimization::LPP::dual(const std::string& basis) const {
    LPP canonical = *this;

    for (algebra::Inequation& constraint : canonical.constraints) {
        if (type == Optimization::MAXIMIZE && constraint.opr == algebra::RelationalOperator::GE ||
            type == Optimization::MINIMIZE && constraint.opr == algebra::RelationalOperator::LE) {
            constraint = constraint.invert();
        }
    }
    const int objective_size = canonical.objective.expression.size(), constraints_size = canonical.constraints.size();
    LPP res;
    ComputationalTable computational_table(canonical);
    std::erase_if(computational_table.coefficient_matrix, [](const std::pair<algebra::Variable, std::vector<algebra::Fraction>>& element) -> bool {
        return element.first.variables[0].name[0] == 'A';
    });
    std::erase_if(computational_table.cost, [](const std::pair<algebra::Variable, algebra::Variable>& element) -> bool {
        return !element.second.variables.empty() && element.second.variables[0].name == "M";
    });
    res.type = canonical.type == Optimization::MAXIMIZE ? Optimization::MINIMIZE : Optimization::MAXIMIZE;
    res.constraints.resize(objective_size);
    res.restrictions.reserve(constraints_size);

    for (int i = 0; i < constraints_size; i++) {
        auto itr = std::next(computational_table.coefficient_matrix.begin()); // B
        algebra::Variable variable(basis + std::to_string(i + 1));
        res.objective += computational_table.coefficient_matrix[B][i] * variable;

        for (int j = 0; j < objective_size; j++) {
            res.constraints[j].lhs += itr->second[i] * variable;
            ++itr;
        }
        res.restrictions.push_back(canonical.constraints[i].opr == algebra::RelationalOperator::EQ ? unrestrict(variable) : variable >= 0);
    }
    auto itr = computational_table.cost.begin();

    for (int i = 0; i < objective_size; i++) {
        res.constraints[i].opr =
            canonical.restrictions[i].lhs.is_fraction() && static_cast<algebra::Fraction>(canonical.restrictions[i].lhs) == algebra::inf ||
                canonical.restrictions[i].rhs.is_fraction() && static_cast<algebra::Fraction>(canonical.restrictions[i].rhs) == algebra::inf
            ? algebra::RelationalOperator::EQ
            : res.type == Optimization::MAXIMIZE ? algebra::RelationalOperator::LE
                                                 : algebra::RelationalOperator::GE;
        res.constraints[i].rhs = itr->second;
        ++itr;
    }
    return res;
}

inline optimization::ComputationalTable optimization::LPP::tabular_optimize(const std::string& method) const {
    assert(method == "simplex" || method == "dual");
    const LPP lpp = method == "simplex" ? standardize() : canonicalize().standardize(true);

    if (GLOBAL_FORMATTING.verbose) {
        *GLOBAL_FORMATTING.out << std::endl << (method == "simplex" ? "Standard" : "Canonical") << " Form: " << std::endl << lpp;
    }
    return ComputationalTable(lpp);
}

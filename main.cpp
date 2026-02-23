#include "lpp.hpp"

using namespace lpp;

void test(LPP&& lpp, const std::string& method = "simplex", Variable var = Variable()) {
    if (method == "simplex" || method == "dual") {
        std::cout << lpp;
        const std::variant<std::vector<std::map<Variable, Fraction>>, std::string> res = lpp.optimize(method, true).get_solutions(method, true);

        if (auto ans = std::get_if<std::vector<std::map<Variable, Fraction>>>(&res)) {
            for (const std::map<Variable, Fraction>& i : *ans) {
                for (const auto& [variable, fraction] : i) {
                    std::cout << variable << '=' << fraction << " ";
                }
                std::cout << std::endl;
            }
        } else {
            std::cout << std::get<std::string>(res);
        }
    } else if (method.starts_with("Var")) {
        std::cout << lpp;
        lpp = lpp.standardize();
        ComputationalTable computational_table(lpp);
        computational_table.optimize_simplex(true);
        std::cout << computational_table;

        for (const Interval& interval : method.ends_with('C') ? computational_table.variation_C() : computational_table.variation_B()) {
            std::cout << interval << std::endl;
        }
    } else if (method.ends_with("add")) {
        std::cout << lpp;
        lpp = lpp.standardize();
        ComputationalTable computational_table(lpp);
        computational_table.optimize_simplex(true);
        std::cout << "BEFORE:" << std::endl << computational_table;
        // computational_table.add_variable(var);
        std::cout << "AFTER:" << std::endl << computational_table;
    } else {
        std::cout << lpp.dual(method);
    }
    std::cout << std::endl << std::endl;
}

void test(std::vector<std::map<Variable, Fraction>>&& res) {
    for (const std::map<Variable, Fraction>& result : res) {
        for (const auto& [variable, fraction] : result) {
            std::cout << variable << '=' << fraction << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
}

int main() {
    const Variable x("x"), y("y"), z("z"), x1("x1"), x2("x2"), x3("x3");

    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 x + y <= 4,
                 x - y <= 2,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, x1 + x2 + 3 * x3,
             {
                 3 * x1 + 2 * x2 + x3 <= 3,
                 2 * x1 + x2 + 2 * x3 <= 2,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MINIMIZE, x1 - 3 * x2 + 2 * x3,
             {
                 3 * x1 - x2 + 2 * x3 <= 7,
                 -2 * x1 + 4 * x2 <= 12,
                 -4 * x1 + 3 * x2 + 8 * x3 <= 10,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 4 * x + 3 * y <= 12,
                 4 * x + y <= 8,
                 4 * x - y <= 8,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 x - y <= 10,
                 2 * x - y <= 40,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 x - y <= 1,
                 3 * x - 2 * y <= 6,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, x1 + 2 * x2 + x3,
             {
                 2 * x1 + x2 - x3 >= -2,
                 -2 * x1 + x2 - 5 * x3 <= 6,
                 4 * x1 + x2 + x3 <= 6,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    // Big M Method
    test(LPP(Optimization::MAXIMIZE, -4 * x1 - x2,
             {
                 3 * x1 + x2 == 3,
                 4 * x1 + 3 * x2 >= 6,
                 x1 + 2 * x2 <= 3,
             },
             {x1 >= 0, x2 >= 0}));
    test(LPP(Optimization::MAXIMIZE, -x - y,
             {
                 3 * x + 2 * y >= 30,
                 -2 * x + 3 * y <= -30,
                 x + y <= 5,
             },
             {x >= 0, y >= 0}));
    // Two Phase Method
    test(LPP(Optimization::MAXIMIZE, -4 * x - y,
             {
                 3 * x + y == 3,
                 4 * x + 3 * y >= 6,
                 x + 2 * y <= 3,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, 3 * x1 + 2 * x2 + x3,
             {
                 -3 * x1 + 2 * x2 + 2 * x3 == 8,
                 -3 * x1 + 4 * x2 + x3 == 7,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    // BFS
    test(basic_feasible_solutions({
        x1 + 2 * x2 + x3 == 4,
        2 * x1 + x2 + 5 * x3 == 5,
    }));
    test(basic_feasible_solutions({
        2 * x1 + x2 - x3 == 2,
        3 * x1 + 2 * x2 + x3 == 3,
    }));
    test(LPP(Optimization::MAXIMIZE, 2 * x1 + 3 * x2 + 10 * x3,
             {
                 x1 + 2 * x3 == 0,
                 x2 + x3 == 1,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x1 + 3 * x2 + 10 * x3,
             {
                 x1 - 2 * x3 == 0,
                 x2 + x3 == 1,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 4 * x + 3 * y <= 12,
                 4 * x + y <= 8,
                 4 * x - y <= 8,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, 5 * x + 3 * y,
             {
                 x + y <= 2,
                 5 * x + 2 * y <= 10,
                 3 * x + 8 * y <= 12,
             },
             {x >= 0, y >= 0}));
    // Surprise Test
    test(LPP(Optimization::MAXIMIZE, 2 * x + 3 * y,
             {
                 x + 2 * y >= 2,
                 3 * x + y >= 3,
                 4 * x + 3 * y <= 6,
             },
             {x >= 0, y >= 0}));
    test(basic_feasible_solutions({
        x + 2 * y - z == 3,
        3 * x - y + 2 * z == 4,
        2 * x + 3 * y - 5 * z == 7,
    }));
    test(LPP(Optimization::MAXIMIZE, x1 + 2 * x2 + x3,
             {
                 2 * x1 + x2 - x3 >= -2,
                 -2 * x1 + x2 - 5 * x3 <= 6,
                 4 * x1 + x2 + x3 <= 6,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}));
    test(LPP(Optimization::MAXIMIZE, 2 * x + 3 * y,
             {
                 x + 2 * y <= 4,
                 x + y == 3,
             },
             {x >= 0, y >= 0}));
    // Dual
    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 4 * x - 3 * y <= 10,
                 x + 2 * y <= 5,
                 3 * x - 5 * y >= 15,
             },
             {x >= 0, y >= 0}),
         "w");
    test(LPP(Optimization::MINIMIZE, 5 * x1 + 4 * x2 - 3 * x3,
             {
                 x1 + x2 + x3 >= 5,
                 2 * x1 + 3 * x2 - 5 * x3 <= 4,
                 x1 + 2 * x2 - 3 * x3 <= 6,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "w");
    test(LPP(Optimization::MAXIMIZE, x + 2 * y,
             {
                 3 * x + 2 * y <= 4,
                 2 * x - 5 * y == 10,
             },
             {x >= 0, y >= 0}),
         "w");
    test(LPP(Optimization::MINIMIZE, 3 * x1 + 4 * x2 + 7 * x3,
             {
                 x1 + 2 * x2 - 5 * x3 >= 4,
                 2 * x1 + 5 * x2 - x3 == 7,
                 2 * x1 + 3 * x2 - x3 <= -8,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "w");
    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 x + y <= 5,
                 2 * x + 3 * y >= 4,
                 x - y <= 2,
             },
             {LPP::unrestrict(x), y >= 0}),
         "w");
    test(LPP(Optimization::MAXIMIZE, x + 2 * y,
             {
                 2 * x + 3 * y <= 4,
                 3 * x + 4 * y == 5,
             },
             {x >= 0, y >= 0}),
         "w");
    // Dual Simplex
    test(LPP(Optimization::MAXIMIZE, -5 * x - 6 * y,
             {
                 x + y >= 2,
                 4 * x + y >= 4,
             },
             {x >= 0, y >= 0}),
         "dual");
    test(LPP(Optimization::MINIMIZE, 10 * x1 + 6 * x2 + 2 * x3,
             {
                 -x1 + x2 + x3 >= 1,
                 3 * x1 + x2 - x3 >= 2,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "dual");
    // Alternate Optimal Solution
    test(LPP(Optimization::MAXIMIZE, 2 * x + 4 * y,
             {
                 x + 2 * y <= 5,
                 x + y <= 4,
             },
             {x >= 0, y >= 0}));
    test(LPP(Optimization::MAXIMIZE, -20 * x - 30 * y,
             {
                 2 * x + 3 * y >= 120,
                 x + y >= 40,
                 2 * x + 3 * y / 2 >= 90,
             },
             {x >= 0, y >= 0}));
    // Variation in C
    test(LPP(Optimization::MAXIMIZE, 3 * x + 5 * y,
             {
                 x + y <= 1,
                 2 * x + 3 * y <= 1,
             },
             {x >= 0, y >= 0}),
         "Var C");
    test(LPP(Optimization::MAXIMIZE, 15 * x + 45 * y,
             {
                 y <= 50,
                 x + 16 * y <= 240,
                 5 * x + 2 * y <= 162,
             },
             {x >= 0, y >= 0}),
         "Var C");
    test(LPP(Optimization::MAXIMIZE, x1 + x2 + 3 * x3,
             {
                 3 * x1 + 2 * x2 + x3 <= 3,
                 2 * x1 + x2 + 2 * x3 <= 2,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "Var C");
    // Variation in B
    test(LPP(Optimization::MAXIMIZE, -x1 + 2 * x2 - x3,
             {
                 3 * x1 + x2 - x3 <= 10,
                 -x1 + 4 * x2 + x3 >= 6,
                 x2 + x3 <= 4,
             },
             {x1 >= 0, x2 >= 0, x3 >= 0}),
         "Var B");
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 3 * x + 5 * y <= 15,
                 6 * x + 2 * y <= 24,
             },
             {x >= 0, y >= 0}),
         "Var B");
    test(LPP(Optimization::MAXIMIZE, 15 * x + 10 * y,
             {
                 4 * x + 6 * y <= 360,
                 3 * x <= 180,
                 5 * y <= 200,
             },
             {x >= 0, y >= 0}));
    // Addition of new variable
    test(LPP(Optimization::MAXIMIZE, 3 * x + 5 * y,
             {
                 x <= 4,
                 3 * x + 2 * y <= 18,
             },
             {x >= 0, y >= 0}), "Var add");
    return 0;
}

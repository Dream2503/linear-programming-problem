#include "lpp.hpp"

using namespace lpp;

void test(LPP&& lpp) {
    std::cout << lpp;
    const std::variant<std::map<Variable, Fraction>, std::string> res = lpp.optimize();

    if (auto ans = std::get_if<std::map<Variable, Fraction>>(&res)) {
        for (const auto& [variable, fraction] : *ans) {
            std::cout << variable << '=' << fraction << " ";
        }
    } else {
        std::cout << std::get<std::string>(res);
    }
    std::cout << std::endl << std::endl;
}

int main() {
    const Variable x("x"), y("y"), x1("x1"), x2("x2"), x3("x3");

    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 x + y <= 4,
                 x - y <= 2,
             }));
    test(LPP(Optimization::MAXIMIZE, x1 + x2 + 3 * x3,
             {
                 3 * x1 + 2 * x2 + x3 <= 3,
                 2 * x1 + x2 + 2 * x3 <= 2,
             }));
    test(LPP(Optimization::MINIMIZE, x1 - 3 * x2 + 2 * x3,
             {
                 3 * x1 - x2 + 2 * x3 <= 7,
                 -2 * x1 + 4 * x2 <= 12,
                 -4 * x1 + 3 * x2 + 8 * x3 <= 10,
             }));
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 4 * x + 3 * y <= 12,
                 4 * x + y <= 8,
                 4 * x - y <= 8,
             }));
    test(LPP(Optimization::MAXIMIZE, 2 * x + y,
             {
                 x - y <= 10,
                 2 * x - y <= 40,
             }));
    test(LPP(Optimization::MAXIMIZE, 3 * x + 2 * y,
             {
                 x - y <= 1,
                 3 * x - 2 * y <= 6,
             }));
    test(LPP(Optimization::MAXIMIZE, x1 + 2 * x2 + x3,
             {
                 2 * x1 + x2 - x3 >= -2,
                 -2 * x1 + x2 - 5 * x3 <= 6,
                 4 * x1 + x2 + x3 <= 6,
             }));
    // Big M Method
    test(LPP(Optimization::MAXIMIZE, -4 * x1 - x2,
             {
                 3 * x1 + x2 == 3,
                 4 * x1 + 3 * x2 >= 6,
                 x1 + 2 * x2 <= 3,
             }));
    test(LPP(Optimization::MAXIMIZE, -x - y,
             {
                 3 * x + 2 * y >= 30,
                 -2 * x + 3 * y <= -30,
                 x + y <= 5,
             }));
    // Two Phase Method
    test(LPP(Optimization::MAXIMIZE, -4 * x - y,
             {
                 3 * x + y == 3,
                 x + 2 * y <= 3,
             }));
    test(LPP(Optimization::MAXIMIZE, 3 * x1 + 2 * x2 + x3,
             {
                 -3 * x1 + 2 * x2 + 2 * x3 == 8,
                 -3 * x1 + 4 * x2 + x3 == 7,
             }));
}

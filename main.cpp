#include "lpp.hpp"

using namespace lpp;

void test(LPP&& lpp) {
    std::cout << lpp;
    const std::map<Variable, Fraction> res = lpp.optimize();

    for (const auto& [variable, fraction] : res) {
        std::cout << variable << '=' << fraction << " ";
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
                 x - y <= 10,
                 2 * x - y <= 40,
             }));
}

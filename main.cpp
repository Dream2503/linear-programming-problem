#include "lpp.hpp"

using namespace lpp;

int main() {
    const Variable x("x"), y("y");

    LPP lpp(Optimization::MAXIMIZE, 3 * x + 2 * y,
            {
                x + y <= 4,
                x - y <= 2,
            });
    std::cout << lpp;
    lpp = lpp.standardize();
    const auto [bv, c, coefficient_matrix] = lpp.prepare_computational_table();
    std::cout << "BV: ";

    for (const Variable& variable : bv) {
        std::cout << variable << " ";
    }
    std::cout << std::endl << "C: ";

    for (const auto& [first, second] : c) {
        std::cout << '(' << first << ',' << second << ')' << " ";
    }
    std::cout << std::endl << "Matrix: " << std::endl;

    for (const auto& [first, second] : coefficient_matrix) {
        std::cout << first << ": ";

        for (const Fraction fraction : second) {
            std::cout << fraction << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << lpp << std::endl;
    lpp.compute({bv, c, coefficient_matrix});
}

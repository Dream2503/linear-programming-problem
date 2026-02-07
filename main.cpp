#include "lpp.hpp"

using namespace lpp;

int main() {
    const Variable x("x"), y("y");

    LPP lpp(Optimization::MAXIMIZE, 2 * x + 3 * y,
            {
                5 * x - 4 * y <= 7,
                9 * x + 4 * y <= 9,
                x + y <= 4,
            });
    std::cout << lpp;
    const auto [bv, c, b, coefficient_matrix, mr] = lpp.prepare_computational_table();
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
    std::cout << std::endl << lpp;
}

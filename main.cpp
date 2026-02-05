#include "lpp.hpp"

int main() {
    const lpp::Variable x("x"), y("y");

    lpp::LPP lpp(lpp::OptimizeType::MAXIMIZE, 2 * x + 3 * y,
            {
                5 * x - 4 * y <= 7,
                9 * x + 4 * y <= 9,
                x + y <= 4,
            });
    // lpp.standardize();
    std::cout << lpp;
}

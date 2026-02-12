#include "algebra.hpp"

using namespace algebra;

int main() {
    Matrix<Fraction> matrix({{1, 2}, {2, 1}});
    const Matrix<Variable> variables({{Variable("X"), Variable("Y")}, {Variable("A"), Variable("B")}});

    // matrix *= matrix;
    // std::cout << matrix << std::endl;
    //
    // const auto res = matrix * variables;
    // std::cout << res << std::endl;

    matrix /= 10;
    std::cout << matrix << std::endl;

    const auto res1 = variables / Fraction(10);
    std::cout << res1 << std::endl;

    std::cout << matrix.inverse();

    Variable v = (Variable("x") ^ 1) * (Variable("y") ^ -1);
    std::cout << v;
    return 0;
}

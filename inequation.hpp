#pragma once

class Inequation : public Polynomial {
public:
    enum class Operator { LT, LE, GT, GE } opr;
    double rhs;

    Inequation(const std::vector<Variable>& terms, const Operator opr, const double rhs) : Polynomial(terms), opr(opr), rhs(rhs) {}

    friend std::ostream& operator<<(std::ostream& out, const Inequation& inequation) {
        out << static_cast<const Polynomial&>(inequation) << ' ';

        switch (inequation.opr) {
        case Operator::LT:
            out << '<';
            break;

        case Operator::LE:
            out << "<=";
            break;

        case Operator::GT:
            out << '>';
            break;

        case Operator::GE:
            out << ">=";
            break;
        }
        return out << ' ' << inequation.rhs;
    }
};

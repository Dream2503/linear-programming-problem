#pragma once

class lpp::Inequation {
public:
    enum class Operator { LT, LE, GT, GE, EQ } opr;
    Fraction rhs;
    Polynomial lhs;

    Inequation(const Polynomial& polynomial, const Operator opr, const Fraction& rhs) : opr(opr), rhs(rhs), lhs(polynomial) {}

    Inequation& operator*=(const Fraction& value) {
        lhs *= value;
        rhs *= value;

        if (value < 0) {
            switch (opr) {
            case Operator::LT:
                opr = Operator::GT;
                break;

            case Operator::LE:
                opr = Operator::GT;
                break;

            case Operator::GT:
                opr = Operator::LT;
                break;

            case Operator::GE:
                opr = Operator::LE;
                break;

            default:
                break;
            }
        }
        return *this;
    }

    Inequation operator*(const Fraction& value) const {
        Inequation inequation = *this;
        inequation *= value;
        return inequation;
    }

    friend std::ostream& operator<<(std::ostream& out, const Inequation& inequation) {
        out << inequation.lhs << ' ';

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

        case Operator::EQ:
            out << '=';
        }
        return out << ' ' << inequation.rhs;
    }
};

class lpp::Equation : public Inequation {
public:
    Equation(const Polynomial& polynomial, const Fraction& rhs) : Inequation(polynomial, Operator::EQ, rhs) {}
};

inline lpp::Inequation operator<(const lpp::Polynomial& polynomial, const lpp::Fraction& value) {
    return lpp::Inequation(polynomial, lpp::Inequation::Operator::LT, value);
}

inline lpp::Inequation operator<=(const lpp::Polynomial& polynomial, const lpp::Fraction& value) {
    return lpp::Inequation(polynomial, lpp::Inequation::Operator::LE, value);
}

inline lpp::Inequation operator>(const lpp::Polynomial& polynomial, const lpp::Fraction& value) {
    return lpp::Inequation(polynomial, lpp::Inequation::Operator::GT, value);
}

inline lpp::Inequation operator>=(const lpp::Polynomial& polynomial, const lpp::Fraction& value) {
    return lpp::Inequation(polynomial, lpp::Inequation::Operator::GE, value);
}

inline lpp::Equation operator==(const lpp::Polynomial& polynomial, const lpp::Fraction& value) { return lpp::Equation(polynomial, value); }

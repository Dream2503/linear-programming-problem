#pragma once

class lpp::Fraction {
    void simplify() {
        assert(denominator != 0);

        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
        if (numerator == 0) {
            denominator = 1;
        } else {
            const int gcd = std::gcd(numerator, denominator);
            numerator /= gcd;
            denominator /= gcd;
        }
    }

public:
    int numerator = 0, denominator = 1;

    Fraction(const int numerator = 0, const int denominator = 1) : numerator(numerator), denominator(denominator) { simplify(); }

    Fraction(double value) {
        while (static_cast<int>(value) != value) {
            denominator *= 10;
            value *= 10;
        }
        numerator = value;
        simplify();
    }

    Fraction operator-() const { return Fraction(-numerator, denominator); }

    Fraction& operator+=(const Fraction& value) {
        numerator = numerator * value.denominator + value.numerator * denominator;
        denominator *= value.denominator;
        simplify();
        return *this;
    }

    Fraction operator+(const Fraction& value) const {
        Fraction fraction = *this;
        fraction += value;
        return fraction;
    }

    Fraction& operator-=(const Fraction& value) {
        *this += -value;
        return *this;
    }

    Fraction operator-(const Fraction& value) const {
        Fraction fraction = *this;
        fraction -= value;
        return fraction;
    }

    Fraction& operator*=(const Fraction& value) {
        numerator *= value.numerator;
        denominator *= value.denominator;
        simplify();
        return *this;
    }

    Fraction operator*(const Fraction& value) const {
        Fraction fraction = *this;
        fraction *= value;
        return fraction;
    }

    Fraction& operator/=(const Fraction& value) {
        numerator *= value.denominator;
        denominator *= value.numerator;
        simplify();
        return *this;
    }

    Fraction operator/(const Fraction& value) const {
        Fraction fraction = *this;
        fraction /= value;
        return fraction;
    }

    Fraction& operator^=(const Fraction& value) {
        const double exponent = static_cast<double>(value);
        *this = Fraction(std::pow(numerator, exponent)) / Fraction(std::pow(denominator, exponent));
        return *this;
    }

    Fraction operator^(const Fraction& value) const {
        Fraction fraction = *this;
        fraction ^= value;
        return fraction;
    }

    std::strong_ordering operator<=>(const Fraction& value) const { return numerator * value.denominator <=> value.numerator * denominator; }
    std::partial_ordering operator<=>(const double value) const { return static_cast<double>(*this) <=> value; }

    bool operator==(const Fraction& value) const = default;
    bool operator==(const double value) const { return static_cast<double>(*this) == value; };

    explicit operator double() const { return static_cast<double>(numerator) / denominator; }

    friend std::ostream& operator<<(std::ostream& out, const Fraction& fraction) {
        out << fraction.numerator;

        if (fraction.denominator != 1) {
            out << '/' << fraction.denominator;
        }
        return out;
    }
};

namespace std {
    inline lpp::Fraction abs(lpp::Fraction fraction) {
        fraction.numerator = std::abs(fraction.numerator);
        return fraction;
    }
} // namespace std

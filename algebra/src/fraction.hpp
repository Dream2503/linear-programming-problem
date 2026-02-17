#pragma once

class algebra::Fraction {
    constexpr void simplify() {
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
    int numerator, denominator;

    constexpr Fraction(const int numerator = 0, const int denominator = 1) : numerator(numerator), denominator(denominator) { simplify(); }

    constexpr Fraction(double value) : numerator(0), denominator(1) {
        while (static_cast<int>(value) != value) {
            denominator *= 10;
            value *= 10;
        }
        numerator = value;
        simplify();
    }

    constexpr Fraction operator-() const { return Fraction(-numerator, denominator); }

    constexpr Fraction& operator+=(const Fraction& value) {
        numerator = numerator * value.denominator + value.numerator * denominator;
        denominator *= value.denominator;
        simplify();
        return *this;
    }

    constexpr Fraction operator+(const Fraction& value) const {
        Fraction fraction = *this;
        fraction += value;
        return fraction;
    }

    constexpr Fraction& operator-=(const Fraction& value) {
        *this += -value;
        return *this;
    }

    constexpr Fraction operator-(const Fraction& value) const {
        Fraction fraction = *this;
        fraction -= value;
        return fraction;
    }

    constexpr Fraction& operator*=(const Fraction& value) {
        numerator *= value.numerator;
        denominator *= value.denominator;
        simplify();
        return *this;
    }

    constexpr Fraction operator*(const Fraction& value) const {
        Fraction fraction = *this;
        fraction *= value;
        return fraction;
    }

    constexpr Fraction& operator/=(const Fraction& value) {
        numerator *= value.denominator;
        denominator *= value.numerator;
        simplify();
        return *this;
    }

    constexpr Fraction operator/(const Fraction& value) const {
        Fraction fraction = *this;
        fraction /= value;
        return fraction;
    }

    constexpr Fraction& operator^=(const Fraction& value) {
        const double exponent = static_cast<double>(value);
        *this = Fraction(std::pow(numerator, exponent)) / Fraction(std::pow(denominator, exponent));
        return *this;
    }

    constexpr Fraction operator^(const Fraction& value) const {
        Fraction fraction = *this;
        fraction ^= value;
        return fraction;
    }

    constexpr std::strong_ordering operator<=>(const Fraction& value) const {
        return numerator * value.denominator <=> value.numerator * denominator;
    }

    constexpr std::partial_ordering operator<=>(const double value) const { return static_cast<double>(*this) <=> value; }

    constexpr bool operator==(const Fraction& value) const = default;

    constexpr bool operator==(const double value) const { return static_cast<double>(*this) == value; }

    constexpr explicit operator double() const { return static_cast<double>(numerator) / denominator; }

    friend std::ostream& operator<<(std::ostream& out, const Fraction& fraction) {
        if (fraction.numerator == INT32_MAX && fraction.denominator == 1) {
            return out << "inf";
        }
        out << fraction.numerator;

        if (fraction.denominator != 1) {
            out << '/' << fraction.denominator;
        }
        return out;
    }
};

inline static constexpr algebra::Fraction inf = INT32_MAX;

namespace std {
    inline algebra::Fraction abs(algebra::Fraction fraction) {
        fraction.numerator = std::abs(fraction.numerator);
        return fraction;
    }
} // namespace std

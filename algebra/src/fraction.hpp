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

    constexpr Fraction& operator-=(const Fraction& value) { return *this += -value; }

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
        return *this = Fraction(std::pow(numerator, exponent)) / Fraction(std::pow(denominator, exponent));
    }

    constexpr Fraction operator^(const Fraction& value) const {
        Fraction fraction = *this;
        fraction ^= value;
        return fraction;
    }

    constexpr std::strong_ordering operator<=>(const Fraction& value) const {
        return static_cast<int64_t>(numerator) * value.denominator <=> static_cast<int64_t>(value.numerator) * denominator;
    }

    constexpr std::partial_ordering operator<=>(const double value) const { return static_cast<double>(*this) <=> value; }

    constexpr bool operator==(const Fraction& value) const = default;

    constexpr bool operator==(const double value) const { return static_cast<double>(*this) == value; }

    constexpr explicit operator double() const { return static_cast<double>(numerator) / denominator; }
};

inline static constexpr algebra::Fraction inf = INT32_MAX;

namespace std {
    inline algebra::Fraction abs(algebra::Fraction fraction) {
        fraction.numerator = abs(fraction.numerator);
        return fraction;
    }

    inline string to_string(const algebra::Fraction& fraction) {
        if (fraction == inf) {
            return "inf";
        }
        if (fraction == -inf) {
            return "-inf";
        }
        std::string res = std::to_string(fraction.numerator);

        if (fraction.denominator != 1) {
            res.push_back('/');
            res.append(to_string(fraction.denominator));
        }
        return res;
    }
} // namespace std

using namespace std;

inline std::ostream& algebra::operator<<(std::ostream& out, const Fraction& fraction) { return out << to_string(fraction); }

#ifndef DUALNUMBERS_CPP
#define DUALNUMBERS_CPP

#include "DualNumbers.hpp"
#include <cmath>

using namespace Math;

template <typename K>
const DualNumber<K> DualNumber<K>::epsilon = DualNumber(0, 1);

template <typename K>
DualNumber<K>::DualNumber() : a(0), b(0)
{

}

template <typename K>
DualNumber<K>::DualNumber(const K a) : a(a), b(0)
{
    
}

template <typename K>
DualNumber<K>::DualNumber(const K a, const K b) : a(a), b(b)
{
    
}

template <typename K>
DualNumber<K> Math::operator+(const DualNumber<K>& x, const DualNumber<K>& y)
{
    return DualNumber<K>(x.a + y.a, x.b + y.b);
}

template <typename K>
DualNumber<K> Math::operator+(const DualNumber<K>& x, const K y)
{
    return x + DualNumber<K>(y);
}

template <typename K>
DualNumber<K> Math::operator+(const K x, const DualNumber<K>& y)
{
    return DualNumber<K>(x) + y;
}

template <typename K>
DualNumber<K> Math::operator-(const DualNumber<K>& x, const DualNumber<K>& y)
{
    return DualNumber<K>(x.a - y.a, x.b - y.b);
}

template <typename K>
DualNumber<K> Math::operator-(const DualNumber<K>& x, const K y)
{
    return x - DualNumber<K>(y);
}

template <typename K>
DualNumber<K> Math::operator-(const K x, const DualNumber<K>& y)
{
    return DualNumber<K>(x) - y;
}

template <typename K>
DualNumber<K> Math::operator*(const DualNumber<K>& x, const DualNumber<K>& y)
{
    return DualNumber<K>(x.a * y.a, x.a * y.b + x.b * y.a);
}

template <typename K>
DualNumber<K> Math::operator*(const DualNumber<K>& x, const K y)
{
    return x * DualNumber<K>(y);
}

template <typename K>
DualNumber<K> Math::operator*(const K x, const DualNumber<K>& y)
{
    return DualNumber<K>(x) * y;
}

template <typename K>
DualNumber<K> Math::operator/(const DualNumber<K>& x, const DualNumber<K>& y)
{
    return DualNumber<K>(x.a / y.a, (x.b * y.a - x.a * y.b) / y.a / y.a);
}

template <typename K>
DualNumber<K> Math::operator/(const DualNumber<K>& x, const K y)
{
    return x / DualNumber<K>(y);
}

template <typename K>
DualNumber<K> Math::operator/(const K x, const DualNumber<K>& y)
{
    return DualNumber<K>(x) / y;
}

template <typename K>
std::string DualNumber<K>::to_string() const
{
    using std::to_string;
    if (b < 0)
    {
        return to_string(a) + " - " + to_string(-b) + "ε";
    }
    else
    {
        return to_string(a) + " + " + to_string(b) + "ε";
    }
}

template <typename K>
DualNumber<K> Math::cos(const DualNumber<K>& x)
{
    using ::cos;
    using ::sin;
    return DualNumber<K>(cos(x.a), - sin(x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::sin(const DualNumber<K>& x)
{
    using ::sin;
    using ::cos;
    return DualNumber<K>(sin(x.a), cos(x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::tan(const DualNumber<K>& x)
{
    using ::cos;
    using ::tan;
    // sec(a)
    K seca = 1./cos(x.a);
    return DualNumber<K>(tan(x.a), seca*seca * x.b);
}

template <typename K>
DualNumber<K> Math::acos(const DualNumber<K>& x)
{
    using ::acos;
    using ::sqrt;
    return DualNumber<K>(acos(x.a), - 1./sqrt(1 - x.a*x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::asin(const DualNumber<K>& x)
{
    using ::asin;
    using ::sqrt;
    return DualNumber<K>(asin(x.a), 1./sqrt(1 - x.a*x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::atan(const DualNumber<K>& x)
{
    using ::atan;
    return DualNumber<K>(atan(x.a), 1./(1 + x.a*x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::cosh(const DualNumber<K>& x)
{
    using ::cosh;
    using ::sinh;
    return DualNumber<K>(cosh(x.a), sinh(x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::sinh(const DualNumber<K>& x)
{
    using ::sinh;
    using ::cosh;
    return DualNumber<K>(sinh(x.a), cosh(x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::tanh(const DualNumber<K>& x)
{
    // sech(a)
    using ::cosh;
    using ::tanh;
    K secha = 1./cosh(x.a);
    return DualNumber<K>(tanh(x.a), secha*secha * x.b);
}

template <typename K>
DualNumber<K> Math::acosh(const DualNumber<K>& x)
{
    using ::acosh;
    using ::sqrt;
    return DualNumber<K>(acosh(x.a), 1. / sqrt(x.a*x.a - 1) * x.b);
}

template <typename K>
DualNumber<K> Math::asinh(const DualNumber<K>& x)
{
    using ::asinh;
    using ::sqrt;
    return DualNumber<K>(asinh(x.a), 1. / sqrt(x.a*x.a + 1) * x.b);
}

template <typename K>
DualNumber<K> Math::atanh(const DualNumber<K>& x)
{
    using ::atanh;
    return DualNumber<K>(atanh(x.a), 1. / (1 - x.a*x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::exp(const DualNumber<K>& x)
{
    using ::exp;
    // exp(a)
    K expa = ::exp(x.a);
    return DualNumber<K>(expa, expa * x.b);
}
template <typename K>
DualNumber<K> Math::log(const DualNumber<K>& x)
{
    using ::log;
    return DualNumber<K>(log(x.a), 1. / x.a * x.b);
}

template <typename K>
DualNumber<K> Math::log10(const DualNumber<K>& x)
{
    using ::log10;
    using ::log;
    return DualNumber<K>(log10(x.a), 1. / (x.a * log(10)) * x.b);
}

template <typename K>
DualNumber<K> Math::exp2(const DualNumber<K>& x)
{
    using ::exp2;
    using ::log;
    // exp2(a)
    K exp2a = exp2(x.a);
    return DualNumber<K>(exp2a, log(2) * exp2a * x.b);
}

template <typename K>
DualNumber<K> Math::expm1(const DualNumber<K>& x)
{
    using ::expm1;
    // expm1(a)
    K expm1a = expm1(x.a);
    return DualNumber<K>(expm1a, (1 + expm1a) * x.b);
}

template <typename K>
DualNumber<K> Math::log1p(const DualNumber<K>& x)
{
    using ::log1p;
    return DualNumber<K>(log1p(x.a), 1. / (1 + x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::log2(const DualNumber<K>& x)
{
    using ::log2;
    using ::log;
    return DualNumber<K>(log2(x.a), 1. / (x.a * log(2)) * x.b);
}

template <typename K>
DualNumber<K> Math::logb(const DualNumber<K>& x)
{
    using ::logb;
    using ::log;
    return DualNumber<K>(logb(x.a), 1. / (x.a * log(2)) * x.b);
}

template <typename K>
DualNumber<K> Math::pow(const DualNumber<K>& x, const DualNumber<K>& y)
{
    using ::pow;
    using ::log;
    K powxy = pow(x.a, y.a);
    return DualNumber<K>(powxy, powxy * (y.b*log(x.a) + y.a*x.b/x.a));
}

template <typename K>
DualNumber<K> Math::pow(const DualNumber<K>& x, const K& y)
{
    using ::pow;
    K powxy = pow(x.a, y);
    return DualNumber<K>(powxy, powxy * y*x.b/x.a);
}

template <typename K>
DualNumber<K> Math::pow(const K& x, const DualNumber<K>& y)
{
    using ::pow;
    using ::log;
    K powxy = pow(x, y.a);
    return DualNumber<K>(powxy, powxy * y.b*log(x));
}

template <typename K>
DualNumber<K> Math::sqrt(const DualNumber<K>& x)
{
    using ::sqrt;
    // sqrt(a)
    K sqrta = sqrt(x.a);
    return DualNumber<K>(sqrta, 1. / (2 * sqrta) * x.b);
}

template <typename K>
DualNumber<K> Math::cbrt(const DualNumber<K>& x)
{
    using ::cbrt;
    return DualNumber<K>(cbrt(x.a), 1. / (3 * cbrt(x.a*x.a)) * x.b);
}

template <typename K>
DualNumber<K> Math::erf(const DualNumber<K>& x)
{
    using ::erf;
    using ::sqrt;
    using ::exp;
    return DualNumber<K>(erf(x.a), 2./sqrt(M_PI)*exp(-x.a*x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::erfc(const DualNumber<K>& x)
{
    using ::erfc;
    using ::sqrt;
    using ::exp;
    return DualNumber<K>(erfc(x.a), -2./sqrt(M_PI)*exp(-x.a*x.a) * x.b);
}

template <typename K>
DualNumber<K> Math::fabs(const DualNumber<K>& x)
{
    using ::fabs;
    return DualNumber<K>(fabs(x.a), (x.a > 0) ? x.b : -x.b);
}


#endif // DUALNUMBERS_CPP
#ifndef RATIONAL_CPP
#define RATIONAL_CPP

#include "Rational.hpp"
using namespace Math;

template <typename K>
Rational<K>::Rational(): num(0), den(1)
{

}

template <typename K>
Rational<K>::Rational(const K num) : num(num), den(1)
{

}

template <typename K>
Rational<K>::Rational(K num, K den)
{
    if (den < 0)
    {
        den = -den;
        num = -num;
    }

    K gcd = std::gcd(num, den);
    this->num = num / gcd;
    this->den = den / gcd;
}

template <typename K>
Rational<K> Rational<K>::operator+ (const Rational<K>& other) const
{
    K lcd = std::lcm(this->den, other.den);
    return Rational<K>(lcd / this->den * this->num + lcd / other.den * other.num, lcd);
}

template <typename K>
Rational<K> Rational<K>::operator- (const Rational<K>& other) const
{
    K lcd = std::lcm(this->den, other.den);
    return Rational<K>(lcd / this->den * this->num - lcd / other.den * other.num, lcd);
}

template <typename K>
Rational<K> Rational<K>::operator* (const Rational<K>& other) const
{
    return Rational<K>(this->num * other.num, this->den * other.den);
}

template <typename K>
Rational<K> Rational<K>::operator/ (const Rational<K>& other) const
{
    return Rational<K>(this->num * other.den, this->den * other.num);
}

template <typename K>
Rational<K> Rational<K>::operator+= (const Rational<K>& other) const
{
    return *this + other;
}

template <typename K>
Rational<K> Rational<K>::operator-= (const Rational<K>& other) const
{
    return *this - other;
}

template <typename K>
Rational<K> Rational<K>::operator*= (const Rational<K>& other) const
{
    return *this * other;
}

template <typename K>
Rational<K> Rational<K>::operator/= (const Rational<K>& other) const
{
    return *this / other;
}

template <typename K>
bool Rational<K>::operator== (const Rational<K>& other) const
{
    return this->num * other.den == this->den * other.num;
}

template <typename K>
bool Rational<K>::operator< (const Rational<K>& other) const
{
    return this->num * other.den < this->den * other.num;
}

template <typename K>
std::string Rational<K>::toString() const
{
    if (this->den == 1)
    {
        return std::to_string(this->num);
    }
    else
    {
        return std::to_string(this->num) + "/" + std::to_string(this->den);
    }    
}

template <typename K>
std::string Rational<K>::toLatex() const
{
    if (this->den == 1)
    {
        return std::to_string(this->num);
    }
    else
    {
        return "\\frac{" + std::to_string(this->num) + "}{" + std::to_string(this->den) + "}";
    }    
}

#endif // RATIONAL_CPP
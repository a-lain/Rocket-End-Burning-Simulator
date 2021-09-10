#ifndef RATIONAL_HPP
#define RATIONAL_HPP

#include <string>


namespace Math
{
    // Forward declaration to make function declaration possible.
    template <typename K>
    class Rational;

    // And we now declare all friend operators.
    template <typename K>
    Rational<K>& operator+(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    Rational<K>& operator-(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    Rational<K>& operator*(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    Rational<K>& operator/(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    bool operator==(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    bool operator!=(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    bool operator<(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    bool operator<=(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    bool operator>(const Rational<K>& p, const Rational<K>& q);
    template <typename K>
    bool operator>=(const Rational<K>& p, const Rational<K>& q);


    template <typename K>
    class Rational
    {
        public:
        
        /*! Numerator. */
        K num;
        /*! Denominator. */
        K den;
        Rational();
        Rational(const K num);
        Rational(const K num, const K den);

        std::string to_string() const;

        friend Rational<K>& operator+<K>(const Rational<K>& p, const Rational<K>& q);
        friend Rational<K>& operator-<K>(const Rational<K>& p, const Rational<K>& q);
        friend Rational<K>& operator*<K>(const Rational<K>& p, const Rational<K>& q);
        friend Rational<K>& operator/<K>(const Rational<K>& p, const Rational<K>& q);
        friend bool operator==<K>(const Rational<K>& p, const Rational<K>& q);
        friend bool operator!=(const Rational<K>& p, const Rational<K>& q);
        friend bool operator<(const Rational<K>& p, const Rational<K>& q);
        friend bool operator<=(const Rational<K>& p, const Rational<K>& q);
        friend bool operator>(const Rational<K>& p, const Rational<K>& q);
        friend bool operator>=(const Rational<K>& p, const Rational<K>& q);      
    };
}


#endif // RATIONAL_HPP
/** 
 * @file DualNumbers.hpp
 * @brief This file implements dual numbers in one variable to allow for automatic differentiation.
 * @author Andrés Laín Sanclemente
 * @version 0.2.0
 * @date 9th September 2021 
 * 
 */


#ifndef DUALNUMBERS_HPP
#define DUALNUMBERS_HPP

#include <string>

namespace Math
{
    // We make a forward declaration of the class in order to be able to declare
    // friend operators.
    template <typename K>
    class DualNumber;

    // Forward declaration of friend operators.
    template <typename K>
    DualNumber<K> operator+(const DualNumber<K>& x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> operator+(const DualNumber<K>& x, const K y);
    template <typename K>
    DualNumber<K> operator+(const K x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> operator-(const DualNumber<K>& x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> operator-(const DualNumber<K>& x, const K y);
    template <typename K>
    DualNumber<K> operator-(const K x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> operator*(const DualNumber<K>& x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> operator*(const DualNumber<K>& x, const K y);
    template <typename K>
    DualNumber<K> operator*(const K x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> operator/(const DualNumber<K>& x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> operator/(const DualNumber<K>& x, const K y);
    template <typename K>
    DualNumber<K> operator/(const K x, const DualNumber<K>& y);

    /**
     * @brief Represents a dual number. A mathematical object written like
     * x = a + bε, where ε satisfies ε²=0.
     * 
     * Every mathematical function on the C++ math library has been extended
     * to act on DualNumbers. This allows automatic differentiation.
     * 
     * @tparam K 
     */
    template <typename K>
    class DualNumber
    {
        public:
        /**
         * @brief Real part.
         * 
         */
        K a;
        
        /**
         * @brief Infinitesimal part.
         * 
         */
        K b;

        /**
         * @brief The pure infinitesimal dual number.
         * 
         */
        static const DualNumber<K> epsilon;

        /**
         * @brief Constructs the dual number zero.
         * 
         */
        DualNumber();

        /**
         * @brief Constructs a dual number whose real part is @a a and whose infinitesimal part equal to zero.
         * 
         * @param a the real part.
         */
        DualNumber(const K a);

        /**
         * @brief Construct a dual number with @a a as real part and @a b as infinitesimal part. 
         * 
         * @param a 
         * @param b 
         */
        DualNumber(const K a, const K b);

        /**
         * @brief Returns the sum of two dual numbers.
         * x + y = (x.a + y.a) + (x.b + y.b)ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator+<K>(const DualNumber<K>& x, const DualNumber<K>& y);

        /**
         * @brief Returns the sum of a dual number and a real number.
         * x + y = (x.a + y) + x.b*ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator+<K>(const DualNumber<K>& x, const K y);

        /**
         * @brief Returns the sum of a dual number and a real number.
         * x + y = (x + y.a) + y.b*ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator+<K>(const K x, const DualNumber<K>& y);

        /**
         * @brief Returns the substraction of two dual numbers.
         * x - y = (x.a - y.a) + (x.b - y.b)ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator-<K>(const DualNumber<K>& x, const DualNumber<K>& y);

        /**
         * @brief Returns the substraction of a dual number and a real number.
         * x - y = (x.a - y) + x.b*ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator-<K>(const DualNumber<K>& x, const K y);

        /**
         * @brief Returns the substraction of a real number and a dual number.
         * x - y = (x - y.a) + y.b*ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator-<K>(const K x, const DualNumber<K>& y);

        /**
         * @brief Returns the product of two dual numbers.
         * x*y = (x.a*y.a) + (x.a*y.b + x.b*y.a)ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator*<K>(const DualNumber<K>& x, const DualNumber<K>& y);

        /**
         * @brief Returns the product of a dual number times a real number.
         * x*y = (x.a*y) + (x.b*y)ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator*<K>(const DualNumber<K>& x, const K y);

        /**
         * @brief Returns the product of a dual number times a real number.
         * x*y = (x*y.a) + (x*y.b)ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator*<K>(const K x, const DualNumber<K>& y);

        /**
         * @brief Performs the division of two dual numbers.
         * x/y = (x.a/y.a) + (x.b*y.a - x.a*y.b)/(y.a*y.a) * ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator/<K>(const DualNumber<K>& x, const DualNumber<K>& y);

        /**
         * @brief Returns the division of a dual number and a real number.
         * x/y = (x.a/y) + (x.b/y)ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator/<K>(const DualNumber<K>& x, const K y);

        /**
         * @brief Returns the division of a real number by a dual number.
         * x/y = (x/y.a) -x*y.b/(y.a*y.a) * ε
         * 
         * @param x 
         * @param y 
         * @return DualNumber<K> 
         */
        friend DualNumber<K> operator/<K>(const K x, const DualNumber<K>& y);

        /**
         * @brief Returns a string representation of the dual number: a + bε.
         * 
         * @return std::string 
         */
        std::string to_string() const;
    };

    /**
     * @brief Returns cos(x.a) - x.b*sin(x.a)ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */
    template <typename K>
    DualNumber<K> cos(const DualNumber<K>& x);

    /**
     * @brief Returns sin(x.a) + x.b*cos(x.a)ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */
    template <typename K>
    DualNumber<K> sin(const DualNumber<K>& x);

    /**
     * @brief Returns tan(x.a) + x.b/cos(x.a)/cos(x.a)*ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */
    template <typename K>
    DualNumber<K> tan(const DualNumber<K>& x);

    /**
     * @brief Returns acos(x.a) - 1/sqrt(1-x.a*x.a)*x.b*ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */
    template <typename K>
    DualNumber<K> acos(const DualNumber<K>& x);

    /**
     * @brief Returns asin(x.a) + 1/sqrt(1-x.a*x.a)*x.b*ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */
    template <typename K>
    DualNumber<K> asin(const DualNumber<K>& x);

    /**
     * @brief Returns atan(x.a) + 1/(1+x.a*x.a)*x.b*ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */
    template <typename K>
    DualNumber<K> atan(const DualNumber<K>& x);

    /**
     * @brief Returns cosh(x.a) + sinh(x.a)*x.b*ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */
    template <typename K>
    DualNumber<K> cosh(const DualNumber<K>& x);

    /**
     * @brief Returns sinh(x.a) + cosh(x.a)*x.b*ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */
    template <typename K>
    DualNumber<K> sinh(const DualNumber<K>& x);

    /**
     * @brief Returns tanh(x.a) + 1/cosh(x.a)/cosh(x.a)*x.b*ε
     * 
     * @tparam K 
     * @param x 
     * @return DualNumber<K> 
     */

    template <typename K>
    DualNumber<K> tanh(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> acosh(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> asinh(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> atanh(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> exp(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> log(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> log10(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> exp2(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> expm1(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> log1p(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> log2(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> logb(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> pow(const DualNumber<K>& x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> pow(const DualNumber<K>& x, const K& y);
    template <typename K>
    DualNumber<K> pow(const K& x, const DualNumber<K>& y);
    template <typename K>
    DualNumber<K> sqrt(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> cbrt(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> erf(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> erfc(const DualNumber<K>& x);
    /* To do. Not implemented. */
    template <typename K>
    DualNumber<K> tgamma(const DualNumber<K>& x);
    /* To do. Not implemented. */
    template <typename K>
    DualNumber<K> lgamma(const DualNumber<K>& x);
    template <typename K>
    DualNumber<K> fabs(const DualNumber<K>& x);
}

#endif // DUALNUMBERS_HPP
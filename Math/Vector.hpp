/** 
 * @file Vector.hpp
 * @brief This file contains all prototypes related to class Vector<K>. 
 * @author Andrés Laín Sanclemente
 * @version 0.9.0
 * @date 9th September 2021 
 * 
 */

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <string>
#include <complex>

/**
 * @brief Namespace that includes mathematical objects such as vectors, matrixes, dual numbers and useful numerical
 * methods such as numerical integrators, ode solvers and solvers for algebraic equations.
 * 
 */
namespace Math
{
    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return char 
     */
    char inline conj(const char& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned char 
     */
    unsigned char inline conj(const unsigned char& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return short int 
     */
    short int inline conj(const short int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned short int 
     */
    unsigned short int inline conj(const unsigned short int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return int 
     */
    int inline conj(const int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned int 
     */
    unsigned int inline conj(const unsigned int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return long int 
     */
    long int inline conj(const long int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned long int 
     */
    unsigned long int inline conj(const unsigned long int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return long long int 
     */
    long long int inline conj(const long long int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return unsigned long long int 
     */
    unsigned long long int inline conj(const unsigned long long int& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return float 
     */
    float inline conj(const float& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return double 
     */
    double inline conj(const double& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return long double 
     */
    long double inline conj(const long double& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns @a val .
     * 
     * @param val 
     * @return wchar_t 
     */
    wchar_t inline conj(const wchar_t& val)
    {
        return val;
    }

    /**
     * @brief This function exists to make the definition of the scalar product more general.
     * It just returns the complex conjugate of @a val .
     * 
     * @tparam K 
     * @param val 
     * @return std::complex<K> 
     */
    template<typename K>
    std::complex<K> inline conj(const std::complex<K>& val)
    {
        return std::conj(val);
    }
    
    template<typename K>
    std::string inline to_string(const std::complex<K>& z)
    {
        using std::to_string;
        return to_string(z.real()) + " + " + to_string(z.imag()) + "i"; 
    }

    // Forward declaration to make function declaration possible.
    template<typename K>
    class Vector;

    // And we now declare all friend operators.
    template<typename K>
    Vector<K> operator+(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    Vector<K> operator+(const Vector<K>& v, const K& alpha);
    template<typename K>
    Vector<K> operator+(const K& alpha, const Vector<K>& w);
    template<typename K>
    Vector<K> operator-(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    Vector<K> operator-(const Vector<K>& v, const K& alpha);
    template<typename K>
    Vector<K> operator-(const K& alpha, const Vector<K>& w);
    template<typename K>
    Vector<K> operator*(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    Vector<K> operator*(const K& alpha, const Vector<K>& v);
    template<typename K>
    Vector<K> operator*(const Vector<K>& v, const K& alpha);
    template<typename K>
    Vector<K> operator/(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    Vector<K> operator/(const Vector<K>& v, const K& alpha);
    template<typename K>
    Vector<K> operator/(const K& alpha, const Vector<K>& w);

    template<typename K>
    bool operator<(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    bool operator<(const Vector<K>& v, const K& alpha);
    template<typename K>
    bool operator<(const K& alpha, const Vector<K>& v); 
    template<typename K>
    bool operator<=(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    bool operator<=(const Vector<K>& v, const K& alpha);
    template<typename K>
    bool operator<=(const K& alpha, const Vector<K>& v);
    template<typename K>
    bool operator>(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    bool operator>(const Vector<K>& v, const K& alpha);
    template<typename K>
    bool operator>(const K& alpha, const Vector<K>& v);
    template<typename K>
    bool operator>=(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    bool operator>=(const Vector<K>& v, const K& alpha);
    template<typename K>
    bool operator>=(const K& alpha, const Vector<K>& v);
    template<typename K>
    bool operator==(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    bool operator==(const Vector<K>& v, const K& alpha);
    template<typename K>
    bool operator==(const K& alpha, const Vector<K>& v);
    template<typename K>
    bool operator!=(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    bool operator!=(const Vector<K>& v, const K& alpha);
    template<typename K>
    bool operator!=(const K& alpha, const Vector<K>& v);

    template<typename K>
    Vector<K> operator&(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    Vector<K> operator&(const Vector<K>& v, const K& w);
    template<typename K>
    Vector<K> operator&(const K& v, const Vector<K>& w);

    template<typename K>
    K operator|(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    K vector_product_2d(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    Vector<K> vector_product_3d(const Vector<K>& v, const Vector<K>& w);

    template<typename K>
    std::ostream& operator<< (std::ostream& os, const Vector<K>& v);

    template<typename K>
    K min(const Vector<K>& v);
    template<typename K>
    K max(const Vector<K>& v);
    template<typename K>
    K sum(const Vector<K>& v);
    template<typename K>
    K multiply(const Vector<K>& v);
    template<typename K>
    Vector<K> cos(const Vector<K>& v);
    template<typename K>
    Vector<K> sin(const Vector<K>& v);
    template<typename K>
    Vector<K> tan(const Vector<K>& v);
    template<typename K>
    Vector<K> acos(const Vector<K>& v);
    template<typename K>
    Vector<K> asin(const Vector<K>& v);
    template<typename K>
    Vector<K> atan(const Vector<K>& v);
    template<typename K>
    Vector<K> atan2(const K v, const Vector<K>& w);
    template<typename K>
    Vector<K> atan2(const Vector<K>& v, const K w);
    template<typename K>
    Vector<K> atan2(const Vector<K>& v, const Vector<K>& w);
    template<typename K>
    Vector<K> cosh(const Vector<K>& v);
    template<typename K>
    Vector<K> sinh(const Vector<K>& v);
    template<typename K>
    Vector<K> tanh(const Vector<K>& v);
    template<typename K>
    Vector<K> acosh(const Vector<K>& v);
    template<typename K>
    Vector<K> asinh(const Vector<K>& v);
    template<typename K>
    Vector<K> atanh(const Vector<K>& v);
    template<typename K>
    Vector<K> exp(const Vector<K>& v);
    template<typename K>
    Vector<K> frexp(const Vector<K>& v, Vector<int>* exp);
    template<typename K>
    Vector<K> ldexp(const Vector<K>& v, const int exp);
    template<typename K>
    Vector<K> ldexp(const Vector<K>& v, const Vector<int>& exp);
    template<typename K>
    Vector<K> log(const Vector<K>& v);
    template<typename K>
    Vector<K> log10(const Vector<K>& v);
    template<typename K>
    Vector<K> modf(const Vector<K>& v, Vector<K>* intpart);
    template<typename K>
    Vector<K> exp2(const Vector<K>& v);
    template<typename K>
    Vector<K> expm1(const Vector<K>& v);
    template<typename K>
    Vector<K> ilogb(const Vector<K>& v);
    template<typename K>
    Vector<K> log1p(const Vector<K>& v);
    template<typename K>
    Vector<K> log2(const Vector<K>& v);
    template<typename K>
    Vector<K> logb(const Vector<K>& v);
    template<typename K>
    Vector<K> scalbn(const Vector<K>& v, const int n);
    template<typename K>
    Vector<K> scalbn(const Vector<K>& v, const Vector<int>& n);
    template<typename K>
    Vector<K> scalbln(const Vector<K>& v, const long int n);
    template<typename K>
    Vector<K> scalbln(const Vector<K>& v, const Vector<long int>& n);
    template<typename K>
    Vector<K> pow(const K v, const Vector<K>& exponent);
    template<typename K>
    Vector<K> pow(const Vector<K>& v, const K exponent);
    template<typename K>
    Vector<K> pow(const Vector<K>& v, const Vector<K>& exponent);
    template<typename K>
    Vector<K> sqrt(const Vector<K>& v);
    template<typename K>
    Vector<K> cbrt(const Vector<K>& v);
    template<typename K>
    Vector<K> hypot(const K x, const Vector<K>& y);
    template<typename K>
    Vector<K> hypot(const Vector<K>& x, const K y);
    template<typename K>
    Vector<K> hypot(const Vector<K>& x, const Vector<K>& y);
    template<typename K>
    Vector<K> erf(const Vector<K>& v);
    template<typename K>
    Vector<K> erfc(const Vector<K>& v);
    template<typename K>
    Vector<K> tgamma(const Vector<K>& v);
    template<typename K>
    Vector<K> lgamma(const Vector<K>& v);
    template<typename K>
    Vector<K> ceil(const Vector<K>& v);
    template<typename K>
    Vector<K> floor(const Vector<K>& v);
    template<typename K>
    Vector<K> fmod(const K numer, const Vector<K>& denom);
    template<typename K>
    Vector<K> fmod(const Vector<K>& numer, const K denom);
    template<typename K>
    Vector<K> fmod(const Vector<K>& numer, const Vector<K>& denom);
    template<typename K>
    Vector<K> trunc(const Vector<K>& v);
    template<typename K>
    Vector<K> round(const Vector<K>& v);
    template<typename K>
    Vector<long int> lround(const Vector<K>& v);
    template<typename K>
    Vector<long long int> llround(const Vector<K>& v);
    template<typename K>
    Vector<K> rint(const Vector<K>& v);
    template<typename K>
    Vector<long int> lrint(const Vector<K>& v);
    template<typename K>
    Vector<long long int> llrint(const Vector<K>& v);
    template<typename K>
    Vector<K> nearbyint(const Vector<K>& v);
    template<typename K>
    Vector<K> remainder(const K numer, const Vector<K>& denom);
    template<typename K>
    Vector<K> remainder(const Vector<K>& numer, const K denom);
    template<typename K>
    Vector<K> remainder(const Vector<K>& numer, const Vector<K>& denom);
    template<typename K>
    Vector<K> remquo(const K numer, const Vector<K>& denom, Vector<int>* quot);
    template<typename K>
    Vector<K> remquo(const Vector<K>& numer, const K denom, Vector<int>* quot);
    template<typename K>
    Vector<K> remquo(const Vector<K>& numer, const Vector<K>& denom, Vector<int>* quot);
    template<typename K>
    Vector<K> copysign(const Vector<K>& x, const K y);
    template<typename K>
    Vector<K> copysign(const Vector<K>& x, const Vector<K>& y);
    template<typename K>
    Vector<K> nan(const unsigned int N, const char* tagp);
    template<typename K>
    Vector<K> nextafter(const Vector<K>& x, const Vector<K>& y);
    template<typename K>
    Vector<K> fdim(const K x, const Vector<K>& y);
    template<typename K>
    Vector<K> fdim(const Vector<K>& x, K y);
    template<typename K>
    Vector<K> fdim(const Vector<K>& x, const Vector<K>& y);
    template<typename K>
    Vector<double> fabs(const Vector<K>& v);
    template<typename K>
    Vector<double> abs(const Vector<K>& v);
    template<typename K>
    Vector<K> fma(const K x, const Vector<K>& y, const Vector<K>& z);
    template<typename K>
    Vector<K> fma(const Vector<K>& x, const K y, const Vector<K>& z);
    template<typename K>
    Vector<K> fma(const Vector<K>& x, const Vector<K>& y, const K z);
    template<typename K>
    Vector<K> fma(const K x, const K y, const Vector<K>& z);
    template<typename K>
    Vector<K> fma(const K x, const Vector<K>& y, const K z);
    template<typename K>
    Vector<K> fma(const Vector<K>& x, const K y, const K z);
    template<typename K>
    Vector<K> fma(const Vector<K>& x, const Vector<K>& y, const Vector<K>& z);
    template<typename K>
    Vector<int> fpclassify(const Vector<K>& v);
    template<typename K>
    bool isfinite(const Vector<K>& v);
    template<typename K>
    bool isinf(const Vector<K>& v);
    template<typename K>
    bool isnan(const Vector<K>& v);
    template<typename K>
    bool isnormal(const Vector<K>& v);

    template<typename K>
    Vector<K> real(const Vector<std::complex<K>>& v);
    template<typename K>
    Vector<K> imag(const Vector<std::complex<K>>& v);
    template<typename K>
    Vector<K> arg(const Vector<std::complex<K>>& v);
    template<typename K>
    Vector<K> conj(const Vector<K>& v);


    /** 
     * @brief Generic Vector over the field K with any number of components.
     *
     * It is implemented through generic programming with a template typename K. The
     * components of the Vector are stored in an array which is dynamically assigned.
     *
     * It implements the following Vector componentwise operations: +, -, *, /.
     * Boolean operators <, <=, >, >=, ==, != apply the corresponding operator
     * component by component and "and" the results together. In other words, two
     * vectors @a v and @a w satisfy <em>v < w = true</em> if and only if
     * <em>a[i] < b[i]</em> for all @a i . Otherwise, <em>a < b = false</em>.
     * In addition, all previously mentioned operators allow a Vector-scalar or
     * a scalar-Vector version, that is: <em>v + 5</em> or <em>5 + v</em>, for
     * example. In that case, the scalar operator is applied to all components of
     * the Vector, in other words: <em>v + 5</em> is equivalent to <em>v[i] + 5
     * </em> for all @a i . In the same way, <em>v < 5 = true</em> is the same as
     * <em>v[i] < 5 = true</em> for all @a i".
     *
     * Moreover, a scalar product operator | is provided. Note that C++ operator
     * precedence for | is low; so it is recommend to always write parenthesis that
     * englobe the scalar product, obtaining a mixture of Dirac's and
     * mathematicians' notation for the scalar product: <em>scalar product of @a v and @a w
     * </em> = <em>(v|w)</em>. The object also counts with a concatenation operator &:
     * <em>v & w</em> places the components of @a v followed by the components of @a w 
     * in a new Vector.
     *
     * The index operator [] can be used to access the Vector components. Bear in mind
     * that indexes start with zero, following the convention of the C programming
     * language.
     * @tparam K Any mathematical field. Normally @a double or @a float. 
     * 
     */
    template<typename K>
    class Vector
    {
        public:
        /**
         * @brief Array where the components of the Vector are stored. 
         * 
         */
        K* components;

        /**
         * @brief Dimension (also called length) of the Vector.
         * 
         */
        unsigned int N;

        public:
        // Constructors and destructor.

        /**
         * @brief Creates a Vector of length zero. A Vector of length zero can be
         * safely added or substracted to another Vector of any length.
         * 
         */
        Vector<K>();

        /**
         * @brief Creates a Vector of length @a N .
         * 
         * @param N Length of the Vector.
         * 
         */
        explicit Vector<K>(const unsigned int N);

        /**
         * @brief Copy constructor. Creates a deep copy of @a w .
         * 
         * @param w 
         */
        Vector<K>(const Vector<K>& w);

        /**
         * @brief Move constructor. Copies the pointer @a components and @a N .
         * 
         * @warning It leaves @a w in a corrupted state.
         * 
         * @param w 
         */
        Vector<K>(Vector<K>&& w);
        
        /**
         * @brief Creates a Vector of length @a N whose components are copied from
         * the array @a components.
         * 
         * @param components An array of numbers.
         * @param N The length of the array.
         * @warning May cause undefined behaviour if the array supplied has not
         * been properly constructed.
         * 
         */
        Vector<K>(const K* components, const unsigned int N);

        /**
         * @brief Allows to create a Vector from a initializer_list. The dimension of
         * the Vector is automatically deduced.
         * 
         * @param list 
         */
        Vector<K>(std::initializer_list<K> list);

        /**
         * @brief Frees the contents of @a components .
         * 
         */
        ~Vector<K>();

        // Assignment operators

        /**
         * @brief Equivalent to v[i] = alpha for all i.
         * 
         * @param alpha 
         * @return Vector<K>& 
         */
        Vector<K>& operator=(const K& alpha);

        /**
         * @brief Assigment operator.
         * 
         * If needed, the @a components array is freed
         * and dynamically reassigned to match the dimension of @a w. Afterwards,
         * @a N and the components of the Vector @a w are copied to the current
         * object.
         * 
         * @param w 
         * @return Vector<K>& 
         */
        Vector<K>& operator=(const Vector<K>& w);

        /**
         * @brief Move-assigment operator.
         * 
         * Copies @a components and @a N instead of
         * creating a new array. Faster than normal assignment operator.
         * 
         * @warning It leaves @a w in a corrupted state.
         * 
         * @param w 
         * @return Vector<K>& 
         */
        Vector<K>& operator=(Vector<K>&& w);

        // Operators defined componentwise

        /**
         * @brief Componentwise addition operator.
         * 
         * <em>v + w</em> is equivalent to <em>v[i] + w[i]</em> for all @a i .
         * If any of the vectors has length zero, then the other one is returned.
         * 
         * @exception A runtime exception is thrown if both vectors have
         * different number of components and their lengths are both non-zero.
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator+<K>(const Vector<K>& v, const Vector<K>& w);

        /**
         * @brief Componentwise Vector + scalar addition operator.
         *  
         * <em>v + alpha</em> is equivalent to <em>v[i] + alpha</em> for all @a i .
         * 
         * @param v 
         * @param alpha 
         * @return Vector<K> 
         */
        friend Vector<K> operator+<K>(const Vector<K>& v, const K& alpha);

        /**
         * @brief Component wise scalar + Vector addition operator.
         *  
         * <em>alpha + v</em> is equivalent to <em>alpha + v[i]</em> for all @a i .
         * 
         * @param alpha 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator+<K>(const K& alpha, const Vector<K>& w);

        /**
         * @brief Componentwise substraction operator.
         * 
         * <em>v - w</em> is equivalent to <em>v[i] - w[i]</em> for all @a i . 
         * If @a v has null length, then @a -w is returned. If @a w has length
         * zero, then @a v is returned.
         * 
         * @exception A runtime exception is thrown if both vectors have
         * different number of components and their lengths are both non-zero.
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator-<K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Componentwise Vector - scalar substraction operator. 
         * 
         * <em>v - alpha</em> is equivalent to <em>v[i] - alpha</em> for all @a i.
         * 
         * @param v 
         * @param alpha 
         * @return Vector<K> 
         */
        friend Vector<K> operator-<K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Componentwise scalar - Vector substraction operator. 
         * 
         * <em>alpha - v</em> is equivalent to <em>alpha - v[i]</em> for all @a i .
         * 
         * @param alpha 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator-<K>(const K& alpha, const Vector<K>& w);
        
        /**
         * @brief Component wise multiplication operator.
         * 
         * <em>v * w</em> is equivalent to <em>v[i] * w[i]</em> for all @a i .
         * 
         * @exception An exception is thrown if @a v and @a w have a different number
         * of components.
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator*<K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Component wise Vector - scalar multiplication operator. 
         * 
         * <em>v * alpha</em> is equivalent to <em>v[i] * alpha</em> for all @a i .
         * 
         * @param alpha 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> operator*<K>(const K& alpha, const Vector<K>& v);
        
        /**
         * @brief Component wise scalar - Vector multiplication operator.
         * 
         * <em>alpha * v</em> is equivalent to <em>alpha * v[i]</em> for all @a i .
         * 
         * @param v 
         * @param alpha 
         * @return Vector<K> 
         */
        friend Vector<K> operator*<K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Componentwise division operator.
         * 
         * <em>v / w</em> is equivalent to <em>v[i] / w[i]</em> for all @a i .
         * 
         * @exception An exception is thrown if @a v and @a w have a different number
         * of components.
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator/<K>(const Vector<K>& v, const Vector<K>& w);

        /**
         * @brief Componentwise Vector - scalar division operator. 
         * 
         * <em>v / alpha</em> is equivalent to <em>v[i] / alpha</em> for all @a i .
         * 
         * @param v 
         * @param alpha 
         * @return Vector<K> 
         */
        friend Vector<K> operator/<K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Componentwise scalar - Vector division operator. 
         * 
         * <em>alpha / v</em> is equivalent to <em>alpha / v[i]</em> for all @a i .
         * 
         * @param alpha 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator/<K>(const K& alpha, const Vector<K>& w);

        // The following operators are done component by component and && together.
        /**
         * @brief Returns true if <em>v[i] \< w[i]</em> for all @a i and false otherwise.
         * 
         * @exception An exception is raised if @a v and @a w have a different number
         * of components.
         * 
         * @param v 
         * @param w 
         * @return true 
         * @return false 
         */
        friend bool operator< <K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Returns @a true if <em>v[i] \< alpha</em> for all @a i and @a false otherwise.
         * 
         * @param v 
         * @param alpha 
         * @return true 
         * @return false 
         */
        friend bool operator< <K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Returns @a true if <em>alpha \< v[i]</em> for all @a i and @a false otherwise.
         * 
         * @param alpha 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool operator< <K>(const K& alpha, const Vector<K>& v);

        /**
         * @brief Returns @a true if <em>v[i] \<= w[i]</em> for all @a i and @a false otherwise.
         * 
         * @exception An exception is raised if @a v and @a w have a different number
         * of components.
         * 
         * @param v 
         * @param w 
         * @return true 
         * @return false 
         */
        friend bool operator<= <K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Returns @a true if <em>v[i] \<= alpha</em> for all @a i and @a false otherwise.
         * 
         * @param v 
         * @param alpha 
         * @return true 
         * @return false 
         */
        friend bool operator<= <K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Returns @a true if <em>alpha \<= v[i]</em> for all @a i and @a false otherwise.
         * 
         * @param alpha 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool operator<= <K>(const K& alpha, const Vector<K>& v);
        
        /**
         * @brief Returns @a true if <em>v[i] \> w[i]</em> for all @a i and @a false otherwise.
         * 
         * @exception An exception is raised if @a v and @a w have a different number
         * of components.
         * 
         * @param v 
         * @param w 
         * @return true 
         * @return false 
         */
        friend bool operator> <K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Returns @a true if <em>v[i] \> alpha</em> for all @a i and @a false otherwise.
         * 
         * @param v 
         * @param alpha 
         * @return true 
         * @return false 
         */
        friend bool operator> <K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Returns @a true if <em>alpha \> v[i]</em> for all @a i and @a false otherwise.
         * 
         * @param alpha 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool operator> <K>(const K& alpha, const Vector<K>& v);

        /**
         * @brief Returns @a true if <em>v[i] \>= w[i]</em> for all @a i and @a false otherwise.
         * 
         * @exception An exception is raised if @a v and @a w have a different number
         * of components.
         * 
         * @param v 
         * @param w 
         * @return true 
         * @return false 
         */
        friend bool operator>= <K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Returns @a true if <em>v[i] \>= alpha</em> for all @a i and @a false otherwise.
         * 
         * @param v 
         * @param alpha 
         * @return true 
         * @return false 
         */
        friend bool operator>= <K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Returns @a true if <em>alpha \>= v[i]</em> for all @a i and @a false otherwise.
         * 
         * @param alpha 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool operator>= <K>(const K& alpha, const Vector<K>& v);
        
        /**
         * @brief If @a v and @a w have a different number of componentes, @a false is returned.
         * If @a v and @a w have the same number of componenets, then @a true is returned if
         * <em>v[i] == w[i]</em> for all @a i and @a false is returned otherwise.
         * 
         * @param v 
         * @param w 
         * @return true 
         * @return false 
         */
        friend bool operator== <K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Returns @a true if <em>v[i] == alpha</em> for all @a i and @a false otherwise.
         * 
         * @param v 
         * @param alpha 
         * @return true 
         * @return false 
         */
        friend bool operator== <K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Returns @a true if <em>alpha == v[i]</em> for all @a i and @a false otherwise.
         * 
         * @param alpha 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool operator== <K>(const K& alpha, const Vector<K>& v);
        
        /**
         * @brief Equivalent to <em>!(v == w)</em> .
         * 
         * @param v 
         * @param w 
         * @return true 
         * @return false 
         */
        friend bool operator!= <K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Returns @a true if <em>v[i] != alpha</em> for all @a i and @a false otherwise.
         * 
         * @param v 
         * @param alpha 
         * @return true 
         * @return false 
         */
        friend bool operator!= <K>(const Vector<K>& v, const K& alpha);
        
        /**
         * @brief Returns @a true if <em>alpha != v[i]</em> for all @a i and @a false otherwise.
         * 
         * @param alpha 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool operator!= <K>(const K& alpha, const Vector<K>& v);

        /**
         * @brief Concatenation operator.
         * 
         * Creates a new Vector whose components
         * are the components of @a v followed by the components of @a w .
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator&<K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Concatenation operator.
         * 
         * Creates a new Vector whose components
         * are the components of @a v followed by the scalar @a w .
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator&<K>(const Vector<K>& v, const K& w);
        
        /**
         * @brief Concatenation operator.
         * 
         * Creates a new Vector whose components
         * are the scalar @a v followed by the components of @a w .
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> operator&<K>(const K& v, const Vector<K>& w);

        /**
         * @brief Scalar product operator. 
         * 
         * Performs the following operation:
         * \f[ (v \vert w) = \sum_{i=0}^{N-1} v_i^{*} w_i\f]
         * which is the canonical scalar product of \f$\mathbb{K}^N\f$.
         * 
         * @param v 
         * @param w 
         * @return K 
         */
        friend K operator|<K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Returns the Vector product of two 2d vectors @a v and @a w.
         * 
         * @exception Throws an exception if either the dimension of @a v or
         * the dimension of @a w is different from two.
         * 
         * @param v 
         * @param w 
         * @return K 
         */
        friend K vector_product_2d<K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Returns the Vector product of two 3d vectors @a v and @a w.
         * 
         * @exception Throws an exception if either the dimension of @a v or
         * the dimension of @a w is different from three.
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> vector_product_3d<K>(const Vector<K>& v, const Vector<K>& w);

        // Unary operators: +, -
        /**
         * @brief Returns the vector unaltered.
         * 
         * @return Vector<K> 
         */
        Vector<K> operator+();

        /**
         * @brief Returns minus the vector.
         * 
         * @return Vector<K> 
         */
        Vector<K> operator-();

        // Basic math-assignment operators: +=, -=, *=, /=
        /**
         * @brief Equivalent to <em>v[i] += alpha</em> for all @a i .
         * 
         * @param alpha 
         * @return Vector<K>& 
         */
        Vector<K>& operator+=(const K& alpha);
        
        /**
         * @brief If @a this has length zero, than @a w is copied into @a this .
         * If @a this and @a w have the same dimension, this is equivalent
         * to <em>this->operator[](i) += w[i]</em> for all @a i . Otherwise,
         * an exception is raised.
         * 
         * @exception An exception is thrown if @a this and @a w have a different
         * number of components and the dimension of @a this is non-zero.
         * 
         * @param w 
         * @return Vector<K>& 
         */
        Vector<K>& operator+=(const Vector<K>& w);
        
        /**
         * @brief Equivalent to <em>v[i] -= alpha</em> for all @a i .
         * 
         * @param alpha 
         * @return Vector<K>& 
         */
        Vector<K>& operator-=(const K& alpha);
        
        /**
         * @brief If @a this has length zero, than @a -w is copied into @a this .
         * If @a this and @a w have the same dimension, this is equivalent
         * to <em>this->operator[](i) -= w[i]</em> for all @a i . Otherwise,
         * an exception is raised.
         * 
         * @exception An exception is thrown if @a this and @a w have a different
         * number of components and the dimension of @a this is non-zero.
         * 
         * @param w 
         * @return Vector<K>& 
         */
        Vector<K>& operator-=(const Vector<K>& w);
        
        /**
         * @brief Equivalent to <em>v[i] *= alpha</em> for all @a i .
         * 
         * @param alpha 
         * @return Vector<K>& 
         */
        Vector<K>& operator*=(const K& alpha);
        
        /**
         * @brief If @a this and @a w have the same dimension, this is equivalent
         * to <em>this->operator[](i) *= w[i]</em> for all @a i . Otherwise,
         * an exception is raised.
         * 
         * @exception An exception is thrown if @a this and @a w have a different
         * number of components.
         * 
         * @param w 
         * @return Vector<K>& 
         */
        Vector<K>& operator*=(const Vector<K>& w);
        
        /**
         * @brief Equivalent to <em>v[i] /= alpha</em> for all @a i .
         * 
         * @param alpha 
         * @return Vector<K>& 
         */
        Vector<K>& operator/=(const K& alpha);
        
        /**
         * @brief If @a this and @a w have the same dimension, this is equivalent
         * to <em>this->operator[](i) /= w[i]</em> for all @a i . Otherwise,
         * an exception is raised.
         * 
         * @exception An exception is thrown if @a this and @a w have a different
         * number of components.
         * 
         * @param w 
         * @return Vector<K>& 
         */
        Vector<K>& operator/=(const Vector<K>& w);

        // Basic unary operators: +, -
        /**
         * @brief Equivalent to <em>v[i] = +v[i]</em> for all @a i .
         * 
         * @return Vector<K> 
         */
        Vector<K> operator+() const;

        /**
         * @brief Equivalent to <em>v[i] = -v[i]</em> for all @a i .
         * 
         * @return Vector<K> 
         */
        Vector<K> operator-() const;

        // Other operators and functions
        /**
         * @brief Indexing operator that allows to use array syntax for accessing the
         * components of a Vector. Checks for out of index errors. Allows for negative
         * indexes.
         * 
         * @exception Throws an exception if the index supplied is greater than @a N .
         * 
         * @param i 
         * @return K& 
         */
        K& operator[] (const int i);
        
        /**
         * @brief Version of the operator[] for const objects. 
         * 
         * @param i 
         * @return K 
         */
        K operator[] (const int i) const;
        
        /**
         * @brief Slicing operator: equivalent to Python v[first: last]. Returns a subvector
         * of @a this whose first component is @a v[first] and whose last component is
         * <em>v[last - 1]</em>.
         * 
         * @param first The first component of the slice will be @a v[first] .
         * @param last The last component of the slice will be <em>v[last - 1]</em>.
         * 
         * @exception Throws an exception if <em>last < first</em>.
         * 
         * @return Vector<K> 
         */
        Vector<K> slice(int first, int last) const;
        
        /**
         * @brief Returns the Vector reversed: it starts with the last component and it finishes
         * with the first.
         * 
         * @return Vector<K> 
         */
        Vector<K> reverse() const;

        /**
         * @brief Returns the length (or dimension) of the Vector.
         * 
         * @return unsigned int 
         */
        unsigned int size() const;

        /**
         * @brief Norm 2 (euclidean norm) of the Vector.
         * 
         * @return K 
         */
        double norm_2() const;
        
        /**
         * @brief Norm 1 (taxi-cab norm) of the Vector.
         * 
         * @return K 
         */
        double norm_1() const;
        
        /**
         * @brief Norm infinity (maximum norm) of the Vector.
         * 
         * @return K 
         */
        double norm_inf() const;

        /**
         * @brief Norm p of the Vector.
         * 
         * @param p 
         * @return K 
         */
        double norm_p(const double p) const;
        
        /**
         * @brief Returns a string representation of the Vector. It will look similar to
         * (1, 0, -3, 7.3).
         * 
         * @return std::string 
         */
        std::string to_string() const;

        /**
         * @brief Allows to print a Vector @a v to screen with <em>std::cout << v</em>.
         * 
         * @param os 
         * @param v 
         * @return std::ostream& 
         */
        friend std::ostream& operator<< <K>(std::ostream& os, const Vector<K>& v);

        /**
         * @brief Reads the Vector from a binary file. It first reads @a N , then deallocates and
         * reallocates @a componentes and finally reads @a components from the file.
         * 
         * @param file 
         */
        void read_from_file(FILE* file);
        
        /**
         * @brief Writes the Vector to a binary file. It first writes @a N and than
         * @a components .
         * 
         * @param file 
         */
        void write_to_file(FILE* file) const;

        /**
         * @brief Returns the minimum element of the vector.
         * 
         * @param v 
         * @return K 
         */
        friend K min<K>(const Vector<K>& v);

        /**
         * @brief Returns the maximal element of the vector.
         * 
         * @param v 
         * @return K 
         */
        friend K max<K>(const Vector<K>& v);

        /**
         * @brief Returns the sum of all the components of the vector.
         * 
         * @param v 
         * @return K 
         */
        friend K sum<K>(const Vector<K>& v);

        /**
         * @brief Returns the product of all the components of the vector.
         * 
         * @param v 
         * @return K 
         */
        friend K multiply<K>(const Vector<K>& v);

        // All math functions component by component
        /**
         * @brief Equivalent to @a cos(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> cos<K>(const Vector<K>& v);

        /**
         * @brief Equivalent to @a sin(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> sin<K>(const Vector<K>& v);

        /**
         * @brief Equivalent to @a tan(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> tan<K>(const Vector<K>& v);

        /**
         * @brief Equivalent to @a acos(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> acos<K>(const Vector<K>& v);

        /**
         * @brief Equivalent to @a asin(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> asin<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a atan(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> atan<K>(const Vector<K>& v);

        /**
         * @brief Equivalent to <em>atan2(v,w[i])</em> for all @a i .
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> atan2<K>(const K v, const Vector<K>& w);
        
        /**
         * @brief Equivalent to <em>atan2(v[i],w)</em> for all @a i .
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> atan2<K>(const Vector<K>& v, const K w);
        
        /**
         * @brief Equivalent to <em>atan2(v[i],w[i])</em> for all @a i .
         * 
         * @exception Throws an exception if @a v and @a w have different
         * lengths.
         * 
         * @param v 
         * @param w 
         * @return Vector<K> 
         */
        friend Vector<K> atan2<K>(const Vector<K>& v, const Vector<K>& w);
        
        /**
         * @brief Equivalent to @a cosh(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> cosh<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a sinh(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> sinh<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a tanh(v[i]) for all @a i . 
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> tanh<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a acosh(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> acosh<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a asinh(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> asinh<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a atanh(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> atanh<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a exp(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> exp<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to <em>frexp(v[i],exp[i])</em> for all @a i .
         * 
         * @param v 
         * @param exp 
         * @return Vector<K> 
         */
        friend Vector<K> frexp<K>(const Vector<K>& v, Vector<int>* exp);
        
        /**
         * @brief Equivalent to <em>ldexp(v[i],exp)</em> for all @a i .
         * 
         * @param v 
         * @param exp 
         * @return Vector<K> 
         */
        friend Vector<K> ldexp<K>(const Vector<K>& v, const int exp);
        
        /**
         * @brief Equivalent to <em>ldexp(v[i],exp[i])</em> for all @a i .
         * 
         * @param v 
         * @param exp 
         * @return Vector<K> 
         */
        friend Vector<K> ldexp<K>(const Vector<K>& v, const Vector<int>& exp);
        
        /**
         * @brief Equivalent to @a log(v[i]) for all @a i . 
         * 
         * @exception Throws an exception if @a v and @a exp have different
         * lengths.
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> log<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a log10(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> log10<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to <em>modf(v[i],intpart[i])</em> for all @a i .
         * 
         * @param v 
         * @param intpart 
         * @return Vector<K> 
         */
        friend Vector<K> modf<K>(const Vector<K>& v, Vector<K>* intpart);
        
        /**
         * @brief Equivalent to @a exp2(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> exp2<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a expm1(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> expm1<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a ilogb(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> ilogb<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a log1p(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> log1p<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a log2(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> log2<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a logb(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> logb<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to <em>sclbn(v[i],n)</em> for all @a i .
         * 
         * @param v 
         * @param n 
         * @return Vector<K> 
         */
        friend Vector<K> scalbn<K>(const Vector<K>& v, const int n);
        
        /**
         * @brief Equivalent to <em>sclbn(v[i],n[i])</em> for all @a i .
         * 
         * @exception Throws an expception if @a v and @a n have different lenghts.
         * 
         * @param v 
         * @param n 
         * @return Vector<K> 
         */
        friend Vector<K> scalbn<K>(const Vector<K>& v, const Vector<int>& n);
        
        /**
         * @brief Equivalent to <em>sclbln(v[i],n[i])</em> for all @a i .
         * 
         * @param v 
         * @param n 
         * @return Vector<K> 
         */
        friend Vector<K> scalbln<K>(const Vector<K>& v, const long int n);
        
        /**
         * @brief Equivalent to <em>sclbln(v[i],n[i])</em> for all @a i .
         * 
         * @exception Throws an expception if @a v and @a n have different lenghts.
         * 
         * @param v 
         * @param n 
         * @return Vector<K> 
         */
        friend Vector<K> scalbln<K>(const Vector<K>& v, const Vector<long int>& n);
        
        /**
         * @brief Equivalent to <em>pow(v,exponent[i])</em> for all @a i .
         * 
         * @param v 
         * @param exponent 
         * @return Vector<K> 
         */
        friend Vector<K> pow<K>(const K v, const Vector<K>& exponent);
        
        /**
         * @brief Equivalent to <em>pow(v[i],exponent)</em> for all @a i .
         * 
         * @param v 
         * @param exponent 
         * @return Vector<K> 
         */
        friend Vector<K> pow<K>(const Vector<K>& v, const K exponent);
        
        /**
         * @brief Equivalent to <em>sclbln(v[i],exponent[i])</em> for all @a i .
         * 
         * @exception Throws an expception if @a v and @a exponent have different lenghts.
         * 
         * @param v 
         * @param exponent 
         * @return Vector<K> 
         */
        friend Vector<K> pow<K>(const Vector<K>& v, const Vector<K>& exponent);
        
        /**
         * @brief Equivalent to @a sqrt(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> sqrt<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a cbrt(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> cbrt<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to <em>hypot(x,y[i])</em> for all @a i .
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> hypot<K>(const K x, const Vector<K>& y);
        
        /**
         * @brief Equivalent to <em>hypot(x[i],y)</em> for all @a i .
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> hypot<K>(const Vector<K>& x, const K y);
        
        /**
         * @brief Equivalent to <em>hypot(x[i],y[i])</em> for all @a i .
         * 
         * @exception Throws an expception if @a v and @a exponent have different lenghts.
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> hypot<K>(const Vector<K>& x, const Vector<K>& y);
        
        /**
         * @brief Equivalent to @a erf(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> erf<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a erfc(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> erfc<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a tgamma(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> tgamma<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a lgamma(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> lgamma<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a ceil(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> ceil<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a floor(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> floor<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to <em>fmod(numer,denom[i])</em> for all @a i .
         * 
         * @param numer 
         * @param denom 
         * @return Vector<K> 
         */
        friend Vector<K> fmod<K>(const K numer, const Vector<K>& denom);
        
        /**
         * @brief Equivalent to <em>fmod(numer[i],denom)</em> for all @a i .
         * 
         * @param numer 
         * @param denom 
         * @return Vector<K> 
         */
        friend Vector<K> fmod<K>(const Vector<K>& numer, const K denom);
        
        /**
         * @brief Equivalent to <em>fmod(numer[i],denom[i])</em> for all @a i . 
         * 
         * @exception Throws an expception if @a v and @a exponent have different lenghts.
         * 
         * @param numer 
         * @param denom 
         * @return Vector<K> 
         */
        friend Vector<K> fmod<K>(const Vector<K>& numer, const Vector<K>& denom);
        
        /**
         * @brief Equivalent to @a trunc(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> trunc<K>(const Vector<K>& v);
       
        /**
         * @brief Equivalent to @a round(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> round<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a lround(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<long int> 
         */
        friend Vector<long int> lround<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a llround(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<long long int> 
         */
        friend Vector<long long int> llround<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a rint(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> rint<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a lrint(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<long int> 
         */
        friend Vector<long int> lrint<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a llrint(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<long long int> 
         */
        friend Vector<long long int> llrint<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a nearbyint(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> nearbyint<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a remainder(numer,denom[i]) for all @a i .
         * 
         * @param numer 
         * @param denom 
         * @return Vector<K> 
         */
        friend Vector<K> remainder<K>(const K numer, const Vector<K>& denom);
        
        /**
         * @brief Equivalent to @a remainder(numer[i],denom) for all @a i .
         * 
         * @param numer 
         * @param denom 
         * @return Vector<K> 
         */
        friend Vector<K> remainder<K>(const Vector<K>& numer, const K denom);
        
        /**
         * @brief Equivalent to @a remainder(numer[i],denom[i]) for all @a i . 
         * 
         * @exception Throws an exception if @a numer and @a denom have different lengths.
         * 
         * @param numer 
         * @param denom 
         * @return Vector<K> 
         */
        friend Vector<K> remainder<K>(const Vector<K>& numer, const Vector<K>& denom);
        
        /**
         * @brief Equivalent to @a remquo(numer,denom[i],quot[i]) for all @a i .
         * 
         * @param numer 
         * @param denom 
         * @param quot 
         * @return Vector<K> 
         */
        friend Vector<K> remquo<K>(const K numer, const Vector<K>& denom, Vector<int>* quot);
        
        /**
         * @brief Equivalent to @a remquo(numer[i],denom,quot[i]) for all @a i .
         * 
         * @param numer 
         * @param denom 
         * @param quot 
         * @return Vector<K> 
         */
        friend Vector<K> remquo<K>(const Vector<K>& numer, const K denom, Vector<int>* quot);
        
        /**
         * @brief Equivalent to @a remquo(numer[i],denom[i],quot[i]) for all @a i . 
         * 
         * @exception Throws an exception if @a numer and @a denom have different lengths.
         * 
         * @param numer 
         * @param denom 
         * @param quot 
         * @return Vector<K> 
         */
        friend Vector<K> remquo<K>(const Vector<K>& numer, const Vector<K>& denom, Vector<int>* quot);
        
        /**
         * @brief Equivalent to @a copysign(x[i],y) for all @a i .
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> copysign<K>(const Vector<K>& x, const K y);
        
        /**
         * @brief Equivalent to @a copysign(x[i],y[i]) for all @a i . 
         * 
         * @exception Throws an exception if @a x and @a y have different lenghts.
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> copysign<K>(const Vector<K>& x, const Vector<K>& y);
        
        /**
         * @brief Equivalent to @a nan(tagp) for all @a i .
         * 
         * @param N 
         * @param tagp 
         * @return Vector<K> 
         */
        friend Vector<K> nan<K>(const unsigned int N, const char* tagp);
        
        /**
         * @brief Equivalent to @a nextafter(x[i],y[i]) for all @a i . 
         * 
         * @exception Throws an exception if @a x and @a y have different lenghts.
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> nextafter<K>(const Vector<K>& x, const Vector<K>& y);
        
        /**
         * @brief Equivalent to @a fdim(x,y[i]) for all @a i .
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> fdim<K>(const K x, const Vector<K>& y);
        
        /**
         * @brief Equivalent to @a fdim(x[i],y) for all @a i .
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> fdim<K>(const Vector<K>& x, K y);
        
        /**
         * @brief Equivalent to @a fdim(x[i],y[i]) for all @a i . 
         * 
         * @exception Throws an exception if @a x and @a y have different lenghts.
         * 
         * @param x 
         * @param y 
         * @return Vector<K> 
         */
        friend Vector<K> fdim<K>(const Vector<K>& x, const Vector<K>& y);

        /**
         * @brief Equivalent to @a fabs(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<double> fabs<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a abs(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<double> abs<K>(const Vector<K>& v);
        
        /**
         * @brief Equivalent to @a fma(x,y[i],z[i]) for all @a i . 
         * 
         * @exception Throws an exception if the length of @a y does not equal the
         * length of @a z .
         * 
         * @param x 
         * @param y 
         * @param z 
         * @return Vector<K> 
         */
        friend Vector<K> fma<K>(const K x, const Vector<K>& y, const Vector<K>& z);
        
        /**
         * @brief Equivalent to @a fma(x[i],y,z[i]) for all @a i . 
         * 
         * @exception Throws an exception if the length of @a x does not equal the
         * length of @a z .
         * 
         * @param x 
         * @param y 
         * @param z 
         * @return Vector<K> 
         */
        friend Vector<K> fma<K>(const Vector<K>& x, const K y, const Vector<K>& z);
        
        /**
         * @brief Equivalent to @a fma(x[i],y[i],z) for all @a i . 
         * 
         * @exception Throws an exception if the length of @a x does not equal the
         * length of @a y .
         * 
         * @param x 
         * @param y 
         * @param z 
         * @return Vector<K> 
         */
        friend Vector<K> fma<K>(const Vector<K>& x, const Vector<K>& y, const K z);
        
        /**
         * @brief Equivalent to @a fma(x,y,z[i]) for all @a i .
         * 
         * @param x 
         * @param y 
         * @param z 
         * @return Vector<K> 
         */
        friend Vector<K> fma<K>(const K x, const K y, const Vector<K>& z);
        
        /**
         * @brief Equivalent to @a fma(x,y[i],z) for all @a i .
         * 
         * @param x 
         * @param y 
         * @param z 
         * @return Vector<K> 
         */
        friend Vector<K> fma<K>(const K x, const Vector<K>& y, const K z);
        
        /**
         * @brief Equivalent to @a fma(x[i],y,z) for all @a i .
         * 
         * @param x 
         * @param y 
         * @param z 
         * @return Vector<K> 
         */
        friend Vector<K> fma<K>(const Vector<K>& x, const K y, const K z);
        
        /**
         * @brief Equivalent to @a fma(x[i],y[i],z[i]) for all @a i .
         * 
         * @exception Throws and exception if the lenghts of @a x , @a y and @a z are not equal.
         * 
         * @param x 
         * @param y 
         * @param z 
         * @return Vector<K> 
         */
        friend Vector<K> fma<K>(const Vector<K>& x, const Vector<K>& y, const Vector<K>& z);
        
        /**
         * @brief Equivalent to @a fpclassify(v[i]) for all @a i .
         * 
         * @param v 
         * @return Vector<int> 
         */
        friend Vector<int> fpclassify<K>(const Vector<K>& v);
        
        /**
         * @brief Returns @a true if all components of @a v are finite and @a false otherwise.
         * 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool isfinite<K>(const Vector<K>& v);
        
        /**
         * @brief Returns @a true if at least one components of @a v is infinite and @a false otherwise.
         * 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool isinf<K>(const Vector<K>& v);

        /**
         * @brief Returns @a true if at least one components of @a v is nan and @a false otherwise.
         * 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool isnan<K>(const Vector<K>& v);
        
        /**
         * @brief Returns @a true if all components of @a v are normal and @a false otherwise.
         * 
         * @param v 
         * @return true 
         * @return false 
         */
        friend bool isnormal<K>(const Vector<K>& v);

        // Complex functions
        /**
         * @brief Equivalent to conj(v[i]) for all i.
         * 
         * @param v 
         * @return Vector<K> 
         */
        friend Vector<K> conj<K>(const Vector<K>& v);
    };
}

#endif // VECTOR_HPP

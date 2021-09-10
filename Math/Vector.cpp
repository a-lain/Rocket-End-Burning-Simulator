#ifndef VECTOR_CPP
#define VECTOR_CPP
#include "Vector.hpp"
using std::to_string;
using std::abs;
#include <cmath>
#include <exception>
#include <iostream>

using namespace Math;

// The following is needed to make writing and reading from file more universal
#include "../Utilities/FileOperations.hpp"

#include "../Utilities/ToString.hpp"

// Constructors and destructor

template<typename K>
Vector<K>::Vector(const unsigned int N) : components(new K[N]), N(N)
{
    for (unsigned int i = 0; i < N; i++)
    {
        components[i] = 0;
    }
}

template<typename K>
Vector<K>::Vector() : Vector<K>(0)
{

}

template<typename K>
Vector<K>::Vector(const Vector<K>& w) : Vector<K>(w.N)
{
    for (unsigned int i = 0; i < w.N; i++)
    {
        this->components[i] = w[i];
    }
}

template<typename K>
Vector<K>::Vector(Vector<K>&& w)
{
    this->N = w.N;
    this->components = w.components;
    w.components = nullptr;
}

template<typename K>
Vector<K>::Vector(const K* components, const unsigned int N) : Vector<K>(N)
{
    for (unsigned int i = 0; i < N; i++)
    {
        this->components[i] = components[i];
    }
}

template<typename K>
Vector<K>::Vector(std::initializer_list<K> list) : Vector<K>(list.size())
{
    const K* it = begin(list);
    for (unsigned int i = 0; i < N; i++)
    {
        this->components[i] = *it;
        it++;
    }
}

template<typename K>
Vector<K>::~Vector()
{
    delete[] this->components;
    this->components = nullptr;
}


// Assignment operators

template<typename K>
Vector<K>& Vector<K>::operator=(const K& alpha)
{
    for (unsigned int i = 0; i < N; i++)
    {
        components[i] = alpha;
    }
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator=(const Vector<K>& w)
{
    // Self-assignment check
    if (this != &w)
    {
        // If needed, we free the memory assigned to components and
        // than we reserve it again according to the new length.
        if (this->N != w.N)
        {
            delete[] components;
            components = new K[w.N];
            N = w.N;
        }
        for (unsigned int i = 0; i < w.N; i++)
        {
            components[i] = w[i];
        }
    }
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator=(Vector<K>&& w)
{
    this->N = w.N;
    delete[] this->components;
    this->components = w.components;
    w.components = nullptr;
    return *this;
}


// Operators defined componentwise

template<typename K>
Vector<K> Math::operator+(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N == 0 && w.N != 0)
    {
        return w;
    }
    else if (v.N != 0 && w.N == 0)
    {
        return v;
    }
    else if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot sum two vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = v[i] + w[i];
        }
        return res;
    }
}

template<typename K>
Vector<K> Math::operator+(const Vector<K>& v, const K& alpha)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = v[i] + alpha;
    }
    return res;
}

template<typename K>
Vector<K> Math::operator+(const K& alpha, const Vector<K>& v)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = alpha + v[i];
    }
    return res;
}

template<typename K>
Vector<K> Math::operator-(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N == 0 && w.N != 0)
    {
        return -w;
    }
    else if (v.N != 0 && w.N == 0)
    {
        return v;
    }
    else if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot sum substract vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = v[i] - w[i];
        }
        return res;
    }    
}

template<typename K>
Vector<K> Math::operator-(const Vector<K>& v, const K& alpha)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = v[i] - alpha;
    }
    return res;
}

template<typename K>
Vector<K> Math::operator-(const K& alpha, const Vector<K>& v)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = alpha - v[i];
    }
    return res;
}

template<typename K>
Vector<K> Math::operator*(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot multiply two vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = v[i] * w[i];
        }
        return res;
    }
}

template<typename K>
Vector<K> Math::operator*(const K& alpha, const Vector<K>& v)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = alpha * v[i];
    }
    return res;
}

template<typename K>
Vector<K> Math::operator*(const Vector<K>& v, const K& alpha)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = v[i] * alpha;
    }
    return res;
}

template<typename K>
Vector<K> Math::operator/(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot divide two vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = v[i] / w[i];
        }
        return res;
    } 
}

template<typename K>
Vector<K> Math::operator/(const Vector<K>& v, const K& alpha)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = v[i] / alpha;
    }
    return res;
}

template<typename K>
Vector<K> Math::operator/(const K& alpha, const Vector<K>& v)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = alpha / v[i];
    }
    return res;
}

template<typename K>
bool Math::operator<(const Vector<K>& v, const K& alpha)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (v[i] >= alpha)
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator<(const K& alpha, const Vector<K>& v)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (alpha >= v[i])
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator<(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot compare two vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        for (unsigned int i = 0; i < v.N; i++)
        {
            if (v[i] >= w[i])
            {
                return false;
            }
        }
        return true;
    } 
}

template<typename K>
bool Math::operator<=(const Vector<K>& v, const K& alpha)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (v[i] > alpha)
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator<=(const K& alpha, const Vector<K>& v)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (alpha > v[i])
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator<=(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot compare two vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        for (unsigned int i = 0; i < v.N; i++)
        {
            if (v[i] > w[i])
            {
                return false;
            }
        }
        return true;
    } 
}

template<typename K>
bool Math::operator>(const Vector<K>& v, const K& alpha)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (v[i] <= alpha)
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator>(const K& alpha, const Vector<K>& v)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (alpha <= v[i])
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator>(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot compare two vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        for (unsigned int i = 0; i < v.N; i++)
        {
            if (v[i] <= w[i])
            {
                return false;
            }
        }
        return true;
    } 
}

template<typename K>
bool Math::operator>=(const Vector<K>& v, const K& alpha)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (v[i] < alpha)
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator>=(const K& alpha, const Vector<K>& v)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (alpha < v[i])
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator>=(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot compare two vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        for (unsigned int i = 0; i < v.N; i++)
        {
            if (v[i] < w[i])
            {
                return false;
            }
        }
        return true;
    } 
}

template<typename K>
bool Math::operator==(const Vector<K>& v, const K& alpha)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (v[i] != alpha)
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator==(const K& alpha, const Vector<K>& v)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (alpha != v[i])
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator==(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != w.N)
    {
        return false;
    }
    else
    {
        for (unsigned int i = 0; i < v.N; i++)
        {
            if (v[i] != w[i])
            {
                return false;
            }
        }
        return true;
    } 
}

template<typename K>
bool Math::operator!=(const Vector<K>& v, const K& alpha)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (v[i] == alpha)
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator!=(const K& alpha, const Vector<K>& v)
{
    for (unsigned int i = 0; i < v.N; i++)
    {
        if (alpha == v[i])
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::operator!=(const Vector<K>& v, const Vector<K>& w)
{
    return !(v == w);
}


// Concatenation operator

template<typename K>
Vector<K> Math::operator&(const Vector<K>& v, const Vector<K>& w)
{
    Vector<K> res(v.N + w.N);
    // We assign the components of v.
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = v[i];
    }
    // And now the components of w.
    for (unsigned int i = 0; i < w.N; i++)
    {
        res[i + v.N] = w[i];
    }
    return res;
}

template<typename K>
Vector<K> Math::operator&(const Vector<K>& v, const K& w)
{
    return v & Vector<K>({w});
}

template<typename K>
Vector<K> Math::operator&(const K& v, const Vector<K>& w)
{
    return Vector<K>({v}) & w;
}

template<typename K>
K Math::operator|(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot scalarly multiply two vectors that have different number of components! ") + 
            ::to_string(v.N) + " and " + ::to_string(w.N));
    }
    else
    {
        K res = 0;
        for (unsigned int i = 0; i < v.N; i++)
        {
            res += Math::conj(v[i]) * w[i];
        }
        return res;
    }
}

template<typename K>
K Math::vector_product_2d(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != 2 || w.N != 2)
    {
        throw std::runtime_error("Cannot vectorly multiply two vectors that are not 2d vectors.");
    }
    return v[0]*w[1] - v[1]*w[0];
}

template<typename K>
Vector<K> Math::vector_product_3d(const Vector<K>& v, const Vector<K>& w)
{
    if (v.N != 3 || w.N != 3)
    {
        throw std::runtime_error("Cannot vectorly multiply two vectors that are not 3d vectors.");
    }
    return Vector<K>({v[1]*w[2] - v[2]*w[1], v[2]*w[0] - v[0]*w[2], v[0]*w[1] - v[1]*w[0]});
}

// Unary operators
template<typename K>
Vector<K> Vector<K>::operator+()
{
    return *this;
}

template<typename K>
Vector<K> Vector<K>::operator-()
{
    Vector<K> res(*this);
    for (unsigned int i = 0; i < N; i++)
    {
        res[i] = -res[i];
    }
    return res;
}

// Basic math-assignment operators: +=, -=, *=, /=

template<typename K>
Vector<K>& Vector<K>::operator+=(const K& alpha)
{
    for (unsigned int i = 0; i < N; i++)
    {
        components[i] += alpha;
    }
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator+=(const Vector<K>& w)
{
    if (N == 0)
    {
        operator=(w);
    }
    else if (N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot sum vectors that have different number of components! ") + 
            ::to_string(N) + " and " + ::to_string(w.N));
    }
    else
    {
        for (unsigned int i = 0; i < w.N; i++)
        {
            components[i] += w[i];
        }
    }
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator-=(const K& alpha)
{
    for (unsigned int i = 0; i < N; i++)
    {
        components[i] -= alpha;
    }
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator-=(const Vector<K>& w)
{
    if (N == 0)
    {
        operator=(-w);
    }
    else if (N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot substract vectors that have different number of components! ") + 
            ::to_string(N) + " and " + ::to_string(w.N));
    }
    else
    {
        for (unsigned int i = 0; i < w.N; i++)
        {
            components[i] -= w[i];
        }
    }    
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator*=(const K& alpha)
{
    for (unsigned int i = 0; i < N; i++)
    {
        components[i] *= alpha;
    }
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator*=(const Vector<K>& w)
{
    if (N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot multiply vectors that have different number of components! ") + 
            ::to_string(N) + " and " + ::to_string(w.N));
    }
    else
    {
        for (unsigned int i = 0; i < N; i++)
        {
            components[i] *= w[i];
        }
    }    
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator/=(const K& alpha)
{
    for (unsigned int i = 0; i < N; i++)
    {
        components[i] /= alpha;
    }
    return *this;
}

template<typename K>
Vector<K>& Vector<K>::operator/=(const Vector<K>& w)
{
    if (N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot divide vectors that have different number of components! ") + 
            ::to_string(N) + " and " + ::to_string(w.N));
    }
    else
    {
        for (unsigned int i = 0; i < N; i++)
        {
            components[i] /= w[i];
        }
    }
    return *this;
}


// Basic unary operators: +, -

template<typename K>
Vector<K> Vector<K>::operator+() const
{
    return *this;
}

template<typename K>
Vector<K> Vector<K>::operator-() const
{
    Vector<K> res(N);
    for (unsigned int i = 0; i < N; i++)
    {
        res[i] = -components[i];
    }
    return res;
}


// Other operators and functions

template<typename K>
K& Vector<K>::operator[] (const int i)
{
    // We check whether the supplied index is vaild.
    if (i >= (int)N || i < -(int)N)
    {
        throw std::runtime_error("Index out of range! " + ::to_string(i));
    }
    else
    {
        return (i >= 0) ? components[i] : components[N - i];
    }
}

template<typename K>
K Vector<K>::operator[] (const int i) const
{
    // We check whether the supplied index is vaild.
    if (i >= (int)N || i < -(int)N)
    {
        throw std::runtime_error("Index out of range! " + ::to_string(i));
    }
    else
    {
        return (i >= 0) ? components[i] : components[N - i];
    }  
}

template<typename K>
Vector<K> Vector<K>::slice(int first, int last) const
{
    if (first < 0)
    {
        first = N - first;
    }
    if (last < 0)
    {
        last = N - last;
    }
    if (last < first)
    {
        throw std::runtime_error("Last index must be larger than first index!");
    }
    else
    {
        Vector<K> res(last - first);
        for (unsigned int i = first; i < (unsigned int)last; i++)
        {
            res[i - first] = components[i];
        }
        return res;
    }
}

template<typename K>
Vector<K> Vector<K>::reverse() const
{
    // Note that components[N-1] is the last element of the Vector
    // and that we are looping through the Vector backwards.
    Vector<K> res(N);
    for (unsigned int i = 0; i < N; i++)
    {
        res[i] = components[N -1 -i];
    }
    return res;
}

template<typename K>
unsigned int Vector<K>::size() const
{
    return N;
}

template<typename K>
double Vector<K>::norm_2() const
{
    using ::sqrt;
    using ::fabs;
    return sqrt(fabs((*this) | (*this)));
}

template<typename K>
double Vector<K>::norm_1() const
{
    using ::fabs;
    return sum(fabs(*this));
}

template<typename K>
double Vector<K>::norm_inf() const
{
    using ::fabs;
    return max(fabs(*this));
}

template<typename K>
double Vector<K>::norm_p(const double p) const
{
    using ::pow;
    using ::fabs;
    double res = 0;
    for (unsigned int i = 0; i < N; i++)
    {
        res += pow(fabs(components[i]), p);
    }
    return pow(res, 1/p);
}

template<typename K>
std::string Vector<K>::to_string() const
{
    using ::to_string;
    using Math::to_string;
    std::string res = "(";
    for (unsigned int i = 0; i < N - 1; i++)
    {
        res += to_string(components[i]) + ", ";
    }
    res += to_string(components[N - 1]) + ")";
    return res;
}

template<typename K>
std::ostream& Math::operator<< (std::ostream& os, const Vector<K>& v)
{
    return os << v.to_string();
}

template<typename K>
void Vector<K>::read_from_file(FILE* file)
{
    using ::read_from_file;
    read_from_file(N, file);
    delete[] components;
    components = new K[N];
    for (unsigned int i = 0; i < N; i++)
    {
        read_from_file(components[i], file);
    }
}

template<typename K>
void Vector<K>::write_to_file(FILE* file) const
{
    using ::write_to_file;
    write_to_file(N, file);
    for (unsigned int i = 0; i < N; i++)
    {
        write_to_file(components[i], file);
    }
}

// Min and max functions
template<typename K>
K Math::min(const Vector<K>& v)
{
    K res = v[0];
    for (unsigned int i = 1; i < v.N; i++)
    {
        if (v[i] < res)
        {
            res = v[i];
        }
    }
    return res;
}

template<typename K>
K Math::max(const Vector<K>& v)
{
    K res = v[0];
    for (unsigned int i = 1; i < v.N; i++)
    {
        if (v[i] > res)
        {
            res = v[i];
        }
    }
    return res;
}

template<typename K>
K Math::sum(const Vector<K>& v)
{
    K res = 0;
    for (unsigned int i = 0; i < v.N; i++)
    {
        res += v[i];
    }
    return res;
}

template<typename K>
K Math::multiply(const Vector<K>& v)
{
    K res = 1;
    for (unsigned int i = 0; i < v.N; i++)
    {
        res *= v[i];
    }
    return res;
}

// Cmath functions

template<typename K>
Vector<K> Math::cos(const Vector<K>& v)
{
    using ::cos;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = cos(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::sin(const Vector<K>& v)
{
    using ::sin;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = sin(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::tan(const Vector<K>& v)
{
    using ::tan;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = tan(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::asin(const Vector<K>& v)
{
    using ::asin;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = asin(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::acos(const Vector<K>& v)
{
    using ::acos;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = acos(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::atan(const Vector<K>& v)
{
    using ::atan;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = atan(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::atan2(const K v, const Vector<K>& w)
{
    using ::atan2;
    Vector<K> res(w.N);
    for (unsigned int i = 0; i < w.N; i++)
    {
        res[i] = atan2(v, w[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::atan2(const Vector<K>& v, const K w)
{
    using ::atan2;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = atan2(v[i], w);
    }
    return res;
}

template<typename K>
Vector<K> Math::atan2(const Vector<K>& v, const Vector<K>& w)
{
    using ::atan2;
    using ::to_string;
    if (v.N != w.N)
    {
        throw std::runtime_error(std::string(
            "Cannot atan2 two vectors that have different number of components! ") + 
            to_string(v.N) + " and " + to_string(w.N));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = atan2(v[i], w[i]);
        }
        return res;   
    }
}

template<typename K>
Vector<K> Math::cosh(const Vector<K>& v)
{
    using ::cosh;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = cosh(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::sinh(const Vector<K>& v)
{
    using ::sinh;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = sinh(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::tanh(const Vector<K>& v)
{
    using ::tanh;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = tanh(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::acosh(const Vector<K>& v)
{
    using ::acosh;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = acosh(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::asinh(const Vector<K>& v)
{
    using ::asinh;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = asinh(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::atanh(const Vector<K>& v)
{
    using ::atanh;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = atanh(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::exp(const Vector<K>& v)
{
    using ::exp;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = exp(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::frexp(const Vector<K>& v, Vector<int>* exp)
{
    using ::frexp;
    exp = new Vector<int>(v.N);
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = frexp(v[i], &exp->operator[](i));
    }
    return res;
}

template<typename K>
Vector<K> Math::ldexp(const Vector<K>& v, const int exp)
{
    using ::ldexp;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = ldexp(v[i], exp);
    }
    return res;
}

template<typename K>
Vector<K> Math::ldexp(const Vector<K>& v, const Vector<int>& exp)
{
    using ::ldexp;
    using ::to_string;
    if (v.N != exp.size())
    {
        throw std::runtime_error(
            "Cannot ldexp two vectores with different number of components. " +
            to_string(v.N) + " and " + to_string(exp.size()));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = ldexp(v[i], exp[i]);
        }
        return res;
    }
}

template<typename K>
Vector<K> Math::log(const Vector<K>& v)
{
    using ::log;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = log(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::log10(const Vector<K>& v)
{
    using ::log10;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = log10(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::modf(const Vector<K>& v, Vector<K>* intpart)
{
    using ::modf;
    intpart = new Vector<K>(v.N);
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = modf(v[i], &intpart->operator[](i));
    }
    return res;
}

template<typename K>
Vector<K> Math::exp2(const Vector<K>& v)
{
    using ::exp2;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = exp2(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::expm1(const Vector<K>& v)
{
    using ::expm1;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = expm1(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::ilogb(const Vector<K>& v)
{
    using ::ilogb;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = ilogb(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::log1p(const Vector<K>& v)
{
    using ::log1p;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = log1p(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::log2(const Vector<K>& v)
{
    using ::log2;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = log2(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::logb(const Vector<K>& v)
{
    using ::logb;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = logb(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::scalbn(const Vector<K>& v, const int n)
{
    using ::scalbn;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = scalbn(v[i], n);
    }
    return res;
}

template<typename K>
Vector<K> Math::scalbn(const Vector<K>& v, const Vector<int>& n)
{
    using ::scalbn;
    using ::to_string;
    if (v.N != n.size())
    {
        throw std::runtime_error(
            "Cannot sclbn two vectores with different number of components. " +
            to_string(v.N) + " and " + to_string(n.size()));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = scalbn(v[i], n[i]);
        }
        return res;
    }
}

template<typename K>
Vector<K> Math::scalbln(const Vector<K>& v, const long int n)
{
    using ::scalbln;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = scalbln(v[i], n);
    }
    return res;
}

template<typename K>
Vector<K> scalbln(const Vector<K>& v, const Vector<long int>& n)
{
    using ::scalbln;
    using ::to_string;
    if (v.N != n.size())
    {
        throw std::runtime_error(
            "Cannot sclbln two vectores with different number of components. " +
            to_string(v.N) + " and " + to_string(n.size()));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = scalbln(v[i], n[i]);
        }
        return res;
    }
}

template<typename K>
Vector<K> Math::pow(const K v, const Vector<K>& exponent)
{
    using ::pow;
    Vector<K> res(exponent.N);
    for (unsigned int i = 0; i < exponent.N; i++)
    {
        res[i] = pow(v, exponent[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::pow(const Vector<K>& v, const K exponent)
{
    using ::pow;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = pow(v[i], exponent);
    }
    return res;
}

template<typename K>
Vector<K> Math::pow(const Vector<K>& v, const Vector<K>& exponent)
{
    using ::pow;
    using ::to_string;
    if (v.N != exponent.N)
    {
        throw std::runtime_error(
            "Cannot pow two vectores with different number of components. " +
            to_string(v.N) + " and " + to_string(exponent.N));
    }
    else
    {
        Vector<K> res(v.N);
        for (unsigned int i = 0; i < v.N; i++)
        {
            res[i] = pow(v[i], exponent[i]);
        }
        return res;
    }
}

template<typename K>
Vector<K> Math::sqrt(const Vector<K>& v)
{
    using ::sqrt;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = sqrt(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::cbrt(const Vector<K>& v)
{
    using ::cbrt;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = cbrt(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::hypot(const K x, const Vector<K>& y)
{
    using ::hypot;
    Vector<K> res(y.N);
    for (unsigned int i = 0; i < y.N; i++)
    {
        res[i] = hypot(x, y[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::hypot(const Vector<K>& x, const K y)
{
    using ::hypot;
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = hypot(x[i], y);
    }
    return res;
}

template<typename K>
Vector<K> Math::hypot(const Vector<K>& x, const Vector<K>& y)
{
    using ::hypot;
    using ::to_string;
    if (x.N != y.N)
    {
        throw std::runtime_error(
            "Cannot hypot two vectores with different number of components. " +
            to_string(x.N) + " and " + to_string(y.N));
    }
    else
    {
        Vector<K> res(x.N);
        for (unsigned int i = 0; i < x.N; i++)
        {
            res[i] = hypot(x[i], y[i]);
        }
        return res;
    }
}

template<typename K>
Vector<K> Math::erf(const Vector<K>& v)
{
    using ::erf;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = erf(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::erfc(const Vector<K>& v)
{
    using ::erfc;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = erfc(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::tgamma(const Vector<K>& v)
{
    using ::tgamma;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = tgamma(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::lgamma(const Vector<K>& v)
{
    using ::lgamma;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = lgamma(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::ceil(const Vector<K>& v)
{
    using ::ceil;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = ceil(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::floor(const Vector<K>& v)
{
    using ::floor;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = floor(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::fmod(const K numer, const Vector<K>& denom)
{
    using ::fmod;
    Vector<K> res(denom.N);
    for (unsigned int i = 0; i < denom.N; i++)
    {
        res[i] = fmod(numer, denom[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::fmod(const Vector<K>& numer, const K denom)
{
    using ::fmod;
    Vector<K> res(numer.N);
    for (unsigned int i = 0; i < numer.N; i++)
    {
        res[i] = fmod(numer[i], denom);
    }
    return res;
}

template<typename K>
Vector<K> Math::fmod(const Vector<K>& numer, const Vector<K>& denom)
{
    using ::fmod;
    if (numer.N != denom.N)
    {
        throw std::runtime_error(
            "Cannot fmod two vectores with different number of components. " +
            ::to_string(numer.N) + " and " + ::to_string(denom.N));
    }
    else
    {
        Vector<K> res(numer.N);
        for (unsigned int i = 0; i < numer.N; i++)
        {
            res[i] = fmod(numer[i], denom[i]);
        }
        return res;
    }
}

template<typename K>
Vector<K> Math::trunc(const Vector<K>& v)
{
    using ::trunc;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = trunc(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::round(const Vector<K>& v)
{
    using ::round;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = round(v[i]);
    }
    return res;
}

template<typename K>
Vector<long int> Math::lround(const Vector<K>& v)
{
    using ::lround;
    Vector<long int> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = lround(v[i]);
    }
    return res;
}

template<typename K>
Vector<long long int> Math::llround(const Vector<K>& v)
{
    using ::llround;
    Vector<long long int> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = llround(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::rint(const Vector<K>& v)
{
    using ::rint;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = rint(v[i]);
    }
    return res;
}

template<typename K>
Vector<long int> Math::lrint(const Vector<K>& v)
{
    using ::lrint;
    Vector<long int> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = lrint(v[i]);
    }
    return res;
}

template<typename K>
Vector<long long int> Math::llrint(const Vector<K>& v)
{
    using ::llrint;
    Vector<long long int> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = llrint(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::nearbyint(const Vector<K>& v)
{
    using ::nearbyint;
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = nearbyint(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::remainder(const K numer, const Vector<K>& denom)
{
    using ::remainder;
    Vector<K> res(denom.N);
    for (unsigned int i = 0; i < denom.N; i++)
    {
        res[i] = remainder(numer, denom[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::remainder(const Vector<K>& numer, const K denom)
{
    using ::remainder;
    Vector<K> res(numer.N);
    for (unsigned int i = 0; i < numer.N; i++)
    {
        res[i] = remainder(numer[i], denom);
    }
    return res;
}

template<typename K>
Vector<K> Math::remainder(const Vector<K>& numer, const Vector<K>& denom)
{
    using ::remainder;
    using ::to_string;
    if (numer.N != denom.N)
    {
        throw std::runtime_error("Cannot remainder two vectors of different lenghs: "
            + to_string(numer.N) + " and " + to_string(denom.N));
    }
    Vector<K> res(numer.N);
    for (unsigned int i = 0; i < numer.N; i++)
    {
        res[i] = remainder(numer[i], denom[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::remquo(const K numer, const Vector<K>& denom, Vector<int>* quot)
{
    using ::remquo;
    Vector<K> res(denom.N);
    quot = new Vector<int>(denom.N);
    for (unsigned int i = 0; i < denom.N; i++)
    {
        res[i] = remquo(numer, denom[i], &quot->operator[](i));
    }
    return res;
}

template<typename K>
Vector<K> Math::remquo(const Vector<K>& numer, const K denom, Vector<int>* quot)
{
    using ::remquo;
    Vector<K> res(numer.N);
    quot = new Vector<int>(numer.N);
    for (unsigned int i = 0; i < numer.N; i++)
    {
        res[i] = remquo(numer[i], denom, &quot->operator[](i));
    }
    return res;
}

template<typename K>
Vector<K> Math::remquo(const Vector<K>& numer, const Vector<K>& denom, Vector<int>* quot)
{
    using ::remquo;
    using ::to_string;
    if (numer.N != denom.N)
    {
        throw std::runtime_error("Cannot remquo two vectors of different lenghs: "
            + to_string(numer.N) + " and " + to_string(denom.N));
    }
    Vector<K> res(numer.N);
    quot = new Vector<int>(numer.N);
    for (unsigned int i = 0; i < numer.N; i++)
    {
        res[i] = remquo(numer[i], denom[i], &quot->operator[](i));
    }
    return res;
}

template<typename K>
Vector<K> Math::copysign(const Vector<K>& x, const K y)
{
    using ::copysign;
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = copysign(x[i], y);
    }
    return res;
}

template<typename K>
Vector<K> Math::copysign(const Vector<K>& x, const Vector<K>& y)
{
    using ::copysign;
    using ::to_string;
    if (x.N != y.N)
    {
        throw std::runtime_error("Cannot copysign two vectors of different lenghs: "
            + to_string(x.N) + " and " + to_string(y.N));
    }
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = copysign(x[i], y[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::nan(const unsigned int N, const char* tagp)
{
    using ::nan;
    Vector<K> res(N);
    for (unsigned int i = 0; i < N; i++)
    {
        res[i] = nan(tagp);
    }
    return res;
}

template<typename K>
Vector<K> Math::nextafter(const Vector<K>& x, const Vector<K>& y)
{
    using ::nextafter;
    using ::to_string;
    if (x.N != y.N)
    {
        throw std::runtime_error("Cannot nextafter two vectors of different lenghs: "
            + to_string(x.N) + " and " + to_string(y.N));
    }
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = nextafter(x[i], y[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::fdim(const K x, const Vector<K>& y)
{
    using ::fdim;
    Vector<K> res(y.N);
    for (unsigned int i = 0; i < y.N; i++)
    {
        res[i] = fdim(x, y[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::fdim(const Vector<K>& x, const K y)
{
    using ::fdim;
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = fdim(x[i], y);
    }
    return res;
}

template<typename K>
Vector<K> Math::fdim(const Vector<K>& x, const Vector<K>& y)
{
    using ::fdim;
    using ::to_string;
    if (x.N != y.N)
    {
        throw std::runtime_error("Cannot fdim two vectors of different lenghs: "
            + to_string(x.N) + " and " + to_string(y.N));
    }
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = fdim(x[i], y[i]);
    }
    return res;
}

template<typename K>
Vector<double> Math::fabs(const Vector<K>& v)
{
    using ::fabs;
    Vector<double> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = fabs(v[i]);
    }
    return res;
}

template<typename K>
Vector<double> Math::abs(const Vector<K>& v)
{
    using ::abs;
    Vector<double> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = abs(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::fma(const K x, const Vector<K>& y, const Vector<K>& z)
{
    using ::fma;
    using ::to_string;
    if (y.N != z.N)
    {
        throw std::runtime_error("Cannot fma two vectors of different lenghs: "
            + to_string(y.N) + " and " + to_string(z.N));
    }
    Vector<K> res(y.N);
    for (unsigned int i = 0; i < y.N; i++)
    {
        res[i] = fma(x, y[i], z[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::fma(const Vector<K>& x, const K y, const Vector<K>& z)
{
    using ::fma;
    using ::to_string;
    if (x.N != z.N)
    {
        throw std::runtime_error("Cannot fma two vectors of different lenghs: "
            + to_string(x.N) + " and " + to_string(z.N));
    }
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = fma(x[i], y, z[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::fma(const Vector<K>& x, const Vector<K>& y, const K z)
{
    using ::fma;
    using ::to_string;
    if (x.N != y.N)
    {
        throw std::runtime_error("Cannot fma two vectors of different lenghs: "
            + to_string(x.N) + " and " + to_string(y.N));
    }
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = fma(x[i], y[i], z);
    }
    return res;
}

template<typename K>
Vector<K> Math::fma(const K x, const K y, const Vector<K>& z)
{
    using ::fma;
    Vector<K> res(z.N);
    for (unsigned int i = 0; i < z.N; i++)
    {
        res[i] = fma(x, y, z[i]);
    }
    return res;
}

template<typename K>
Vector<K> Math::fma(const K x, const Vector<K>& y, const K z)
{
    using ::fma;
    Vector<K> res(y.N);
    for (unsigned int i = 0; i < y.N; i++)
    {
        res[i] = fma(x, y[i], z);
    }
    return res;
}

template<typename K>
Vector<K> Math::fma(const Vector<K>& x, const K y, const K z)
{
    using ::fma;
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = fma(x[i], y, z);
    }
    return res;
}

template<typename K>
Vector<K> Math::fma(const Vector<K>& x, const Vector<K>& y, const Vector<K>& z)
{
    using ::fma;
    using ::to_string;
    if (x.N != y.N || y.N != z.N)
    {
        throw std::runtime_error("Cannot fma three vectors of different lenghs: "
            + to_string(x.N) + ", " + to_string(y.N) + " and " + to_string(z.N));
    }
    Vector<K> res(x.N);
    for (unsigned int i = 0; i < x.N; i++)
    {
        res[i] = fma(x[i], y[i], z[i]);
    }
    return res;
}

template<typename K>
Vector<int> Math::fpclassify(const Vector<K>& v)
{
    using std::fpclassify;
    Vector<int> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = fpclassify(v[i]);
    }
    return res;
}

template<typename K>
bool Math::isfinite(const Vector<K>& v)
{
    using std::isfinite;
    for (unsigned int i = 0; i < v.N; i++)
    {
        if(!isfinite(v[i]))
        {
            return false;
        }
    }
    return true;
}

template<typename K>
bool Math::isinf(const Vector<K>& v)
{
    using std::isinf;
    for (unsigned int i = 0; i < v.N; i++)
    {
        if(isinf(v[i]))
        {
            return true;
        }
    }
    return false;
}

template<typename K>
bool Math::isnan(const Vector<K>& v)
{
    using std::isnan;
    for (unsigned int i = 0; i < v.N; i++)
    {
        if(isnan(v[i]))
        {
            return true;
        }
    }
    return false;
}

template<typename K>
bool Math::isnormal(const Vector<K>& v)
{
    using std::isnormal;
    for (unsigned int i = 0; i < v.N; i++)
    {
        if(!isnormal(v[i]))
        {
            return false;
        }
    }
    return true;
}


// Complex functions

template<typename K>
Vector<K> real(const Vector<std::complex<K>>& v)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = std::real(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> imag(const Vector<std::complex<K>>& v)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = std::imag(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> arg(const Vector<std::complex<K>>& v)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = std::arg(v[i]);
    }
    return res;
}

template<typename K>
Vector<K> conj(const Vector<K>& v)
{
    Vector<K> res(v.N);
    for (unsigned int i = 0; i < v.N; i++)
    {
        res[i] = conj(v[i]);
    }
    return res;
}

#endif // VECTOR_CPP

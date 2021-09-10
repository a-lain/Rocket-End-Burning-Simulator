#ifndef MATRIX_CPP
#define MATRIX_CPP

#include "Matrix.hpp"
#include <limits>
#include <cmath>
#include <iostream>

using namespace Math;

template <typename K>
Matrix<K>::Matrix(unsigned int m, unsigned int n, std::vector<K> elements) : Matrix(m, n)
{
    if (elements.size() != m*n)
    {
        throw std::runtime_error("Unable to create the matrix.");
    }
    else
    {
        for (unsigned int i = 0; i < m; i++)
        {
            for (unsigned int j = 0; j < n; j++)
            {
                this->operator()(i,j) = elements[i*n + j];
            }
        }
    }
}

template <typename K>
Matrix<K>::Matrix(unsigned int m, unsigned int n)
{
    _m = m;
    vec = new Vector<K>[m];
    for (unsigned int i = 0; i < m; i++)
    {
        vec[i] = Vector<K>(n);
    }
}

template <typename K>
Matrix<K>::Matrix(const Matrix<K>& B) : Matrix(B.m(), B.n())
{
    for (unsigned int i = 0; i < m(); i++)
    {
        vec[i] = B.vec[i];
    }
}

template <typename K>
Matrix<K>& Matrix<K>::operator=(const Matrix<K>& B)
{
    // If needed, we free the memory assigned to components and
    // than we reserve it again according to the new length.
    if (m() != B.m() || n() != B.n())
    {
        delete[] vec;
        vec = new Vector<K>[B.m()];
        _m = B.m();
    }

    for (unsigned int i = 0; i < m(); i++)
    {
        vec[i] = B.vec[i];
    }
    return *this;
}

template <typename K>
Matrix<K>::~Matrix()
{
    delete[] vec;
}

template <typename K>
Matrix<K>::operator Matrix<std::complex<K>>() const
{
    Matrix<std::complex<K>> res(m(), n());
    for (unsigned int i = 0; i < m(); i++)
    {
        for (unsigned int j = 0; j < n(); j++)
        {
            res(i,j) = std::complex<K>((*this)(i,j));
        }   
    }
    return res;
}

template <typename K>
Matrix<K> Math::operator+(const Matrix<K>& A, const Matrix<K>& B)
{
    if (A.m() != B.m() || A.n() != B.n())
    {
        throw std::runtime_error("Cannot multiply matrices of different sizes!");
    }
    else
    {
        Matrix<K> res(A.m(), A.n());
        for (unsigned int i = 0; i < A.m(); i++)
        {
            res.vec[i] = A.vec[i] + B.vec[i];
        }
        return res;
    }
}

template <typename K>
Matrix<K> Math::operator+(const Matrix<K>& A, const K& alpha)
{
    Matrix<K> res = A;
    for (unsigned int i = 0; i < A.m(); i++)
    {
        res.vec[i] += alpha;
    }
    return res;
}

template <typename K>
Matrix<K> Math::operator+(const K& alpha, const Matrix<K>& A)
{
    Matrix<K> res = A;
    for (unsigned int i = 0; i < A.m(); i++)
    {
        res.vec[i] += alpha;
    }
    return res;
}

template <typename K>
Matrix<K> Math::operator-(const Matrix<K>& A, const Matrix<K>& B)
{
    if (A.m() != B.m() || A.n() != B.n())
    {
        throw std::runtime_error("Cannot multiply matrices of different sizes!");
    }
    else
    {
        Matrix<K> res(A.m(), A.n());
        for (unsigned int i = 0; i < A.m(); i++)
        {
            res.vec[i] = A.vec[i] - B.vec[i];
        }
        return res;
    }
}

template <typename K>
Matrix<K> Math::operator-(const Matrix<K>& A, const K& alpha)
{
    Matrix<K> res = A;
    for (unsigned int i = 0; i < A.m(); i++)
    {
        res.vec[i] -= alpha;
    }
    return res;
}

template <typename K>
Matrix<K> Math::operator-(const K& alpha, const Matrix<K>& A)
{
    Matrix<K> res = A;
    for (unsigned int i = 0; i < A.m(); i++)
    {
        res.vec[i] = alpha - res.vec[i];
    }
    return res;
}

template <typename K>
Matrix<K> Math::operator*(const K& alpha, const Matrix<K>& A)
{
    Matrix<K> res = A;
    for (unsigned int i = 0; i < A.m(); i++)
    {
        res.vec[i] *= alpha;
    }
    return res;
}

template <typename K>
Matrix<K> Math::operator*(const Matrix<K>& A, const K& alpha)
{
    Matrix<K> res = A;
    for (unsigned int i = 0; i < A.m(); i++)
    {
        res.vec[i] *= alpha;
    }
    return res;
}

template <typename K>
Matrix<K> Math::operator/(const Matrix<K>& A, const K& alpha)
{
    Matrix<K> res = A;
    for (unsigned int i = 0; i < A.m(); i++)
    {
        res.vec[i] /= alpha;
    }
    return res;
}

template <typename K>
Matrix<K> Math::operator|(const Matrix<K>& A, const Matrix<K>& B)
{
    Matrix<K> res = Matrix<K>::zero(A.m(), B.n());
    if (A.n() != B.m())
    {
        throw std::runtime_error("Cannot multiply matrices!");
    }
    else
    {
        for (unsigned int i = 0; i < A.m(); i++)
        {
            for (unsigned int j = 0; j < B.n(); j++)
            {
                res(i,j) = A(i,"")|B("",j);
            }
        }
    }
    return res;
}

template <typename K>
Vector<K> Math::operator|(const Matrix<K>& A, const Vector<K>& v)
{
    Vector<K> res = Vector<K>(A.m());
    if (A.n() != v.size())
    {
        throw std::runtime_error("Cannot multiply the matrix and the vector!");
    }
    else
    {
        for (unsigned int i = 0; i < A.m(); i++)
        {
            res[i] = A(i,"")|v;
        }
    }
    return res;
}

template <typename K>
Vector<K> Math::operator|(const Vector<K>& v, const Matrix<K>& A)
{
    Vector<K> res = Vector<K>(A.n());
    if (A.m() != v.size())
    {
        throw std::runtime_error("Cannot multiply the matrix and the vector!");
    }
    else
    {
        for (unsigned int j = 0; j < A.n(); j++)
        {
            res[j] = v|A("",j);
        }
    }
    return res;
}

template <typename K>
Matrix<K>& Matrix<K>::operator+=(const K& alpha)
{
    for (unsigned int i = 0; i < m(); i++)
    {
        vec[i] += alpha;
    }
    return *this;
}

template <typename K>
Matrix<K>& Matrix<K>::operator-=(const K& alpha)
{
    for (unsigned int i = 0; i < m(); i++)
    {
        vec[i] -= alpha;
    }
    return *this;
}

template <typename K>
Matrix<K>& Matrix<K>::operator*=(const K& alpha)
{
    for (unsigned int i = 0; i < m(); i++)
    {
        vec[i] *= alpha;
    }
    return *this;
}

template <typename K>
Matrix<K>& Matrix<K>::operator/=(const K& alpha)
{
    for (unsigned int i = 0; i < m(); i++)
    {
        vec[i] /= alpha;
    }
    return *this;
}

template <typename K>
Vector<K>& Matrix<K>::operator()(const unsigned int i, const std::string& j)
{
    return vec[i];
}

template <typename K>
Vector<K> Matrix<K>::operator()(const unsigned int i, const std::string& j) const
{
    return vec[i];
}

template <typename K>
Vector<K> Matrix<K>::operator()(const std::string& str, const unsigned int j) const
{
    Vector<K> res(m());
    for (unsigned int i = 0; i < m(); i++)
    {
        res[i] = this->operator()(i,j);
    }
    return res;
}

template <typename K>
K& Matrix<K>::operator()(const unsigned int i, const unsigned int j)
{
    return vec[i][j];
}

template <typename K>
K Matrix<K>::operator()(const unsigned int i, const unsigned int j) const
{
    return vec[i][j];
}

template <typename K>
Matrix<K> Matrix<K>::zero(const unsigned int m, const unsigned int n)
{
    return Matrix<K>(m, n);
}

template <typename K>
Matrix<K> Matrix<K>::zero(const unsigned int n)
{
    return Matrix<K>::zero(n, n);
}

template <typename K>
Matrix<K> Matrix<K>::identity(const unsigned int n)
{
    Matrix<K> A = Matrix<K>::zero(n);
    for (unsigned int i = 0; i < n; i++)
    {
        A(i,i) = 1;
    }
    return A;
}

template <typename K>
unsigned int Matrix<K>::m() const
{
    return _m;
}

template <typename K>
unsigned int Matrix<K>::n() const
{
    return vec[0].size();
}

template <typename K>
bool Matrix<K>::is_square() const
{
    return m() == n();
}

template <typename K>
void Matrix<K>::LU(Matrix<K>& L, Matrix<K>& U, std::vector<unsigned int>& P) const
{
    if (!is_square())
    {
        throw std::runtime_error("Cannot calculate LU factorization of non-square matrix.\n");
    }
    else
    {
        L = Matrix<K>(n(), n());
        U = *this;
        P = std::vector<unsigned int>(n());
        for (unsigned int i = 0; i < n(); i++)
        {
            P[i] = i;
        }

        for (unsigned k = 0; k < n() - 1; k++)
        {
            // We search for the pivot.
            unsigned int i_pivot = k;
            double max_value = std::numeric_limits<double>::min();
            for (unsigned int i = k; i < n(); i++)
            {
                if (::fabs(U(i,k)) > max_value)
                {
                    max_value = ::fabs(U(i,k));
                    i_pivot = i;
                }
            }

            // We exchange the rows.
            if (i_pivot != k)
            {
                Vector<K> aux_row_U = U(k, "");
                Vector<K> aux_row_L = L(k, "");
                unsigned int aux_index = P[k];
                U(k, "") = U(i_pivot, "");
                L(k, "") = L(i_pivot, "");
                P[k] = P[i_pivot];
                U(i_pivot, "") = aux_row_U;
                L(i_pivot, "") = aux_row_L;
                P[i_pivot] = aux_index;
            }

            // We apply Gaussian Elimination.
            L(k,k) = 1;
            for (unsigned int i = k+1; i < n(); i++)
            {
                L(i,k) = U(i,k) / U(k,k);
                U(i,k) = 0;
                for (unsigned int j = k+1; j < n(); j++)
                {
                    U(i,j) = U(i,j) - L(i,k)*U(k,j);
                }
            }
        }
        L(n() - 1, n() - 1) = 1;
    }
}

template <typename K>
K Matrix<K>::det() const
{
    if (!is_square())
    {
        throw std::runtime_error("Cannot calculate determinant of non-square matrix.");
    }
    else
    {
        const Matrix<K>& a = *this;
        unsigned int dim = n();
        if (dim == 1)
        {
            return a(0,0);
        }
        else if (dim == 2)
        {
            return a(0,0)*a(1,1) - a(0,1)*a(1,0);
        }
        else if (dim == 3) // Sarrus rule.
        {
            return a(0,0)*a(1,1)*a(2,2) + a(0,1)*a(1,2)*a(2,0) + a(0,2)*a(1,0)*a(2,1)
                - a(2,0)*a(1,1)*a(0,2) - a(2,1)*a(1,2)*a(0,0) - a(2,2)*a(1,0)*a(0,1);
        }
        else // LU factorization.
        {
            // TODO
        }
    }
}

template <typename K>
K Matrix<K>::tr() const
{
    if (!is_square())
    {
        throw std::runtime_error("Cannot calculate trace of a non-square matrix.");
    }
    else
    {
        K res = 0;
        for (unsigned int i = 0; i < n(); i++)
        {
            res+= (*this)(i,i);
        }
        return res;
    }
}

template <typename K>
std::vector<std::complex<double>> Matrix<K>::eigenvalues() const
{
    if (!is_square())
    {
        throw std::runtime_error("Cannot calculate eigenvalues of a non-square matrix.");
    }
    else
    {
        Matrix<K> a = *this;
        unsigned int dim = n();
        if (dim == 1)
        {
            return std::vector<std::complex<double>>({std::complex<double>(0,0)});
        }
        else if (dim == 2)
        {
            K trA = tr();
            K detA = det();
            std::complex<double> s = std::sqrt(std::complex<double>(trA*trA - 4.*detA))/2.;
            std::complex<double> p = trA/2.;
            return std::vector<std::complex<double>>({p - s, p + s});
        }
        else if (dim == 3)
        {
            K q = tr()/3.;
            Matrix<K> B = *this - q*Matrix<K>::identity(3);
            K p = std::sqrt((B|B).tr() / 6.);
            B /= p;
            std::complex<double> b = std::acos(std::complex<double>(B.det()/2.)) / 3.;
            std::complex<double> beta1 = 2.*std::cos(b);
            std::complex<double> beta2 = 2.*std::cos(b + 2*M_PI/3);
            std::complex<double> beta3 = 2.*std::cos(b + 4*M_PI/3);
            return std::vector<std::complex<double>>({
                p*beta1 + q, p*beta2 + q, p*beta3 + q});
        }
        else
        {
            // TODO
        }
    }
}

template <typename K>
std::vector<Vector<std::complex<double>>> Matrix<K>::eigenvectors(
        std::vector<std::complex<double>>& eigenvalues, 
        const bool calculate_eigenvalues) const
{
    std::vector<Vector<std::complex<double>>> res;
    if (!is_square())
    {
        throw std::runtime_error("Cannot calculate eigenvectors of a non-square matrix.");
    }
    else
    {
        if (n() <= 3) // We use the Carley-Hamilton theorem.
        {
            // If we haven't computed the eigenvalues, we calculate them.
            if (calculate_eigenvalues) 
            {
                eigenvalues = this->eigenvalues();
            }

            for (unsigned int i = 0; i < n(); i++)
            {
                Matrix<std::complex<double>> B = Matrix<std::complex<double>>::identity(n());
                for (unsigned j = 0; j < n(); j++)
                {
                    if (j == i)
                    {
                        continue;
                    }
                    else
                    {
                        B = B | (Matrix<std::complex<double>>(*this)
                            - eigenvalues[j]*Matrix<std::complex<double>>::identity(n()));
                    }
                }
                // We search for non-zero columns of B.
                double max = 0;
                unsigned int index = 0;
                for (unsigned j = 0; j < n(); j++)
                {
                    double temp = B("", j).norm_2();
                    if (temp > max)
                    {
                        max = temp;
                        index = j;
                    }
                }
                res.push_back(B("", index)/std::complex<double>(B("", index).norm_2()));
            }
            return res;
        }
        else // We use other methods to recompute eigenvalues.
        {

        }
    }
}

template <typename K>
std::string Matrix<K>::to_string() const
{
    std::string res = "M(";
    for (unsigned int i = 0; i < m() - 1; i++)
    {
        res += vec[i].to_string() + ",\n";
    }
    res += vec[m() - 1].to_string() + ")";
    return res;
}

#endif // MATRIX_CPP
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Vector.hpp"
#include <vector>

double inline fabs(std::complex<double> z)
{
    return std::abs(z);
}

namespace Math
{
    // Forward declaration to make function declaration possible.
    template<typename K>
    class Matrix;

    // We declare all friend operators.
    template <typename K>
    Matrix<K> operator+(const Matrix<K>& A, const Matrix<K>& B);
    template <typename K>
    Matrix<K> operator+(const Matrix<K>& A, const K& alpha);
    template <typename K>
    Matrix<K> operator+( const K& alpha, const Matrix<K>& A);
    template <typename K>
    Matrix<K> operator-(const Matrix<K>& A, const Matrix<K>& B);
    template <typename K>
    Matrix<K> operator-(const Matrix<K>& A, const K& alpha);
    template <typename K>
    Matrix<K> operator-( const K& alpha, const Matrix<K>& A);
    template <typename K>
    Matrix<K> operator*(const Matrix<K>& A, const K& alpha);
    template <typename K>
    Matrix<K> operator*( const K& alpha, const Matrix<K>& A);
    template <typename K>
    Matrix<K> operator/(const Matrix<K>& A, const K& alpha);
    template <typename K>
    Matrix<K> operator|(const Matrix<K>& A, const Matrix<K>& B);
    template <typename K>
    Vector<K> operator|(const Matrix<K>& A, const Vector<K>& v);
    template <typename K>
    Vector<K> operator|(const Vector<K>& v, const Matrix<K>& A);

    // In the future, define exp, log, sqrt, ...

    template <typename K>
    class Matrix
    {
        public:
        Matrix(unsigned int m = 0, unsigned int n = 0);
        Matrix(unsigned int m, unsigned int n, std::vector<K> elements);

        Matrix(const Matrix<K>& B);
        Matrix& operator=(const Matrix<K>& B);
        ~Matrix();

        operator Matrix<std::complex<K>>() const;

        friend Matrix<K> operator+<K>(const Matrix<K>& A, const Matrix<K>& B);
        friend Matrix<K> operator+<K>(const Matrix<K>& A, const K& alpha);
        friend Matrix<K> operator+<K>(const K& alpha, const Matrix<K>& A);
        friend Matrix<K> operator-<K>(const Matrix<K>& A, const Matrix<K>& B);
        friend Matrix<K> operator-<K>(const Matrix<K>& A, const K& alpha);
        friend Matrix<K> operator-<K>(const K& alpha, const Matrix<K>& A);
        friend Matrix<K> operator*<K>(const Matrix<K>& A, const K& alpha);
        friend Matrix<K> operator*<K>(const K& alpha, const Matrix<K>& A);
        friend Matrix<K> operator/<K>(const Matrix<K>& A, const K& alpha);
        friend Matrix<K> operator|<K>(const Matrix<K>& A, const Matrix<K>& B);
        friend Vector<K> operator|<K>(const Matrix<K>& A, const Vector<K>& v);
        friend Vector<K> operator|<K>(const Vector<K>& v, const Matrix<K>& A);

        Matrix<K>& operator+=(const K& alpha);
        Matrix<K>& operator-=(const K& alpha);
        Matrix<K>& operator*=(const K& alpha);
        Matrix<K>& operator/=(const K& alpha);

        /*! Returns row i.*/
        Vector<K>& operator()(const unsigned int i, const std::string& j);

        /*! Returns row i.*/
        Vector<K> operator()(const unsigned int i, const std::string& j) const;

        /*! Returns column j.*/
        Vector<K> operator()(const std::string& i, const unsigned int j) const;

        K& operator()(const unsigned int i, const unsigned int j);
        K operator()(const unsigned int i, const unsigned int j) const;

        static Matrix<K> zero(const unsigned int m, const unsigned int n);
        static Matrix<K> zero(const unsigned int n);
        static Matrix<K> identity(const unsigned int n);

        /**
         * @brief Horizontal dimension.
         * 
         * @return unsigned int 
         */
        unsigned int m() const;

        /**
         * @brief Vertical dimension.
         * 
         * @return unsigned int 
         */
        unsigned int n() const;
        
        /**
         * @brief Returns true if m==n.
         * 
         * @return true 
         * @return false 
         */
        bool is_square() const;

        void LU(Matrix<K>& L, Matrix<K>& U, std::vector<unsigned int>& P) const;

        /**
         * @brief Calculates determinant.
         * 
         * @return double 
         */
        K det() const;

        /**
         * @brief Calculates trace.
         * 
         * @return double 
         */
        K tr() const;

        /**
         * @brief 
         * 
         * @return std::vector<std::complex<double>> 
         */
        std::vector<std::complex<double>> eigenvalues() const;

        /**
         * @brief Returns eigenvectors.
         * 
         * May also return eigenvalues through argument.
         * 
         * @param calculate_eigenvalues 
         * @param eigenvalues 
         * @return std::vector<Vector<std::complex<double>>> 
         */
        std::vector<Vector<std::complex<double>>> eigenvectors(
            std::vector<std::complex<double>>& eigenvalues, 
            const bool calculate_eigenvalues = true) const;

        std::string to_string() const;

        protected:
        Vector<K>* vec;
        unsigned int _m;
    };
}

#endif // MATRIX_HPP
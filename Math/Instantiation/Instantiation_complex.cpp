#include "../Vector.cpp"
#include "../Matrix.cpp"

#include <complex>

namespace Math
{
    // Vector.
    template class Vector<std::complex<double>>;

    template Vector<std::complex<double>> operator+(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template Vector<std::complex<double>> operator+(const Vector<std::complex<double>>& v, const std::complex<double>& alpha);
    template Vector<std::complex<double>> operator+(const std::complex<double>& alpha, const Vector<std::complex<double>>& v);
    template Vector<std::complex<double>> operator-(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template Vector<std::complex<double>> operator-(const Vector<std::complex<double>>& v, const std::complex<double>& alpha);
    template Vector<std::complex<double>> operator-(const std::complex<double>& alpha, const Vector<std::complex<double>>& v);
    template Vector<std::complex<double>> operator*(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template Vector<std::complex<double>> operator*(const Vector<std::complex<double>>& v, const std::complex<double>& alpha);
    template Vector<std::complex<double>> operator*(const std::complex<double>& alpha, const Vector<std::complex<double>>& v);
    template Vector<std::complex<double>> operator/(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template Vector<std::complex<double>> operator/(const Vector<std::complex<double>>& v, const std::complex<double>& alpha);
    template Vector<std::complex<double>> operator/(const std::complex<double>& alpha, const Vector<std::complex<double>>& v);
    template bool operator==(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template bool operator==(const Vector<std::complex<double>>& v, const std::complex<double>& alpha);
    template bool operator==(const std::complex<double>& alpha, const Vector<std::complex<double>>& v);
    template bool operator!=(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template bool operator!=(const Vector<std::complex<double>>& v, const std::complex<double>& alpha);
    template bool operator!=(const std::complex<double>& alpha, const Vector<std::complex<double>>& v);
    template Vector<std::complex<double>> operator&(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template Vector<std::complex<double>> operator&(const Vector<std::complex<double>>& v, const std::complex<double>& alpha);
    template Vector<std::complex<double>> operator&(const std::complex<double>& alpha, const Vector<std::complex<double>>& v);
    template std::complex<double> operator|(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);    
    template std::complex<double> vector_product_2d(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template Vector<std::complex<double>> vector_product_3d(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& w);
    template std::ostream& Math::operator<< <std::complex<double>>(std::ostream&os, const Vector<std::complex<double>>& v);    

    template std::complex<double> sum(const Vector<std::complex<double>>& v);
    template std::complex<double> multiply(const Vector<std::complex<double>>& v);

    template Vector<std::complex<double>> cos(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> sin(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> tan(const Vector<std::complex<double>>& v);
    template Vector<std::complex<double>> acos(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> asin(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> atan(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> cosh(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> sinh(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> tanh(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> acosh(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> asinh(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> atanh(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> exp(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> log(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> log10(const Vector<std::complex<double>>& v);       
    template Vector<std::complex<double>> pow(const std::complex<double> v, const Vector<std::complex<double>>& exponent);    
    template Vector<std::complex<double>> pow(const Vector<std::complex<double>>& v, const std::complex<double> exponent);    
    template Vector<std::complex<double>> pow(const Vector<std::complex<double>>& v, const Vector<std::complex<double>>& exponent);    
    template Vector<std::complex<double>> sqrt(const Vector<std::complex<double>>& v);    
    template Vector<std::complex<double>> nan(const unsigned int N, const char* tagp);
    template Vector<double> fabs(const Vector<std::complex<double>>& v);    

    // Matrix.
    template class Matrix<std::complex<double>>;
    template Matrix<std::complex<double>> operator+(const Matrix<std::complex<double>>& A, const Matrix<std::complex<double>>& B);
    template Matrix<std::complex<double>> operator+(const Matrix<std::complex<double>>& A, const std::complex<double>& alpha);
    template Matrix<std::complex<double>> operator+(const std::complex<double>& alpha, const Matrix<std::complex<double>>& A);
    template Matrix<std::complex<double>> operator-(const Matrix<std::complex<double>>& A, const Matrix<std::complex<double>>& B);
    template Matrix<std::complex<double>> operator-(const Matrix<std::complex<double>>& A, const std::complex<double>& alpha);
    template Matrix<std::complex<double>> operator-(const std::complex<double>& alpha, const Matrix<std::complex<double>>& A);
    template Matrix<std::complex<double>> operator*(const Matrix<std::complex<double>>& A, const std::complex<double>& alpha);
    template Matrix<std::complex<double>> operator*(const std::complex<double>& alpha, const Matrix<std::complex<double>>& A);
    template Matrix<std::complex<double>> operator/(const Matrix<std::complex<double>>& A, const std::complex<double>& alpha);
    template Matrix<std::complex<double>> operator|(const Matrix<std::complex<double>>& A, const Matrix<std::complex<double>>& B);
    template Vector<std::complex<double>> operator|(const Matrix<std::complex<double>>& A, const Vector<std::complex<double>>& v);
    template Vector<std::complex<double>> operator|(const Vector<std::complex<double>>& v, const Matrix<std::complex<double>>& A);
}
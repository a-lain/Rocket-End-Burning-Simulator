#include "../Vector.cpp"

namespace Math
{
    template class Vector<float>;

    template Vector<float> operator+(const Vector<float>& v, const Vector<float>& w);
    template Vector<float> operator+(const Vector<float>& v, const float& alpha);
    template Vector<float> operator+(const float& alpha, const Vector<float>& v);
    template Vector<float> operator-(const Vector<float>& v, const Vector<float>& w);
    template Vector<float> operator-(const Vector<float>& v, const float& alpha);
    template Vector<float> operator-(const float& alpha, const Vector<float>& v);
    template Vector<float> operator*(const Vector<float>& v, const Vector<float>& w);
    template Vector<float> operator*(const Vector<float>& v, const float& alpha);
    template Vector<float> operator*(const float& alpha, const Vector<float>& v);
    template Vector<float> operator/(const Vector<float>& v, const Vector<float>& w);
    template Vector<float> operator/(const Vector<float>& v, const float& alpha);
    template Vector<float> operator/(const float& alpha, const Vector<float>& v);
    template bool operator<(const Vector<float>& v, const float& alpha);
    template bool operator<(const float& alpha, const Vector<float>& v);
    template bool operator<(const Vector<float>& v, const Vector<float>& w); 
    template bool operator<=(const Vector<float>& v, const float& alpha);
    template bool operator<=(const float& alpha, const Vector<float>& v);
    template bool operator<=(const Vector<float>& v, const Vector<float>& w);
    template bool operator>(const Vector<float>& v, const float& alpha);
    template bool operator>(const float& alpha, const Vector<float>& v);
    template bool operator>(const Vector<float>& v, const Vector<float>& w);
    template bool operator>=(const Vector<float>& v, const float& alpha);
    template bool operator>=(const float& alpha, const Vector<float>& v);
    template bool operator>=(const Vector<float>& v, const Vector<float>& w);
    template bool operator==(const Vector<float>& v, const Vector<float>& w);
    template bool operator==(const Vector<float>& v, const float& alpha);
    template bool operator==(const float& alpha, const Vector<float>& v);
    template bool operator!=(const Vector<float>& v, const Vector<float>& w);
    template bool operator!=(const Vector<float>& v, const float& alpha);
    template bool operator!=(const float& alpha, const Vector<float>& v);
    template Vector<float> operator&(const Vector<float>& v, const Vector<float>& w);
    template Vector<float> operator&(const Vector<float>& v, const float& alpha);
    template Vector<float> operator&(const float& alpha, const Vector<float>& v);
    template float operator|(const Vector<float>& v, const Vector<float>& w);    
    template float vector_product_2d(const Vector<float>& v, const Vector<float>& w);
    template Vector<float> vector_product_3d(const Vector<float>& v, const Vector<float>& w);
    template std::ostream& Math::operator<< <float>(std::ostream&os, const Vector<float>& v);    

    template float min(const Vector<float>& v);
    template float max(const Vector<float>& v);
    template float sum(const Vector<float>& v);
    template float multiply(const Vector<float>& v);

    template Vector<float> cos(const Vector<float>& v);    
    template Vector<float> sin(const Vector<float>& v);    
    template Vector<float> tan(const Vector<float>& v);
    template Vector<float> acos(const Vector<float>& v);    
    template Vector<float> asin(const Vector<float>& v);    
    template Vector<float> atan(const Vector<float>& v);    
    template Vector<float> atan2(const float v, const Vector<float>& w);    
    template Vector<float> atan2(const Vector<float>& v, const float w);    
    template Vector<float> atan2(const Vector<float>& v, const Vector<float>& w);    
    template Vector<float> cosh(const Vector<float>& v);    
    template Vector<float> sinh(const Vector<float>& v);    
    template Vector<float> tanh(const Vector<float>& v);    
    template Vector<float> acosh(const Vector<float>& v);    
    template Vector<float> asinh(const Vector<float>& v);    
    template Vector<float> atanh(const Vector<float>& v);    
    template Vector<float> exp(const Vector<float>& v);    
    template Vector<float> frexp(const Vector<float>& v, Vector<int>* exp);    
    template Vector<float> ldexp(const Vector<float>& v, const int exp);    
    template Vector<float> ldexp(const Vector<float>& v, const Vector<int>& exp);    
    template Vector<float> log(const Vector<float>& v);    
    template Vector<float> log10(const Vector<float>& v);    
    // template Vector<float> modf(const Vector<float>& v, Vector<float>* intpart);    
    template Vector<float> exp2(const Vector<float>& v);    
    template Vector<float> expm1(const Vector<float>& v);    
    template Vector<float> ilogb(const Vector<float>& v);    
    template Vector<float> log1p(const Vector<float>& v);    
    template Vector<float> log2(const Vector<float>& v);    
    template Vector<float> logb(const Vector<float>& v);    
    template Vector<float> scalbn(const Vector<float>& v, const int n);    
    template Vector<float> scalbn(const Vector<float>& v, const Vector<int>& n);    
    template Vector<float> scalbln(const Vector<float>& v, const long int n);    
    template Vector<float> pow(const float v, const Vector<float>& exponent);    
    template Vector<float> pow(const Vector<float>& v, const float exponent);    
    template Vector<float> pow(const Vector<float>& v, const Vector<float>& exponent);    
    template Vector<float> sqrt(const Vector<float>& v);    
    template Vector<float> cbrt(const Vector<float>& v);    
    template Vector<float> hypot(const float x, const Vector<float>& y);    
    template Vector<float> hypot(const Vector<float>& x, const float y);    
    template Vector<float> hypot(const Vector<float>& x, const Vector<float>& y);    
    template Vector<float> erf(const Vector<float>& v);    
    template Vector<float> erfc(const Vector<float>& v);    
    template Vector<float> tgamma(const Vector<float>& v);    
    template Vector<float> lgamma(const Vector<float>& v);    
    template Vector<float> ceil(const Vector<float>& v);    
    template Vector<float> floor(const Vector<float>& v);    
    template Vector<float> fmod(const float numer, const Vector<float>& denom);    
    template Vector<float> fmod(const Vector<float>& numer, const float denom);    
    template Vector<float> fmod(const Vector<float>& numer, const Vector<float>& denom);    
    template Vector<float> trunc(const Vector<float>& v);    
    template Vector<float> round(const Vector<float>& v);    
    template Vector<long int> lround(const Vector<float>& v);    
    template Vector<long long int> llround(const Vector<float>& v);    
    template Vector<float> rint(const Vector<float>& v);    
    template Vector<long int> lrint(const Vector<float>& v);
    template Vector<long long int> llrint(const Vector<float>& v);    
    template Vector<float> nearbyint(const Vector<float>& v);    
    template Vector<float> remainder(const float numer, const Vector<float>& denom);    
    template Vector<float> remainder(const Vector<float>& numer, const float denom);    
    template Vector<float> remainder(const Vector<float>& numer, const Vector<float>& denom);    
    template Vector<float> remquo(const float numer, const Vector<float>& denom, Vector<int>* quot);    
    template Vector<float> remquo(const Vector<float>& numer, const float denom, Vector<int>* quot);    
    template Vector<float> remquo(const Vector<float>& numer, const Vector<float>& denom, Vector<int>* quot);    
    template Vector<float> copysign(const Vector<float>& x, const float y);    
    template Vector<float> copysign(const Vector<float>& x, const Vector<float>& y);    
    template Vector<float> nan(const unsigned int N, const char* tagp);
    template Vector<float> nextafter(const Vector<float>& x, const Vector<float>& y);
    template Vector<float> fdim(const float x, const Vector<float>& y);    
    template Vector<float> fdim(const Vector<float>& x, float y);    
    template Vector<float> fdim(const Vector<float>& x, const Vector<float>& y);    
    template Vector<float> fabs(const Vector<float>& v);    
    template Vector<float> abs(const Vector<float>& v);    
    template Vector<float> fma(const float x, const Vector<float>& y, const Vector<float>& z);    
    template Vector<float> fma(const Vector<float>& x, const float y, const Vector<float>& z);    
    template Vector<float> fma(const Vector<float>& x, const Vector<float>& y, const float z);    
    template Vector<float> fma(const float x, const float y, const Vector<float>& z);    
    template Vector<float> fma(const float x, const Vector<float>& y, const float z);    
    template Vector<float> fma(const Vector<float>& x, const float y, const float z);    
    template Vector<float> fma(const Vector<float>& x, const Vector<float>& y, const Vector<float>& z);    
    template Vector<int> fpclassify(const Vector<float>& v);    
    template bool isfinite(const Vector<float>& v);    
    template bool isinf(const Vector<float>& v);    
    template bool isnan(const Vector<float>& v);    
    template bool isnormal(const Vector<float>& v);    
}

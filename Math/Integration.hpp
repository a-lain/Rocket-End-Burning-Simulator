#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include <functional>

namespace Math
{
    namespace Integrators
    {
        double rectangle_rule(std::function<double(double x)> f, const double a, const double b,
            const unsigned int N);

        double trapezoidal_rule(std::function<double(double x)> f, const double a, const double b,
            const unsigned int N);

        double Gauss_Konrad_G7_K15(std::function<double(double x)> f, const double a, const double b,
            const double tol = 1e-6);
    }
}

#endif // INTEGRATION_HPP
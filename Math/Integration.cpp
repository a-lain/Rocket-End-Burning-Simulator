#ifndef INTEGRATION_CPP
#define INTEGRATION_CPP

#include "Integration.hpp"
#include <cmath>

using namespace Math::Integrators;

double Math::Integrators::rectangle_rule(std::function<double(double x)> f, const double a, const double b,
    const unsigned int N)
{
    double delta_x = (b - a) / N;
    double I = 0;
    for (unsigned int i = 0; i < N; i++)
    {
        I += f(a + delta_x*(1./2 + i));
    }
    return I * delta_x;
}

double Math::Integrators::trapezoidal_rule(std::function<double(double x)> f, const double a, const double b,
    const unsigned int N)
{
    double delta_x = (b - a) / N;
    double I = 0;
    for (unsigned int i = 1; i < N; i++)
    {
        I += f(a + i*delta_x);
    }
    return (I + (f(a) + f(b)) / 2) * delta_x;
}

double Math::Integrators::Gauss_Konrad_G7_K15(std::function<double(double x)> f, const double a, const double b, const double tol)
{
    static const double nodes[15] = {
                        -0.9914553711208126392069,
                        -0.9491079123427585245262,
                        -0.8648644233597690727897,
                        -0.7415311855993944398639,
                        -0.5860872354676911302941,
                        -0.4058451513773971669066,
                        -0.2077849550078984676007,
                        0,
                        0.2077849550078984676007,
                        0.4058451513773971669066,
                        0.5860872354676911302941,
                        0.7415311855993944398639,
                        0.8648644233597690727897,
                        0.9491079123427585245262,
                        0.9914553711208126392069
                        };
    static const double weights_K_15[15] = {
                        0.0229353220105292249637,
                        0.0630920926299785532907,
                        0.1047900103222501838399,
                        0.140653259715525918745,
                        0.1690047266392679028266,
                        0.1903505780647854099133,
                        0.2044329400752988924142,
                        0.209482141084727828013,
                        0.2044329400752988924142,
                        0.190350578064785409913,
                        0.1690047266392679028266,
                        0.1406532597155259187452,
                        0.10479001032225018384,
                        0.0630920926299785532907,
                        0.02293532201052922496373
                        };
    const double weights_G_7[15] = {
                        0,
                        0.129484966168869693271,
                        0,
                        0.279705391489276667901,
                        0,
                        0.3818300505051189449504,
                        0,
                        0.4179591836734693877551,
                        0,
                        0.3818300505051189449504,
                        0,
                        0.2797053914892766679015,
                        0,
                        0.1294849661688696932706,
                        0
                        };

    double G7 = 0;
    double K15 = 0;
    double half_sum = (b + a) / 2;
    double half_substract = (b - a) / 2;
    double fx = 0;
    for (unsigned int i = 0; i < 15; i++)
    {
        fx = f(half_sum + half_substract*nodes[i]);
        G7 += weights_G_7[i] * fx;
        K15 += weights_K_15[i] * fx;
    }
    G7 *= half_substract;
    K15 *= half_substract;
    
    if(fabs(K15 - G7) > tol)
    {
        return Gauss_Konrad_G7_K15(f, a, half_sum, tol / 2) +
            Gauss_Konrad_G7_K15(f, half_sum, b, tol / 2);
    }
    else
    {
        return K15;
    }
}

#endif // INTEGRATION_CPP
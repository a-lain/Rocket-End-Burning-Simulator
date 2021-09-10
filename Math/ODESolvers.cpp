#ifndef ODE_SOLVERS_CPP
#define ODE_SOLVERS_CPP

#include "ODESolvers.hpp"
#include <iostream>
#include <cmath>

using namespace Math;

Vector<double> Math::ODESolvers::euler_explicit(const Vector<double>& x, const double t,
    std::function<Vector<double>(const Vector<double>&, const double)> f, const double dt)
{
    return x + f(x, t) * dt;
}

Vector<double> Math::ODESolvers::euler_explicit(const Vector<double>& x,
    std::function<Vector<double>(const Vector<double>&)> f, const double dt)
{
    return x + f(x) * dt;
}

Vector<double> Math::ODESolvers::runge_kutta_order_4(const Vector<double>& x, const double t,
    std::function<Vector<double>(const Vector<double>&, const double)> f, const double dt)
{
    Vector<double> k1, k2, k3, k4;
    k1 = f(x, t);
    k2 = f(x + 1./2*k1*dt, t + dt/2);
    k3 = f(x + 1./2*k2*dt, t + dt/2);
    k4 = f(x + k3*dt, t + dt);
    return x + 1./6*dt*(k1 + 2.*k2 + 2.*k3 + k4);
}

Vector<double> Math::ODESolvers::runge_kutta_order_4(const Vector<double>& x,
    std::function<Vector<double>(const Vector<double>&)> f, const double dt)
{
    Vector<double> k1, k2, k3, k4;
    k1 = f(x);
    k2 = f(x + 1./2*k1*dt);
    k3 = f(x + 1./2*k2*dt);
    k4 = f(x + k3*dt);
    return x + 1./6*dt*(k1 + 2.*k2 + 2.*k3 + k4);
}

void Math::ODESolvers::runge_kutta_fehlberg(Vector<double>& x, double& t, double& dt, 
    std::function<Vector<double>(const Vector<double>&, const double)> f,
    const Vector<double>& abs_tol)
{
    Vector<double> x_old = x;
    double t_old = t;
    // We define the constantes of the method.
    const double a1 = 0, a2 = 1./4, a3 = 3./8, a4 = 12./13, a5 = 1., a6 = 1./2;
    const double b21 = 1./4, b31 = 3./32, b41 = 1932./2197, b51 = 439./216, b61 = -8./27;
    const double b32 = 9./32, b42 = -7200./2197, b52 = -8., b62 = 2.;
    const double b43 = 7296./2197, b53 = 3680./513, b63 = -3544./2565;
    const double b54 = -845./4104, b64 = 1859./4104;
    const double b65 = -11./40;
    const double ch1 = 16./135, ch2 = 0., ch3 = 6656./12825, ch4 = 28561./56430, ch5 = -9./50, ch6 = 2./55;
    const double ct1 = 1./360, ct2 = 0., ct3 = -128./4275, ct4 = -2197./75240, ct5 = 1./50, ct6 = 2./55;

    Vector<double> k1, k2, k3, k4, k5, k6;
    /*! Truncation error.*/
    Vector<double> TE;
    k1 = dt * f(x, t + a1*dt);
    k2 = dt * f(x + b21*k1, t + a2*dt);
    k3 = dt * f(x + b31*k1 + b32*k2, t + a3*dt);
    k4 = dt * f(x + b41*k1 + b42*k2 + b43*k3, t + a4*dt);
    k5 = dt * f(x + b51*k1 + b52*k2 + b53*k3 + b54*k4, t + a5*dt);
    k6 = dt * f(x + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5, t + a6*dt);

    x = x_old + ch1*k1 + ch2*k2 + ch3*k3 + ch4*k4 + ch5*k5 + ch6*k6;
    TE = fabs(ct1*k1 + ct2*k2 + ct3*k3 + ct4*k4 + ct5*k5 + ct6*k6);
    t += dt;
    double quotient_min = min(abs_tol / TE);
    dt = 0.9 * dt * ::pow(quotient_min, 1./5);

    // If the truncation error exceeds are absolute tolerance, 
    if (quotient_min < 1)
    {
        x = x_old;
        t = t_old;
        runge_kutta_fehlberg(x, t, dt, f, abs_tol);
    }
}

void Math::ODESolvers::runge_kutta_fehlberg(Vector<double>& x, double& t, double& dt, 
    std::function<Vector<double>(const Vector<double>&)> f,
    const Vector<double>& abs_tol)
{
    Vector<double> x_old = x;
    double t_old = t;
    // We define the constantes of the method.
    const double b21 = 1./4, b31 = 3./32, b41 = 1932./2197, b51 = 439./216, b61 = -8./27;
    const double b32 = 9./32, b42 = -7200./2197, b52 = -8., b62 = 2.;
    const double b43 = 7296./2197, b53 = 3680./513, b63 = -3544./2565;
    const double b54 = -845./4104, b64 = 1859./4104;
    const double b65 = -11./40;
    const double ch1 = 16./135, ch2 = 0., ch3 = 6656./12825, ch4 = 28561./56430, ch5 = -9./50, ch6 = 2./55;
    const double ct1 = 1./360, ct2 = 0., ct3 = -128./4275, ct4 = -2197./75240, ct5 = 1./50, ct6 = 2./55;

    Vector<double> k1, k2, k3, k4, k5, k6;
    /*! Truncation error.*/
    Vector<double> TE;
    k1 = dt * f(x);
    k2 = dt * f(x + b21*k1);
    k3 = dt * f(x + b31*k1 + b32*k2);
    k4 = dt * f(x + b41*k1 + b42*k2 + b43*k3);
    k5 = dt * f(x + b51*k1 + b52*k2 + b53*k3 + b54*k4);
    k6 = dt * f(x + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5);

    x = x_old + ch1*k1 + ch2*k2 + ch3*k3 + ch4*k4 + ch5*k5 + ch6*k6;
    TE = fabs(ct1*k1 + ct2*k2 + ct3*k3 + ct4*k4 + ct5*k5 + ct6*k6);
    t += dt;
    double quotient_min = min(abs_tol / TE);
    dt = 0.9 * dt * ::pow(quotient_min, 1./5);

    // If the truncation error exceeds are absolute tolerance, 
    if (quotient_min < 1)
    {
        x = x_old;
        t = t_old;
        runge_kutta_fehlberg(x, t, dt, f, abs_tol);
    }
}

#endif // ODE_SOLVERS_CPP
#ifndef ALGEBRAIC_SOLVERS_CPP
#define ALGEBRAIC_SOLVERS_CPP

#include "AlgebraicSolvers.hpp"
#include <cmath>
#include <iostream>
#include <limits>

using namespace Math;
using namespace AlgebraicSolvers;

double AlgebraicSolvers::bisection(std::function<double(const double x)> f,
        double a, double b,
        const unsigned int iter_max, const double abs_tol, const double rel_tol)
{
    double fa = f(a);
    double old_p = a;
    double p;
    double fp;
    unsigned int N_iter = 0;

    // We start the algorithm.
    p = a + (b - a) / 2;
    fp = f(p);
    
    // We check whether we have to continue.
    while((b - a)/2 > abs_tol && ::fabs(p - old_p) / ::fabs(p) > rel_tol && N_iter < iter_max && fp != 0)
    {
        if (fa * fp > 0)
        {
            a = p;
            fa = fp;
        }
        else
        {
            b = p;
        }

        old_p = p;
        p = a + (b - a) / 2;
        fp = f(p);
        N_iter++;
    }
    return p;
}

double AlgebraicSolvers::secant(std::function<double(const double x)> f,
        double x1, double x2,
        const unsigned int iter_max, const double abs_tol, const double rel_tol)
{
    double f1 = f(x1);
    double f2 = f(x2);
    double x;
    unsigned int N_iter = 0;

    // We start the algorithm.
    x = x2 - f2*(x2 - x1) / (f2 - f1);

    // We check whether we have to continue.
    while (::fabs(x - x2) > abs_tol && ::fabs(x - x2) / ::fabs(x) > rel_tol && N_iter < iter_max)
    {        
        x1 = x2;
        f1 = f2;
        x2 = x;
        f2 = f(x);
        x = x2 - f2*(x2 - x1) / (f2 - f1);
        N_iter++;
    }
    return x;
}

double AlgebraicSolvers::newton_raphson(std::function<double(const double x)>f,
        std::function<double(const double x)>dfdx, const double x0,
        const unsigned int iter_max, const double abs_tol, const double rel_tol)
{
    // We start the algorithm.
    double x_old = x0;
    double x = x_old - f(x_old) / dfdx(x_old);
    unsigned int N_iter = 0;

    while(::fabs(x - x_old) > abs_tol && ::fabs(x - x_old) / ::fabs(x) > rel_tol && N_iter < iter_max)
    {
        x_old = x;
        x = x_old - f(x_old) / dfdx(x_old);
        N_iter++;
    }
    return x;
}

double AlgebraicSolvers::newton_raphson(std::function<DualNumber<double>(const DualNumber<double> x)>f,
        const double x0,
        const unsigned int iter_max, const double abs_tol, const double rel_tol)
{
    // We start the algorithm.
    double x_old = x0;
    DualNumber<double> fx_old = f(x_old + DualNumber<double>::epsilon);
    double x = x_old - fx_old.a / fx_old.b;
    unsigned int N_iter = 0;

    while(::fabs(x - x_old) > abs_tol && ::fabs(x - x_old) / ::fabs(x) > rel_tol && N_iter < iter_max)
    {
        x_old = x;
        fx_old = f(x_old + DualNumber<double>::epsilon);
        x = x_old - fx_old.a / fx_old.b;
        N_iter++;
    }
    return x;
}

LinearSystemSolver::LinearSystemSolver(const Matrix<double>& A)
{
    A.LU(L, U, P);
}

Vector<double> LinearSystemSolver::solve(const Vector<double>& b) const
{
    if (b.size() != U.n())
    {
        throw std::runtime_error("Coefficient matrix and independent vector must have the same number of rows.\n");
    }
    else
    {
        if (b.norm_2() < 1e-12)
        {
            return Vector<double>(b.size());
        }
        else
        {
            // Step I: We permute the rows of b. This is the same as doing b=Pb.
            Vector<double> b_new = Vector<double>(b.size());
            for (unsigned int i = 0; i < b.size(); i++)
            {
                b_new[i] = b[P[i]];
            }

            // Step II: We apply progressive substitution. We are solving Lz=Pb.
            Vector<double> z = Vector<double>(b.size());
            for (unsigned int i = 0; i < b.size(); i++)
            {
                double sum = 0;
                for (unsigned int j = 0; j < i; j++)
                {
                    sum += L(i,j)*z[j];
                }
                z[i] = b_new[i] - sum;
            }

            // Step III: We apply regressive substitution. We are solving Ux=z.
            Vector<double> x = Vector<double>(b.size());
            for (unsigned int i = b.size() - 1; i != std::numeric_limits<unsigned int>::max(); i--)
            {
                double sum = 0;
                for (unsigned int j = i+1; j < b.size(); j++)
                {
                    sum += U(i,j)*x[j];
                }
                x[i] = (z[i] - sum) / U(i,i);
            }

            return x;
        }        
    }
}

#endif // ALGEBRAIC_SOLVERS_CPP
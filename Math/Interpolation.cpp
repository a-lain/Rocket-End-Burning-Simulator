#ifndef INTERPOLATION_CPP
#define INTERPOLATION_CPP

#include "Interpolation.hpp"
#include <stdexcept>
#include "AlgebraicSolvers.hpp"
#include <limits>

#include <iostream>

using namespace Math::Interpolation;

AverageLinearInterpolation::AverageLinearInterpolation(double* averages,
    double *partition, unsigned int N_cells)
{
    this->N_cells = N_cells;
    x_part = new double[N_cells + 1];
    for (unsigned int i = 0; i < N_cells + 1; i++)
    {
        x_part[i] = partition[i];
    }

    b = new double[N_cells - 1];
    c = new double[N_cells - 1];

    for (unsigned int i = 0; i < N_cells - 1; i++)
    {
        b[i] = 2*(averages[i+1] - averages[i]) / (x_part[i+2] - x_part[i]);
        c[i] = averages[i] - (x_part[i] + x_part[i+1]) / 2 * b[i];
    }
}

AverageLinearInterpolation::AverageLinearInterpolation(const std::vector<double>& averages,
    const std::vector<double>& partition)
{
    this->N_cells = averages.size();
    x_part = new double[N_cells + 1];
    for (unsigned int i = 0; i < N_cells + 1; i++)
    {
        x_part[i] = partition[i];
    }

    b = new double[N_cells - 1];
    c = new double[N_cells - 1];

    for (unsigned int i = 0; i < N_cells - 1; i++)
    {
        b[i] = 2*(averages[i+1] - averages[i]) / (x_part[i+2] - x_part[i]);
        c[i] = averages[i] - (x_part[i] + x_part[i+1]) / 2 * b[i];
    }
}

AverageLinearInterpolation::AverageLinearInterpolation(const AverageLinearInterpolation& f)
{
    N_cells = f.N_cells;
    x_part = new double[f.N_cells + 1];
    b = new double[f.N_cells];
    c = new double[f.N_cells];

    for (unsigned int i = 0; i < N_cells; i++)
    {
        x_part[i] = f.x_part[i];
        b[i] = f.b[i];
        c[i] = f.c[i];
    }
    x_part[N_cells] = f.x_part[N_cells];
}

AverageLinearInterpolation& AverageLinearInterpolation::operator=(const AverageLinearInterpolation& f)
{
    if (&f != this)
    {
        delete[] x_part;
        delete[] b;
        delete[] c;
        N_cells = f.N_cells;
        x_part = new double[f.N_cells + 1];
        b = new double[f.N_cells];
        c = new double[f.N_cells];

        for (unsigned int i = 0; i < N_cells; i++)
        {
            x_part[i] = f.x_part[i];
            b[i] = f.b[i];
            c[i] = f.c[i];
        }
        x_part[N_cells] = f.x_part[N_cells];
    }
    return *this;
}

AverageLinearInterpolation::~AverageLinearInterpolation()
{
    delete[] x_part;
    delete[] b;
    delete[] c;
}

double AverageLinearInterpolation::operator()(const double x) const
{
    if (x < (x_part[0] + x_part[1]) / 2)
    {
        return b[0]*x + c[0];
    }
    else if (x > (x_part[N_cells - 1] + x_part[N_cells]) / 2)
    {
        return b[N_cells - 2]*x + c[N_cells - 2];
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = N_cells - 1;
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (x > (x_part[i+1] + x_part[i+2]) / 2)
            {
                left = i+1;
            }
            else if (x < (x_part[i] + x_part[i+1]) / 2)
            {
                right = i-1;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }
        return b[i]*x + c[i];
    }
}

AverageQuadraticInterpolation::AverageQuadraticInterpolation(double* averages,
    double *partition, unsigned int N_cells)
{
    this->N_cells = N_cells;
    x_part = new double[N_cells + 1];
    for (unsigned int i = 0; i < N_cells + 1; i++)
    {
        x_part[i] = partition[i];
    }

    a = new double[N_cells];
    b = new double[N_cells];
    c = new double[N_cells];
    double* alpha = new double[3*N_cells - 2];
    double* beta = new double[3*N_cells - 2];
    double* gamma = new double[3*N_cells - 2];
    double* delta = new double[3*N_cells - 2];
    double* epsilon = new double[3*N_cells - 2];
    double* indep = new double[3*N_cells - 2];

    // We calculate the matrix diagonals.
    alpha[0] = 0;
    alpha[1] = 1;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        alpha[2+3*i] = (2*x_part[i+2]*x_part[i+1] -x_part[i+1]*x_part[i+1] -x_part[i+2]*x_part[i+2]) / 3;
        alpha[3+3*i] = 1;
        alpha[4+3*i] = 1;
    }
    alpha[3*N_cells - 4] = 0;
    alpha[3*N_cells - 3] = 0;

    beta[0] = (x_part[0] + x_part[1]) / 2;
    beta[1] = -1;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        beta[2+3*i] = 2*x_part[i+1];
        beta[3+3*i] = (x_part[i+2] + x_part[i+1]) / 2;
        beta[4+3*i] = -1;
    }
    beta[3*N_cells - 4] = 1;
    beta[3*N_cells - 3] = 1;

    gamma[0] = -(x_part[2] + x_part[1])/2;
    gamma[1] = 0;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        gamma[2+3*i] = (x_part[i+2]*x_part[i+2] + x_part[i+1]*x_part[i+1] + x_part[i+2]*x_part[i+1]) / 3;
        gamma[3+3*i] = -(x_part[i+3] + x_part[i+2]) / 2;
        gamma[4+3*i] = 0;
    }
    gamma[3*N_cells - 4] = (x_part[N_cells] + x_part[N_cells - 1]) / 2;
    gamma[3*N_cells - 3] = 0;

    delta[0] = -1;
    delta[1] = 0;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        delta[2+3*i] = -x_part[i+3]*x_part[i+2];
        delta[3+3*i] = -1;
        delta[4+3*i] = 0;
    }
    delta[3*N_cells - 4] = 0;
    delta[3*N_cells - 3] = 0;

    epsilon[0] = 0;
    epsilon[1] = 0;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        epsilon[2+3*i] = -2*x_part[i+2];
        epsilon[3+3*i] = 0;
        epsilon[4+3*i] = 0;
    }
    epsilon[3*N_cells - 4] = 0;
    epsilon[3*N_cells - 3] = 0;

    for(unsigned int i = 0; i < N_cells - 1; i++)
    {
        indep[3*i] = averages[i];
        indep[3*i+1] = -averages[i+1];
        indep[3*i+2] = 0;
    }
    indep[3*N_cells - 3] = averages[N_cells - 1];

    // We apply Gaussian Elimination.
    double l = 1;
    for (unsigned int j = 3*N_cells - 5; j >= 3; j--)
    {
        l = alpha[j] / beta[j];
        alpha[j] = 0;
        beta[j-1] = beta[j-1] - l*gamma[j-1];
        gamma[j-2] = gamma[j-2] - l*delta[j-2];
        delta[j-3] = delta[j-3] - l*epsilon[j-3];
        indep[j-1] = indep[j-1] - l*indep[j];
    }
    l = alpha[2] / beta[2];
    alpha[2] = 0;
    beta[1] = beta[1] - l*gamma[1];
    gamma[0] = gamma[0] - l*delta[0];
    indep[1] = indep[1] - l*indep[2];
    l = alpha[1] / beta[1];
    alpha[1] = 0;
    beta[0] = beta[0] - l*gamma[0];
    indep[0] = indep[0] - l*indep[1];

    // We apply regressive substitution.
    a[0] = 0;
    b[0] = indep[0] / beta[0];
    c[0] = (indep[1] - gamma[0]*b[0]) / beta[1];
    a[1] = (indep[2] - gamma[1]*c[0] - delta[0]*b[0]) / beta[2];
    b[1] = (indep[3] - gamma[2]*a[1] - delta[1]*c[0] - epsilon[0]*b[0]) / beta[3];
    c[1] = (indep[4] - gamma[3]*b[1] - delta[2]*a[1] - epsilon[1]*c[0]) / beta[4];
    unsigned int j = 5;
    for (unsigned int i = 2; i < N_cells - 1; i++)
    {
        a[i] = (indep[j] - gamma[j-1]*c[i-1] - delta[j-2]*b[i-1] - epsilon[j-3]*a[i-1]) / beta[j];
        j++;
        b[i] = (indep[j] - gamma[j-1]*a[i] - delta[j-2]*c[i-1] - epsilon[j-3]*b[i-1]) / beta[j];
        j++; 
        c[i] = (indep[j] - gamma[j-1]*b[i] - delta[j-2]*a[i] - epsilon[j-3]*c[i-1]) / beta[j];
        j++; 
    }
    a[N_cells - 1] = 0;
    b[N_cells - 1] = 2*x_part[N_cells - 1]*a[N_cells - 2] + b[N_cells - 2];
    c[N_cells - 1] = averages[N_cells - 1] - (x_part[N_cells] + x_part[N_cells-1]) / 2 * b[N_cells - 1];

    delete[] alpha;
    delete[] beta;
    delete[] gamma;
    delete[] delta;
    delete[] epsilon;
    delete[] indep;
}

AverageQuadraticInterpolation::AverageQuadraticInterpolation(const std::vector<double>& averages,
    const std::vector<double>& partition)
{
    this->N_cells = averages.size();
    x_part = new double[N_cells + 1];
    for (unsigned int i = 0; i < N_cells + 1; i++)
    {
        x_part[i] = partition[i];
    }

    a = new double[N_cells];
    b = new double[N_cells];
    c = new double[N_cells];
    double* alpha = new double[3*N_cells - 2];
    double* beta = new double[3*N_cells - 2];
    double* gamma = new double[3*N_cells - 2];
    double* delta = new double[3*N_cells - 2];
    double* epsilon = new double[3*N_cells - 2];
    double* indep = new double[3*N_cells - 2];

    // We calculate the matrix diagonals.
    alpha[0] = 0;
    alpha[1] = 1;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        alpha[2+3*i] = (2*x_part[i+2]*x_part[i+1] -x_part[i+1]*x_part[i+1] -x_part[i+2]*x_part[i+2]) / 3;
        alpha[3+3*i] = 1;
        alpha[4+3*i] = 1;
    }
    alpha[3*N_cells - 4] = 0;
    alpha[3*N_cells - 3] = 0;

    beta[0] = (x_part[0] + x_part[1]) / 2;
    beta[1] = -1;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        beta[2+3*i] = 2*x_part[i+1];
        beta[3+3*i] = (x_part[i+2] + x_part[i+1]) / 2;
        beta[4+3*i] = -1;
    }
    beta[3*N_cells - 4] = 1;
    beta[3*N_cells - 3] = 1;

    gamma[0] = -(x_part[2] + x_part[1])/2;
    gamma[1] = 0;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        gamma[2+3*i] = (x_part[i+2]*x_part[i+2] + x_part[i+1]*x_part[i+1] + x_part[i+2]*x_part[i+1]) / 3;
        gamma[3+3*i] = -(x_part[i+3] + x_part[i+2]) / 2;
        gamma[4+3*i] = 0;
    }
    gamma[3*N_cells - 4] = (x_part[N_cells] + x_part[N_cells - 1]) / 2;
    gamma[3*N_cells - 3] = 0;

    delta[0] = -1;
    delta[1] = 0;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        delta[2+3*i] = -x_part[i+3]*x_part[i+2];
        delta[3+3*i] = -1;
        delta[4+3*i] = 0;
    }
    delta[3*N_cells - 4] = 0;
    delta[3*N_cells - 3] = 0;

    epsilon[0] = 0;
    epsilon[1] = 0;
    for (unsigned int i = 0; i < N_cells - 2; i++)
    {
        epsilon[2+3*i] = -2*x_part[i+2];
        epsilon[3+3*i] = 0;
        epsilon[4+3*i] = 0;
    }
    epsilon[3*N_cells - 4] = 0;
    epsilon[3*N_cells - 3] = 0;

    for(unsigned int i = 0; i < N_cells - 1; i++)
    {
        indep[3*i] = averages[i];
        indep[3*i+1] = -averages[i+1];
        indep[3*i+2] = 0;
    }
    indep[3*N_cells - 3] = averages[N_cells - 1];

    // We apply Gaussian Elimination.
    double l = 1;
    for (unsigned int j = 3*N_cells - 5; j >= 3; j--)
    {
        l = alpha[j] / beta[j];
        alpha[j] = 0;
        beta[j-1] = beta[j-1] - l*gamma[j-1];
        gamma[j-2] = gamma[j-2] - l*delta[j-2];
        delta[j-3] = delta[j-3] - l*epsilon[j-3];
        indep[j-1] = indep[j-1] - l*indep[j];
    }
    l = alpha[2] / beta[2];
    alpha[2] = 0;
    beta[1] = beta[1] - l*gamma[1];
    gamma[0] = gamma[0] - l*delta[0];
    indep[1] = indep[1] - l*indep[2];
    l = alpha[1] / beta[1];
    alpha[1] = 0;
    beta[0] = beta[0] - l*gamma[0];
    indep[0] = indep[0] - l*indep[1];

    // We apply regressive substitution.
    a[0] = 0;
    b[0] = indep[0] / beta[0];
    c[0] = (indep[1] - gamma[0]*b[0]) / beta[1];
    a[1] = (indep[2] - gamma[1]*c[0] - delta[0]*b[0]) / beta[2];
    b[1] = (indep[3] - gamma[2]*a[1] - delta[1]*c[0] - epsilon[0]*b[0]) / beta[3];
    c[1] = (indep[4] - gamma[3]*b[1] - delta[2]*a[1] - epsilon[1]*c[0]) / beta[4];
    unsigned int j = 5;
    for (unsigned int i = 2; i < N_cells - 1; i++)
    {
        a[i] = (indep[j] - gamma[j-1]*c[i-1] - delta[j-2]*b[i-1] - epsilon[j-3]*a[i-1]) / beta[j];
        j++;
        b[i] = (indep[j] - gamma[j-1]*a[i] - delta[j-2]*c[i-1] - epsilon[j-3]*b[i-1]) / beta[j];
        j++; 
        c[i] = (indep[j] - gamma[j-1]*b[i] - delta[j-2]*a[i] - epsilon[j-3]*c[i-1]) / beta[j];
        j++; 
    }
    a[N_cells - 1] = 0;
    b[N_cells - 1] = 2*x_part[N_cells - 1]*a[N_cells - 2] + b[N_cells - 2];
    c[N_cells - 1] = averages[N_cells - 1] - (x_part[N_cells] + x_part[N_cells-1]) / 2 * b[N_cells - 1];

    delete[] alpha;
    delete[] beta;
    delete[] gamma;
    delete[] delta;
    delete[] epsilon;
    delete[] indep;
}

AverageQuadraticInterpolation::AverageQuadraticInterpolation(const AverageQuadraticInterpolation& f)
{
    N_cells = f.N_cells;
    x_part = new double[f.N_cells + 1];
    a = new double[f.N_cells];
    b = new double[f.N_cells];
    c = new double[f.N_cells];

    for (unsigned int i = 0; i < N_cells; i++)
    {
        x_part[i] = f.x_part[i];
        a[i] = f.a[i];
        b[i] = f.b[i];
        c[i] = f.c[i];
    }
    x_part[N_cells] = f.x_part[N_cells];
}

AverageQuadraticInterpolation& AverageQuadraticInterpolation::operator=(const AverageQuadraticInterpolation& f)
{
    if (&f != this)
    {
        delete[] x_part;
        delete[] a;
        delete[] b;
        delete[] c;
        N_cells = f.N_cells;
        x_part = new double[f.N_cells + 1];
        a = new double[f.N_cells];
        b = new double[f.N_cells];
        c = new double[f.N_cells];

        for (unsigned int i = 0; i < N_cells; i++)
        {
            x_part[i] = f.x_part[i];
            a[i] = f.a[i];
            b[i] = f.b[i];
            c[i] = f.c[i];
        }
        x_part[N_cells] = f.x_part[N_cells];
    }
    return *this;
}

AverageQuadraticInterpolation::~AverageQuadraticInterpolation()
{
    delete[] x_part;
    delete[] a;
    delete[] b;
    delete[] c;
}

double AverageQuadraticInterpolation::operator()(const double x) const
{
    if (x < x_part[0])
    {
        return b[0]*x + c[0];
    }
    else if (x > x_part[N_cells])
    {
        return b[N_cells - 1]*x + c[N_cells - 1];
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = N_cells;
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (x > x_part[i+1])
            {
                left = i+1;
            }
            else if (x < x_part[i])
            {
                right = i-1;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }
        return a[i]*x*x + b[i]*x + c[i];
    }
}

LinearInterpolation::LinearInterpolation(const std::vector<double>& Y, const std::vector<double>& X)
{
    this->N = Y.size();
    b = new double[N-1];
    c = new double[N-1];
    x_part = new double[N];

    for (unsigned int i = 0; i < N - 1; i++)
    {
        x_part[i] = X[i];
        b[i] = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
        c[i] = Y[i] - b[i]*X[i];
    }
    x_part[N - 1] = X[N - 1];
}

LinearInterpolation::LinearInterpolation(const double* Y, const double* X, const unsigned int N)
{
    this->N = N;
    b = new double[N-1];
    c = new double[N-1];
    x_part = new double[N];

    for (unsigned int i = 0; i < N - 1; i++)
    {
        x_part[i] = X[i];
        b[i] = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
        c[i] = Y[i] - b[i]*X[i];
    }
    x_part[N - 1] = X[N - 1];
}

LinearInterpolation::LinearInterpolation(const LinearInterpolation& I)
{
    this->N = I.N;
    b = new double[N-1];
    c = new double[N-1];
    x_part = new double[N];

    for (unsigned int i = 0; i < N - 1; i++)
    {
        x_part[i] = I.x_part[i];
        b[i] = I.b[i];
        c[i] = I.c[i];
    }
    x_part[N - 1] = I.x_part[N - 1];
}

LinearInterpolation& LinearInterpolation::operator=(const LinearInterpolation& I)
{
    if (this != &I)
    {
        delete[] b;
        delete[] c;
        delete[] x_part;

        this->N = I.N;
        b = new double[N-1];
        c = new double[N-1];
        x_part = new double[N];

        for (unsigned int i = 0; i < N - 1; i++)
        {
            x_part[i] = I.x_part[i];
            b[i] = I.b[i];
            c[i] = I.c[i];
        }
        x_part[N - 1] = I.x_part[N - 1];
    }   
    return *this; 
}

LinearInterpolation::~LinearInterpolation()
{
    delete[] x_part;
    delete[] b;
    delete[] c;
}

double LinearInterpolation::operator()(const double x) const
{
    if (x < x_part[0])
    {
        return b[0]*x + c[0];
    }
    else if (x > x_part[N - 1])
    {
        return b[N - 2]*x + c[N - 2];
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = N - 1;
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (x > x_part[i+1])
            {
                left = i+1;
            }
            else if (x < x_part[i])
            {
                right = i-1;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }
        return b[i]*x + c[i];
    }
}

#endif // INTERPOLATION_CPP
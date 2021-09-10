#ifndef EXACT_STEADY_SOLVER_CPP
#define EXACT_STEADY_SOLVER_CPP

#include "ExactSteadySolver.hpp"
#include <cmath>
#include "../../Math/AlgebraicSolvers.hpp"

using namespace Solvers::Gas;

double ExactSteadySolver::T(const double A) const
{
    double _M = M(A);
    return T_0 / (G1*_M*_M + 1);
}

double ExactSteadySolver::P(const double A) const
{
    double _M = M(A);
    return P_0 / pow(G1*_M*_M + 1, G2);
}

double ExactSteadySolver::rho(const double A) const
{
    double _M = M(A);
    return rho_0 / pow(G1*_M*_M + 1, G3);
}

double ExactSteadySolver::v(const double A) const
{
    return M(A)*sqrt(gamma * R * T(A));
}

double ExactSteadySolver::M(const double A) const
{
    auto f = [=](const double x)
    {
        return x / pow(1 + G1*x*x, G4) * A - C;
    };

    if (solution_type == SolutionType::SUBSONIC)
    {
        return Math::AlgebraicSolvers::bisection(f, 0, 1, 1000, 1e-6, 1e-6);
    }

    else
    {
        if (f(sqrt(3)) > 0)
        {
            auto dfdx = [=](const double x)
            {
                return -A*(x*x - 1) / (pow(G1*x*x + 1, G4) * (G1*x*x + 1));
            };
            return Math::AlgebraicSolvers::newton_raphson(f, dfdx, sqrt(3), 1000, 1e-6, 1e-6);
        }
        else
        {
            return Math::AlgebraicSolvers::bisection(f, 1, sqrt(3), 1000, 1e-6, 1e-6);
        }
    }
}

ExactSteadySolver::ExactSteadySolver(const double rho_0, const double P_0, const double T_0,
    const double M_x, const double A_x, const double gamma, const SolutionType solution_type):
        solution_type(solution_type), rho_0(rho_0), P_0(P_0), T_0(T_0), gamma(gamma)
{
    R = P_0 / (rho_0 * T_0);

    G1 = (gamma - 1) / 2;
    G2 = gamma / (gamma - 1);
    G3 = 1 / (gamma - 1);
    G4 = (gamma + 1) / (2*(gamma - 1));

    C = M_x / pow(1 + G1*M_x*M_x, G4) * A_x;
}

#endif // EXACT_STEADY_SOLVER_CPP
#ifndef ROCKET_SOLVER_CPP
#define ROCKET_SOLVER_CPP

#include "RocketSolver.hpp"
#include <cmath>
#include "../../Math/AlgebraicSolvers.hpp"

using namespace Solvers::Rocket;

SteadySolver::SteadySolver(const double T_c, const double M_catm, const double n, 
    const double gamma, const double R, const double A_c):
    A_c(A_c), M_catm(M_catm), gamma(gamma), R(R), n(n), T_c(T_c)
{
    G1 = (gamma - 1) / 2;
    G2 = gamma / (gamma - 1);
    G3 = (gamma + 1) / 2;
    G4 = (gamma + 1) / (2*(gamma - 1));
}

void SteadySolver::solve_for_exit_area(const double A_e)
{
    this->A_e = A_e;
    auto g = [=](double P_c)
    {
        // We calculate M_c.
        double M_c = M_catm * pow(P_c/P_atm, n - 1);

        // We calculate M_e.
        double C = M_c / pow(1 + G1*M_c*M_c, G4) * A_c;
        auto f = [=](const double x)
        {
            return x / pow(1 + G1*x*x, G4) * A_e - C;
        };
        double M_e = Math::AlgebraicSolvers::bisection(f, 0, 1, 1000, 1e-6, 1e-6);
        double P_e = P_c * pow((G1*M_c*M_c+ 1) / (G1*M_e*M_e + 1), G2);
        return P_e - P_atm;
    };

    // We calculate chamber values.
    P_c = Math::AlgebraicSolvers::secant(g, P_atm, 1.2*P_atm, 1000, 1e-6, 1e-6);
    M_c = M_catm * pow(P_c/P_atm, n - 1);
    v_c = M_c * sqrt(gamma*R*T_c);
    rho_c = P_c / (R*T_c);

    double C = M_c / pow(1 + G1*M_c*M_c, G4) * A_c;
    auto f = [=](const double x)
    {
        return x / pow(1 + G1*x*x, G4) * A_e - C;
    };

    // We calculate exit values.
    M_e = Math::AlgebraicSolvers::bisection(f, 0, 1, 1000, 1e-6, 1e-6);
    double mach_factor = (G1*M_c*M_c+ 1) / (G1*M_e*M_e + 1);
    P_e = P_atm;
    T_e = T_c * mach_factor;
    rho_e = P_e / (R*T_e);
    v_e = M_e * sqrt(gamma*R*T_e);
    m_dot = rho_e*v_e*A_e;
    thrust = m_dot*v_e;
}

void SteadySolver::calculate_optimum_parameters()
{
    auto g = [=](double P_c)
    {
        // We calculate M_c.
        double M_c = M_catm * pow(P_c/P_atm, n - 1);
        double P_e = P_c * pow((G1*M_c*M_c+ 1) / G3, G2);
        return P_e - P_atm;
    };

    // We calculate chamber values.
    P_c = Math::AlgebraicSolvers::secant(g, P_atm, 1.2*P_atm, 1000, 1e-6, 1e-6);
    M_c = M_catm * pow(P_c/P_atm, n - 1);
    v_c = M_c * sqrt(gamma*R*T_c);
    rho_c = P_c / (R*T_c);

    // We calculate exit values.
    A_e = A_c * pow(G3, G4) * M_c / pow(1 + G1*M_c*M_c, G4);
    M_e = 1;
    double mach_factor = (G1*M_c*M_c+ 1) / G3;
    P_e = P_atm;
    T_e = T_c * mach_factor;
    rho_e = P_e / (R*T_e);
    v_e = M_e * sqrt(gamma*R*T_e);
    m_dot = rho_e*v_e*A_e;
    thrust = m_dot*v_e;
}

std::string SteadySolver::to_string() const
{
    return std::string("Begin Steady Rocket Solver.\n")
        + "\t— Input:\n"
        + "\t\t• chamber Mach number at atmospheric pressure = " + std::to_string(M_catm) + "\n"
        + "\t\t• adiabatic gas constant = " + std::to_string(gamma) + "\n"
        + "\t\t• gas constant = " + std::to_string(R) + " J/(kg·K)\n"
        + "\t\t• power law exponent n = " + std::to_string(n) + "\n"
        + "\t\t• chamber temperature = " + std::to_string(T_c) + " K\n"
        + "\t— Chamber data:\n"
        + "\t\t• chamber area = " + std::to_string(A_c) + " m²\n"
        + "\t\t• chamber pressure = " + std::to_string(P_c) + " Pa\n"
        + "\t\t• chamber speed = " + std::to_string(v_c) + " m/s\n"
        + "\t\t• chamber Mach number = " + std::to_string(M_c) + "\n"
        + "\t\t• chamber gas density = " + std::to_string(rho_c) + " kg/m³\n"
        + "\t— Exit data\n"
        + "\t\t• exit area = " + std::to_string(A_e) + " m²\n"
        + "\t\t• exit pressure = " + std::to_string(P_e) + " Pa\n"
        + "\t\t• exit temperature = " + std::to_string(T_e) + " K\n"
        + "\t\t• exit speed = " + std::to_string(v_e) + " m/s\n"
        + "\t\t• exit Mach number = " + std::to_string(M_e) + "\n"
        + "\t\t• exit density = " + std::to_string(rho_e) + " kg/m³\n"
        + "\t\t• mass flow = " + std::to_string(m_dot) + " kg/s\n"
        + "\t\t• thrust = " + std::to_string(thrust) + " N\n"
        + "End Steady Rocket Solver.";
}

#endif // ROCKET_SOLVER_CPP
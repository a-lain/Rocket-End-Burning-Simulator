#ifndef EXACT_RIEMANN_SOLVER_CPP
#define EXACT_RIEMANN_SOLVER_CPP

#include "ExactRiemannSolver.hpp"
#include "../../Math/AlgebraicSolvers.hpp"
#include <iostream>
#include <cmath>

using namespace Solvers::Gas;

ExactRiemannSolver::ExactRiemannSolver(const double rho_L, const double v_L, const double P_L,
        const double rho_R, const double v_R, const double P_R, const double gamma,
        const double tol)
{
    // We save every input parameter.
    this->gamma = gamma;
    this->rho_L = rho_L;
    this->v_L = v_L;
    this->P_L = P_L;
    this->rho_R = rho_R;
    this->v_R = v_R;
    this->P_R = P_R;

    // We calculate gamma constants.
    G1 = gamma + 1;
    G2 = gamma - 1;
    G3 = G2 / G1;
    G4 = G2 / (2*gamma);
    G5 = G1 / (2*gamma);
    G6 = 1. / gamma;
    G7 = 2. / G1;
    G8 = 2. / G2;

    // We calculate the constants involved in the estimation of the pressure in the star region.
    A_L = G7/rho_L;
    A_R = G7/rho_R;
    B_L = G3*P_L;
    B_R = G3*P_R;    

    // We calculate the left and right sound speeds.
    c_L = sqrt(gamma * P_L / rho_L);
    c_R = sqrt(gamma * P_R / rho_R);

    
    if (rho_R == 0 && P_R == 0) // If we a vacuum right state.
    {
        vacuum_state = VacuumState::VACUUM_RIGHT_STATE;
        S_starL = v_L + 2*c_L / G2;
    }    
    else if (rho_L == 0 && P_L == 0) // If we have a vacuum left state.
    {
        vacuum_state = VacuumState::VACUUM_LEFT_STATE;
        S_starR = v_R - 2*c_R / G2;
    }    
    else
    {
        // We check for the critical condition for vacuum generation.
        if (v_R - v_L >= 2./G2 * (c_L + c_R)) // We generate vacuum.
        {
            vacuum_state = VacuumState::GENERATED_VACUUM;
            S_starL = v_L + 2*c_L / G2;
            S_starR = v_R - 2*c_R / G2;
        }
        
        else // This is the normal case (no vacuum).
        {
            vacuum_state = VacuumState::NO_VACUUM;
            // We have to find the pressure in the star region.
            // We use the easiest estimate: the average.
            P_star = Math::AlgebraicSolvers::newton_raphson(
                [&] (double P) {return this->f(P);},
                [&] (double P) {return this->dfdP(P);},
                1./2*(P_L + P_R), 1000, 0, tol);
            
            // If we have obtained a pressure that doesn't make sense,
            // we try again with a different initial guess.
            if (std::isnan(P_star))
            {
                P_star =  Math::AlgebraicSolvers::newton_raphson(
                [&] (double P) {return this->f(P);},
                [&] (double P) {return this->dfdP(P);},
                tol, 1000, 0, tol);
            }
            v_star = 1./2*(v_L + v_R + f_R(P_star) - f_L(P_star)); 

            // We deal with the left side first.
            double P_star_by_P_L = P_star / P_L;
            if (P_star > P_L) // left shock
            {                
                rho_starL = rho_L * (P_star_by_P_L + G3) / (G3*P_star_by_P_L + 1);
                S_L = v_L - c_L*sqrt(G5 * P_star_by_P_L + G4);
            }
            else // left rarefaction
            {
                rho_starL = rho_L * pow(P_star_by_P_L, G6);
                S_HL = v_L - c_L;
                S_TL = v_star - c_L*pow(P_star_by_P_L, G4);
            }

            // We proceed with the right side.
            double P_star_by_P_R = P_star / P_R;
            if (P_star > P_R) // right shock
            {
                rho_starR = rho_R * (P_star_by_P_R + G3) / (G3*P_star_by_P_R + 1);
                S_R = v_R + c_R*sqrt(G5 * P_star_by_P_R + G4);
            }
            else // right rarefaction
            {
                rho_starR = rho_R * pow(P_star_by_P_R, G6);
                S_HR = v_R + c_R;
                S_TR = v_star + c_R*pow(P_star_by_P_R, G4);
            }
        }
    }
}

double ExactRiemannSolver::f_L(const double P) const
{
    if (P > P_L) // left shock
    {
        return (P - P_L)*sqrt(A_L / (P + B_L));
    }
    else // left rarefaction
    {
        return 2*c_L/G2*(pow(P/P_L, G4) - 1);
    }
}

double ExactRiemannSolver::f_R(const double P) const
{
    if (P > P_R) // right shock
    {
        return (P - P_R)*sqrt(A_R / (P + B_R));
    }
    else // right rarefaction
    {
        return 2*c_R/G2*(pow(P/P_R, G4) - 1);
    }
}

double ExactRiemannSolver::f(const double P) const
{
    return f_L(P) + f_R(P) + v_R - v_L;
}

double ExactRiemannSolver::df_LdP(const double P) const
{
    if (P > P_L) // left shock
    {
        return sqrt(A_L / (B_L + P)) * (1 - (P - P_L) / (2*(B_L + P)));
    }
    else // left rarefaction
    {
        return 1./rho_L/c_L / pow(P/P_L, G5);
    }
}

double ExactRiemannSolver::df_RdP(const double P) const
{
    if (P > P_R) // right shock
    {
        return sqrt(A_R / (B_R + P)) * (1 - (P - P_R) / (2*(B_R + P)));
    }
    else // left rarefaction
    {
        return 1./rho_R/c_R / pow(P/P_R, G5);
    }
}

double ExactRiemannSolver::dfdP(const double P) const
{
    return df_LdP(P) + df_RdP(P);
}

double ExactRiemannSolver::rho_Lrf(const double S) const
{
    return rho_L * pow(G7 + G3/c_L*(v_L - S), G8);
}

double ExactRiemannSolver::v_Lrf(const double S) const
{
    return G7 * (c_L + v_L/G8 + S);
}

double ExactRiemannSolver::P_Lrf(const double S) const
{
    return P_L * pow(G7 + G3/c_L*(v_L - S), 1./G4);
}

double ExactRiemannSolver::rho_Rrf(const double S) const
{
    return rho_R * pow(G7 - G3/c_R*(v_R - S), G8);
}

double ExactRiemannSolver::v_Rrf(const double S) const
{
    return G7 * (-c_R + v_R/G8 + S);
}

double ExactRiemannSolver::P_Rrf(const double S) const
{
    return P_R * pow(G7 - G3/c_R*(v_R - S), 1./G4);
}

double ExactRiemannSolver::rho_L0(const double S) const
{
    if (S <= v_L - c_L)
    {
        return rho_L;
    }
    else if (S >= S_starL)
    {
        return 0;
    }
    else
    {
        return rho_Lrf(S);
    }
}

double ExactRiemannSolver::v_L0(const double S) const
{
    if (S <= v_L - c_L)
    {
        return v_L;
    }
    else if (S >= S_starL)
    {
        return 0;
    }
    else
    {
        return v_Lrf(S);
    }
}

double ExactRiemannSolver::P_L0(const double S) const
{
    if (S <= v_L - c_L)
    {
        return P_L;
    }
    else if (S >= S_starL)
    {
        return 0;
    }
    else
    {
        return P_Lrf(S);
    }
}

double ExactRiemannSolver::rho_R0(const double S) const
{
    if (S <= S_starR)
    {
        return 0;
    }
    else if (S >= v_R + c_R)
    {
        return rho_R;
    }
    else
    {
        return rho_Rrf(S);
    }
}

double ExactRiemannSolver::v_R0(const double S) const
{
    if (S <= S_starR)
    {
        return 0;
    }
    else if (S >= v_R + c_R)
    {
        return v_R;
    }
    else
    {
        return v_Rrf(S);
    }
}

double ExactRiemannSolver::P_R0(const double S) const
{
    if (S <= S_starR)
    {
        return 0;
    }
    else if (S >= v_R + c_R)
    {
        return P_R;
    }
    else
    {
        return P_Rrf(S);
    }
}

double ExactRiemannSolver::rho(const double x, const double t) const
{
    double res = 0;
    if (t == 0)
    {
        if (x >= 0)
        {
            res = rho_R;
        }
        else
        {
            res = rho_L;
        }
    }
    else
    {
        double S = x / t;
        switch (vacuum_state)
        {
            case VacuumState::VACUUM_RIGHT_STATE:
                res = rho_L0(S);
                break;
            
            case VacuumState::VACUUM_LEFT_STATE:
                res = rho_R0(S);
                break;
            case VacuumState::GENERATED_VACUUM:
                if (S <= S_starL)
                {
                    res = rho_L0(S);
                }
                else if (S >= S_starR)
                {
                    res = rho_R0(S);
                }
                else
                {
                    res = 0;
                }
                break;
            case VacuumState::NO_VACUUM:
                if (S <= v_star) // to the left of the contact discontinuity
                {
                    if (P_star > P_L) // left shock
                    {
                        if (S <= S_L)
                        {
                            res = rho_L;
                        }
                        else
                        {
                            res = rho_starL;
                        }
                    }
                    else // left rarefaction
                    {
                        if (S <= S_HL)
                        {
                            res = rho_L;
                        }
                        else if (S >= S_TL)
                        {
                            res = rho_starL;
                        }
                        else
                        {
                            res = rho_Lrf(S);
                        }
                    }
                }
                else // to the right of the contact discontinuity
                {
                    if (P_star > P_R) // right shock
                    {
                        if (S >= S_R)
                        {
                            res = rho_R;
                        }
                        else
                        {
                            res = rho_starR;
                        }
                    }
                    else // right rarefaction
                    {
                        if (S >= S_HR)
                        {
                            res = rho_R;
                        }
                        else if (S <= S_TR)
                        {
                            res = rho_starR;
                        }
                        else
                        {
                            res = rho_Rrf(S);
                        }
                    }
                }
                break;
        }
    }

    return res;
}

double ExactRiemannSolver::v(const double x, const double t) const
{
    double res = 0;

    if (t == 0)
    {
        if (x >= 0)
        {
            res = v_R;
        }
        else
        {
            res = v_L;
        }
    }
    else
    {
        double S = x / t;
        switch (vacuum_state)
        {
            case VacuumState::VACUUM_RIGHT_STATE:
                res = v_L0(S);
                break;
            
            case VacuumState::VACUUM_LEFT_STATE:
                res = v_R0(S);
                break;
            case VacuumState::GENERATED_VACUUM:
                if (S <= S_starL)
                {
                    res = v_L0(S);
                }
                else if (S >= S_starR)
                {
                    res = v_R0(S);
                }
                else
                {
                    res = 0;
                }
                break;
            case VacuumState::NO_VACUUM:
                if (S <= v_star) // to the left of the contact discontinuity
                {
                    if (P_star > P_L) // left shock
                    {
                        if (S <= S_L)
                        {
                            res = v_L;
                        }
                        else
                        {
                            res = v_star;
                        }
                    }
                    else // left rarefaction
                    {
                        if (S <= S_HL)
                        {
                            res = v_L;
                        }
                        else if (S >= S_TL)
                        {
                            res = v_star;
                        }
                        else
                        {
                            res = v_Lrf(S);
                        }
                    }
                }
                else // to the right of the contact discontinuity
                {
                    if (P_star > P_R) // right shock
                    {
                        if (S >= S_R)
                        {
                            res = v_R;
                        }
                        else
                        {
                            res = v_star;
                        }
                    }
                    else // right rarefaction
                    {
                        if (S >= S_HR)
                        {
                            res = v_R;
                        }
                        else if (S <= S_TR)
                        {
                            res = v_star;
                        }
                        else
                        {
                            res = v_Rrf(S);
                        }
                    }
                }
                break;
        }
    }

    return res;
}

double ExactRiemannSolver::P(const double x, const double t) const
{
    double res = 0;

    if (t == 0)
    {
        if (x >= 0)
        {
            res = P_R;
        }
        else
        {
            res = P_L;
        }
    }
    else
    {
        double S = x / t;
        switch (vacuum_state)
        {
            case VacuumState::VACUUM_RIGHT_STATE:
                res = P_L0(S);
                break;
            
            case VacuumState::VACUUM_LEFT_STATE:
                res = P_R0(S);
                break;
            case VacuumState::GENERATED_VACUUM:
                if (S <= S_starL)
                {
                    res = P_L0(S);
                }
                else if (S >= S_starR)
                {
                    res = P_R0(S);
                }
                else
                {
                    res = 0;
                }
                break;
            case VacuumState::NO_VACUUM:
                if (S <= v_star) // to the left of the contact discontinuity
                {
                    if (P_star > P_L) // left shock
                    {
                        if (S <= S_L)
                        {
                            res = P_L;
                        }
                        else
                        {
                            res = P_star;
                        }
                    }
                    else // left rarefaction
                    {
                        if (S <= S_HL)
                        {
                            res = P_L;
                        }
                        else if (S >= S_TL)
                        {
                            res = P_star;
                        }
                        else
                        {
                            res = P_Lrf(S);
                        }
                    }
                }
                else // to the right of the contact discontinuity
                {
                    if (P_star > P_R) // right shock
                    {
                        if (S >= S_R)
                        {
                            res = P_R;
                        }
                        else
                        {
                            res = P_star;
                        }
                    }
                    else // right rarefaction
                    {
                        if (S >= S_HR)
                        {
                            res = P_R;
                        }
                        else if (S <= S_TR)
                        {
                            res = P_star;
                        }
                        else
                        {
                            res = P_Rrf(S);
                        }
                    }
                }
                break;
        }
    }

    return res;
}

double ExactRiemannSolver::S_max() const
{
    double _S_L = 0;
    double _S_R = 0;
    switch (vacuum_state)
    {
        case VacuumState::VACUUM_RIGHT_STATE:
            _S_L = v_L - c_L;
            _S_R = S_starL;
            break;
        
        case VacuumState::VACUUM_LEFT_STATE:
            _S_L = S_starR;
            _S_R = v_R + c_R;
            break;
        case VacuumState::GENERATED_VACUUM:
            _S_L = S_HL;
            _S_R = S_HR;
            break;
        case VacuumState::NO_VACUUM:
            if (P_star > P_L) // left shock
            {
                _S_L = S_L;
            }
            else // left rarefaction
            {
                _S_L = S_HL;
            }
            if (P_star > P_R) // right shock
            {
                _S_R = S_R;
            }
            else // right rarefaction
            {
                _S_R = S_HR;
            }
            break;
    }
    return fmax(fabs(_S_L), fabs(_S_R));
}

#endif // EXACT_RIEMANN_SOLVER_CPP
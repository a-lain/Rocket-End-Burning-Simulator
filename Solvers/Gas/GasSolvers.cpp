#ifndef GAS_SOLVERS_CPP
#define GAS_SOLVERS_CPP

#include "GasSolvers.hpp"
#include "ExactRiemannSolver.hpp"
#include <cmath>
#include <iostream>

using namespace Math;
using namespace Solvers;

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}


void Gas::exact(const GasCell& L, const GasCell& R, Math::Vector<double>* F, double* S_max)
{
    Gas::ExactRiemannSolver RS(L.rho, L.v, L.P,
                                R.rho, R.v, R.P,
        L.QR->gamma, EXACT_RIEMANN_SOLVER_RELATIVE_TOLERANCE);
    double rho = RS.rho(0, 1);
    double v = RS.v(0, 1);
    double P = RS.P(0, 1);
    double A = L.A;
    *F = Vector<double>({rho*v*A, (rho*v*v + P)*A,
        v*(1./2*rho*v*v + L.QR->gamma/(L.QR->gamma - 1)*P)*A});
    *S_max = RS.S_max();
}

void Gas::HLL(const GasCell& L, const GasCell& R, Math::Vector<double>* F, double* S_max)
{
    double S_L, S_R;
    Vector<double> F_L, F_R;

    // Step I: We do a Roe estimation of the wave speeds.
    double beta_L = sqrt(L.rho*L.A);
    double beta_R = sqrt(R.rho*R.A);
    double v_tilde = (beta_L * L.v + beta_R * R.v) / (beta_L + beta_R);
    beta_L = sqrt(L.rho);
    beta_R = sqrt(R.rho);
    double H_tilde = (beta_L * L.H + beta_R * R.H) / (beta_L + beta_R);
    double c_tilde = sqrt((L.QR->gamma - 1) * (H_tilde - 1./2*v_tilde*v_tilde));
    S_L = v_tilde - c_tilde;
    S_R = v_tilde + c_tilde;

    // Step II: We obtain the flux.
    if (0 <= S_L)
    {
        *F = L.F();
    }
    else if (0 >= S_R)
    {
        *F = R.F();
    }
    else
    {
        F_L = L.F();
        F_R = R.F();
        *F = (S_R*F_L - S_L*F_R + S_L*S_R*(R.U - L.U)) / (S_R - S_L);
    }
    *S_max = fmax(fabs(S_L), fabs(S_R));
}

void Gas::HLLC(const GasCell& L, const GasCell& R, Math::Vector<double>* F, double* S_max)
{
    double S_L, S_R, S_star;

    // Step I: We do a Roe estimation of the wave speeds.
    double beta_L = sqrt(L.rho*L.A);
    double beta_R = sqrt(R.rho*R.A);
    double v_tilde = (beta_L * L.v + beta_R * R.v) / (beta_L + beta_R);
    beta_L = sqrt(L.rho);
    beta_R = sqrt(R.rho);
    double H_tilde = (beta_L * L.H + beta_R * R.H) / (beta_L + beta_R);
    double c_tilde = sqrt((L.QR->gamma - 1) * (H_tilde - 1./2*v_tilde*v_tilde));
    S_L = v_tilde - c_tilde;
    S_R = v_tilde + c_tilde;
    S_star = (R.P*R.A - L.P*L.A + L.rho*L.v*L.A*(S_L - L.v) - R.rho*R.v*R.A*(S_R - R.v))
        / (L.rho*L.A*(S_L - L.v) - R.rho*R.A*(S_R - R.v));

    // Step II: We obtain the flux.
    if (0 <= S_L)
    {
        *F = L.F();
    }
    else if (0 >= S_R)
    {
        *F = R.F();
    }
    else if (0 <= S_star)
    {
        Vector<double> U_starL = L.rho*L.A * (S_L - L.v)/(S_L - S_star)
            * Vector<double>({1, S_star, L.E/L.rho + (S_star - L.v)*(S_star + L.P/(L.rho*(S_L - L.v)))});
        *F = L.F() + S_L*(U_starL - L.U);
    }
    else
    {
        Vector<double> U_starR = R.rho*R.A * (S_R - R.v)/(S_R - S_star)
            * Vector<double>({1, S_star, R.E/R.rho + (S_star - R.v)*(S_star + R.P/(R.rho*(S_R - R.v)))});
        *F =  R.F() + S_R*(U_starR - R.U);
    }
    *S_max = fmax(fabs(S_L), fmax(fabs(S_R), fabs(S_star)));
}

void Gas::Roe(const GasCell& L, const GasCell& R, Math::Vector<double>* F, double* S_max)
{
    // Step I: we compute Roe average values.
    // double rhoA_tilde = sqrt(L.rho*L.A*R.rho*R.A);
    double beta_L = sqrt(L.rho*L.A);
    double beta_R = sqrt(R.rho*R.A);
    double v_tilde = (beta_L * L.v + beta_R * R.v) / (beta_L + beta_R);
    beta_L = sqrt(L.rho);
    beta_R = sqrt(R.rho);
    double H_tilde = (beta_L * L.H + beta_R * R.H) / (beta_L + beta_R);
    double c_tilde = sqrt((L.QR->gamma - 1) * (H_tilde - 1./2*v_tilde*v_tilde));

    // Step II: we compute the eigenvalues.
    double lambda1 = v_tilde - c_tilde;
    double lambda2 = v_tilde;
    double lambda3 = v_tilde + c_tilde;

    // Step III: we compute the eigenvectors.
    Vector<double> K1({1., v_tilde - c_tilde, H_tilde - v_tilde*c_tilde});
    Vector<double> K2({1., v_tilde, 1./2*v_tilde*v_tilde});
    Vector<double> K3({1., v_tilde + c_tilde, H_tilde + v_tilde*c_tilde});

    // Step IV: we compute the wave strengths.
    // Eq 11.110 from Toro.
    double delta_U_1 = R.U[0] - L.U[0];
    double delta_U_2 = R.U[1] - L.U[1];
    double delta_U_3 = R.U[2] - L.U[2];
    double alpha_2 = (L.QR->gamma - 1)/c_tilde/c_tilde *(delta_U_1*(H_tilde - v_tilde*v_tilde) + v_tilde*delta_U_2
        -delta_U_3);
    double alpha_1 = (delta_U_1*(v_tilde + c_tilde) - delta_U_2 - c_tilde*alpha_2) / (2*c_tilde);
    double alpha_3 = delta_U_1 - (alpha_1 + alpha_2);

    // Step V: we calculate the Roe averaged states.
    double rhoA_starL = L.rho*L.A + alpha_1;
    double v_starL = (L.rho*L.A*L.v + alpha_1*(v_tilde - c_tilde)) / (L.rho*L.A + alpha_1);
    double pA_starL = (L.QR->gamma - 1)*(L.E*L.A + alpha_1*(H_tilde - v_tilde*c_tilde) - 1./2*rhoA_starL*v_starL*v_starL);
    double c_starL = sqrt(L.QR->gamma*pA_starL/rhoA_starL);
    double lambda1_L = L.v - L.c;
    double lambda1_R = v_starL - c_starL;
    
    double rhoA_starR = R.rho*R.A - alpha_3;
    double v_starR = (R.rho*R.A*R.v - alpha_3*(v_tilde + c_tilde)) / (R.rho*R.A - alpha_3);
    double pA_starR = (R.QR->gamma - 1)*(R.E*L.A - alpha_3*(H_tilde + v_tilde*c_tilde) - 1./2*rhoA_starR*v_starR*v_starR);
    double c_starR = sqrt(R.QR->gamma*pA_starR/rhoA_starR);
    double lambda3_L = v_starR + c_starR;
    double lambda3_R = R.v + R.c;

    // Step VI: we compute the flux.
    if (lambda1_L < 0 && lambda1_R > 0) // left entropy fix required
    {
        double lambda1_bar = lambda1_L * (lambda1_R - lambda1) / (lambda1_R - lambda1_L);
        *F = L.F() + lambda1_bar*alpha_1*K1;
        *S_max = fmax(fabs(lambda1_bar), fmax(fabs(lambda2), fabs(lambda3)));
    }
    else if (lambda3_L < 0 && lambda3_R > 0) // right entropy fix required
    {
        double lambda3_bar = lambda3_R * (lambda3 - lambda3_L) / (lambda3_R - lambda3_L);
        *F = R.F() - lambda3_bar*alpha_3*K3;
        *S_max = fmax(fabs(lambda1), fmax(fabs(lambda2), fabs(lambda3_bar)));
    }
    else
    {
        *F = 1./2 * (L.F() + R.F() - alpha_1*fabs(lambda1)*K1 - alpha_2*fabs(lambda2)*K2 - alpha_3*fabs(lambda3)*K3);
        *S_max = fmax(fabs(lambda1), fmax(fabs(lambda2), fabs(lambda3)));
    }
}

#endif // GAS_SOLVERS_CPP
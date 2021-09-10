#ifndef REACTION_CPP
#define REACTION_CPP

#include "Reaction.hpp"
#include "../Utilities/FileOperations.hpp"
#include <cmath>

using namespace Chemistry;

SolidGasReaction::SolidGasReaction(const double rho_s, const double k_s, const double cV_s, const double cV_g, const double R_g,
    const double a, const double P_ref, const double n, const double delta_H):
        rho_s(rho_s), k_s(k_s), cV_s(cV_s), alpha(k_s / (rho_s*cV_s)), cV_g(cV_g), R_g(R_g), gamma((R_g + cV_g) / cV_g),
        a(a), P_ref(P_ref), n(n), delta_H(delta_H)
{
    
}

double SolidGasReaction::v_q(const double P) const
{
    return a * pow(P/P_ref, n);
}

SolidGasReaction::SolidGasReaction(FILE* file)
{
    read_from_file(file);
}

void SolidGasReaction::read_from_file(FILE* file)
{
    using ::read_from_file;
    read_from_file(rho_s, file);
    read_from_file(k_s, file);
    read_from_file(cV_s, file);
    read_from_file(alpha, file);
    read_from_file(cV_g, file);
    read_from_file(R_g, file);
    read_from_file(gamma, file);
    read_from_file(a, file);
    read_from_file(P_ref, file);
    read_from_file(n, file);
    read_from_file(delta_H, file);
}

void SolidGasReaction::write_to_file(FILE* file) const
{
    using ::write_to_file;
    write_to_file(rho_s, file);
    write_to_file(k_s, file);
    write_to_file(cV_s, file);
    write_to_file(alpha, file);
    write_to_file(cV_g, file);
    write_to_file(R_g, file);
    write_to_file(gamma, file);
    write_to_file(a, file);
    write_to_file(P_ref, file);
    write_to_file(n, file);
    write_to_file(delta_H, file);
}

std::string SolidGasReaction::to_string() const
{
    std::string text = "Begin SolidGasReaction\n\tρ_s = " + std::to_string(rho_s)
        + "\n\tk_s = " + std::to_string(k_s)
        + "\n\tcV_s = " + std::to_string(cV_s)
        + "\n\tα = " + std::to_string(alpha)
        + "\n\tcV_g = " + std::to_string(cV_g)
        + "\n\tR_g = " + std::to_string(R_g)
        + "\n\tγ = " + std::to_string(gamma)
        + "\n\ta = " + std::to_string(a)
        + "\n\tP_ref = " + std::to_string(P_ref)
        + "\n\tn = " + std::to_string(n)
        + "\n\tΔH = " + std::to_string(delta_H)
        + "\nEnd SolidGasReaction";
    return text;
}

#endif // REACTION_CPP
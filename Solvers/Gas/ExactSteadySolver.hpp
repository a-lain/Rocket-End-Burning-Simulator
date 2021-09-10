#ifndef EXACT_STEADY_SOLVER_HPP
#define EXACT_STEADY_SOLVER_HPP

#include "../../Math/Interpolation.hpp"
#include <functional>
#include <vector>

namespace Solvers
{
    namespace Gas
    {
        /**
         * @brief This enum class is used to determine whether the subsonic
         * or the supersonic solution should be obtained.
         * 
         */
        enum class SolutionType
        {
            SUBSONIC,
            SUPERSONIC,
        };

        class ExactSteadySolver
        {
            public:
            /**
             * @brief Construct a new Exact Steady Solver object
             * 
             * @param rho_0 rest density.
             * @param P_0 rest pressure.
             * @param T_0 rest temperature.
             * @param M_x Mach number at a point x.
             * @param A_x area at a point x.
             * @param gamma gas adiabatic expansion coefficient.
             * @param solution_type whether to compute the subsonic or the supersonic solution.
             */
            ExactSteadySolver(const double rho_0, const double P_0, const double T_0,
                const double M_x, const double A_x, const double gamma,
                const SolutionType solution_type);

            /**
             * @brief Returns Mach number as a function of area.
             * 
             * @param A area.
             * @return double 
             */
            double M(const double A) const;

            /**
             * @brief Returns density as a function of area.
             * 
             * @param A area.
             * @return double 
             */
            double rho(const double A) const;

            /**
             * @brief Returns speed as a function of area.
             * 
             * @param A area.
             * @return double 
             */
            double v(const double A) const;

            /**
             * @brief Returns pressure as a function of area.
             * 
             * @param A area.
             * @return double 
             */
            double P(const double A) const;

            /**
             * @brief Returns temperature as a function of area.
             * 
             * @param A area.
             * @return double 
             */
            double T(const double A) const;

            protected:
            /**
             * @brief Stores whether the subsonic or the supersonic solution is found.
             * 
             */
            SolutionType solution_type;

            /**
             * @brief The constant f(M(x))A(x).
             * 
             */
            double C;

            /**
             * @brief Rest density of the gas.
             * 
             */
            double rho_0;

            /**
             * @brief Rest pressure of the gas.
             * 
             */
            double P_0;

            /**
             * @brief Rest temperature of the gas.
             * 
             */
            double T_0;

            /**
             * @brief Gas constant.
             * 
             */
            double R;

            /**
             * @brief Gas adiabatic expansion coefficient.
             * 
             */
            double gamma;
            
            /**
             * @brief (gamma - 1) / 2
             * 
             */
            double G1;

            /**
             * @brief gamma / (gamma - 1)
             * 
             */
            double G2;

            /**
             * @brief 1 / (gamma - 1)
             * 
             */
            double G3;

            /**
             * @brief (gamma + 1) / (2*(gamma - 1))
             * 
             */
            double G4;
        };
    }
}


#endif // EXACT_STEADY_SOLVER_HPP
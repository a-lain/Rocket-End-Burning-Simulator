#ifndef ROCKET_SOLVER_HPP
#define ROCKET_SOLVER_HPP

#include <string>

namespace Solvers
{
    namespace Rocket
    {
        /**
         * @brief This object can be used to compute chamber and exit values for a specific
         * rocket geometry.
         * 
         */
        class SteadySolver
        {
            public:
            /**
             * @brief Construct a new Steady Solver object
             * 
             * @param T_c chamber temperature, also known as propellant burning temperature.
             * @param M_catm chamber Mach number when the chamber pressure is the atmospheric pressure.
             * @param n Vieille's power law exponent.
             * @param gamma gas adiabatic expansion coefficient.
             * @param R gas constant.
             * @param A_c chamber area.
             */
            SteadySolver(const double T_c, const double M_catm, const double n, 
                const double gamma, const double R, const double A_c);

            /**
             * @brief Determines all chamber and exit values.
             * 
             * @param A_e exit area.
             */
            void solve_for_exit_area(const double A_e);

            /**
             * @brief Calculates the chamber and exit values that provide the maximum thrust,
             * that is, when the exit Mach number is set to one. 
             * 
             */
            void calculate_optimum_parameters();

            /**
             * @brief Chamber area.
             * 
             */
            double A_c;

            /**
             * @brief Exit area.
             * 
             */
            double A_e;

            /**
             * @brief Chamber Mach number when chamber pressure is the atmospheric pressure.
             * 
             */
            double M_catm;

            /**
             * @brief Gas adiabatic expansion coefficient.
             * 
             */
            double gamma;

            /**
             * @brief Gas constant.
             * 
             */
            double R;

            /**
             * @brief Vieille's power law exponent.
             * 
             */
            double n;

            /**
             * @brief Chamber temperature. Also, the burning temperature of the propellant.
             * 
             */
            double T_c;

            /**
             * @brief Chamber gas speed.
             * 
             */
            double v_c;

            /**
             * @brief Chamber Mach number.
             * 
             */
            double M_c;

            /**
             * @brief Chamber pressure.
             * 
             */
            double P_c;

            /**
             * @brief Chamber gas density.
             * 
             */
            double rho_c;

            /**
             * @brief Exit speed.
             * 
             */
            double v_e;

            /**
             * @brief Exit Mach number.
             * 
             */
            double M_e;

            /**
             * @brief Exit pressure.
             * 
             */
            double P_e;

            /**
             * @brief Exit temperature.
             * 
             */
            double T_e;

            /**
             * @brief Exist density.
             * 
             */
            double rho_e;

            /**
             * @brief Mass flux.
             * 
             */
            double m_dot;

            /**
             * @brief Thrust generated by the rocket.
             * 
             */
            double thrust;

            /**
             * @brief Returns a string representation of the object.
             * 
             * @return std::string 
             */
            std::string to_string() const;

            protected:
            /**
             * @brief The atmospheric pressure.
             * 
             */
            constexpr static double P_atm = 101325;

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
             * @brief (gamma + 1) / 2
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

#endif // ROCKET_SOLVER_HPP
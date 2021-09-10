#ifndef EXACT_RIEMANN_SOLVER_HPP
#define EXACT_RIEMANN_SOLVER_HPP

namespace Solvers
{
    namespace Gas
    {
        /**
         * @brief This enum class is used to store whether there is vacuum and, if that is
         * the case, which problem we have to solve.
         * 
         */
        enum class VacuumState
        {
            VACUUM_LEFT_STATE,
            VACUUM_RIGHT_STATE,
            GENERATED_VACUUM,
            NO_VACUUM,
        };

        /**
         * @brief Solves the Riemann problem exactly for the 1D Euler Equations.
         * 
         * Only valid when the area doesn't change. Computes the solution as a
         * function of x and t.
         * 
         * The algorithm works by first computing the pressure in
         * the star region numerically. Afterwards, the speed in the star region is calculated.
         * Next, the character (shock wave or rarefaction wave) of the left and right waves
         * is found out. Lastly, depending on the value of the quotient x/t, the correct
         * zone is selected and the associated formulae are used to calculate the variables.
         * 
         */
        class ExactRiemannSolver
        {
            public:
            /**
             * @brief Gas density as a function of space and time.
             * 
             * @param x spatial coordinate.
             * @param t time coordinate.
             * @return double 
             */
            double rho(const double x, const double t) const;

            /**
             * @brief Gas speed as a function of space and time.
             * 
             * @param x spatial coordinate.
             * @param t time coordinate.
             * @return double 
             */
            double v(const double x, const double t) const;

            /**
             * @brief Gas pressure as a function of space and time.
             * 
             * @param x spatial coordinate.
             * @param t time coordinate.
             * @return double 
             */
            double P(const double x, const double t) const;

            /**
             * @brief Returns the maximum absolute value of all wave speeds of the solution.
             * 
             * @return double 
             */
            double S_max() const;

            /**
             * @brief Construct a new Exact Riemann Solver object.
             * 
             * @param rho_L the gas density to the left of x=0.
             * @param v_L the gas speed to the left of x=0.
             * @param P_L the gas pressure to the left of x=0.
             * @param rho_R the gas density to the right of x=0.
             * @param v_R the gas speed to the right of x=0.
             * @param P_R the gas pressure to the right of x=0.
             * @param gamma the gas adiabatic expansion coefficient.
             * @param tol the precision with which the pressure in the
             * star region is determined.
             */
            ExactRiemannSolver(const double rho_L, const double v_L, const double P_L,
                const double rho_R, const double v_R, const double P_R, const double gamma,
                const double tol);

            protected:

            /**
             * @brief Gas adiabatic expansion coefficient.
             * 
             */
            double gamma;

            /**
             * @brief Stores what type of problem we have to solve.
             * 
             */
            VacuumState vacuum_state;

            /**
             * @brief The gas density to the left of x=0.
             * 
             */
            double rho_L;

            /**
             * @brief The gas speed to the left of x=0.
             * 
             */
            double v_L;

            /**
             * @brief The gas pressure to the left of x=0.
             * 
             */
            double P_L;

            /**
             * @brief The gas density to the right of x=0.
             * 
             */
            double rho_R;

            /**
             * @brief The gas speed to the right of x=0.
             * 
             */
            double v_R;

            /**
             * @brief The gas pressure to the right of x=0.
             * 
             */
            double P_R;            

            // Gamma constants.
            /**
             * @brief gamma + 1
             * 
             */
            double G1;

            /**
             * @brief gamma - 1
             * 
             */
            double G2;

            /**
             * @brief (gamma - 1)/(gamma + 1)
             * 
             */
            double G3;

            /**
             * @brief (gamma - 1)/(2*gamma)
             * 
             */
            double G4;

            /**
             * @brief (gamma + 1)/(2*gamma)
             * 
             */
            double G5;

            /**
             * @brief 1/gamma
             * 
             */
            double G6;

            /**
             * @brief 2 / (gamma + 1)
             * 
             */
            double G7;

            /**
             * @brief 2 / (gamma - 1)
             * 
             */
            double G8;

            /**
             * @brief 2 / (gamma + 1) / rho_L
             * 
             */
            double A_L;

            /**
             * @brief 2 / (gamma + 1) / rho_R
             * 
             */
            double A_R;

            /**
             * @brief (gamma - 1)/(gamma + 1) * P_L
             * 
             */
            double B_L;

            /**
             * @brief (gamma - 1)/(gamma + 1) * P_R
             * 
             */
            double B_R;

            /**
             * @brief The gas sound speed to the left of x=0.
             * 
             */
            double c_L;

            /**
             * @brief The gas sound speed to the right of x=0.
             * 
             */
            double c_R;

            /**
             * @brief Function used in the computation of the pressure equation (left part).
             * 
             * @param P 
             * @return double 
             */
            double f_L(const double P) const;

            /**
             * @brief Function used in the computation of the pressure equation (right part).
             * 
             * @param P 
             * @return double 
             */
            double f_R(const double P) const;

            /**
             * @brief The pressure equation. P_star is the only solution of f(P)=0.
             * 
             * @param P 
             * @return double 
             */
            double f(const double P) const;

            /**
             * @brief The derivative of f_L with respect to P.
             * 
             * @param P 
             * @return double 
             */
            double df_LdP(const double P) const;

            /**
             * @brief The derivative of f_R with respect to P.
             * 
             * @param P 
             * @return double 
             */
            double df_RdP(const double P) const;

            /**
             * @brief The derivative of f with respecto to P.
             * 
             * @param P 
             * @return double 
             */
            double dfdP(const double P) const;

            /**
             * @brief Speed of the left shock wave.
             * 
             */
            double S_L;

            /**
             * @brief Speed of the head of the left rarefaction wave.
             * 
             */
            double S_HL;

            /**
             * @brief Speed of the tail of the left rarefaction wave.
             * 
             */
            double S_TL;

            /**
             * @brief Speed of the vacuum front to the right of x=0.
             * 
             */
            double S_starL;

            /**
             * @brief Density in the star region to the left of the contact discontinuity.
             * 
             */
            double rho_starL;

            /**
             * @brief Speed in the star region.
             * 
             */
            double v_star;

            /**
             * @brief Pressure in the star region.
             * 
             */
            double P_star;

            /**
             * @brief Speed of the right shock wave.
             * 
             */
            double S_R;

            /**
             * @brief Speed of the head of the right rarefaction wave.
             * 
             */
            double S_HR;

            /**
             * @brief Speed of the tail of the right rarefaction wave.
             * 
             */
            double S_TR;

            /**
             * @brief Speed of the vacuum front to the left of x=0.
             * 
             */
            double S_starR;

            /**
             * @brief Density in the star region to the right of the contact discontinuity.
             * 
             */
            double rho_starR;

            /**
             * @brief Returns the density in the left rarefaction region (between the head and the tail).
             * 
             * @param S 
             * @return double 
             */
            double rho_Lrf(const double S) const;

            /**
             * @brief Returns the speed in the left rarefaction region (between the head and the tail).
             * 
             * @param S 
             * @return double 
             */
            double v_Lrf(const double S) const;

            /**
             * @brief Returns the pressure in the left rarefaction region (between the head and the tail).
             * 
             * @param S 
             * @return double 
             */
            double P_Lrf(const double S) const;

            /**
             * @brief Returns the density in the right rarefaction region (between the head and the tail).
             * 
             * @param S 
             * @return double 
             */
            double rho_Rrf(const double S) const;

            /**
             * @brief Returns the speed in the right rarefaction region (between the head and the tail).
             * 
             * @param S 
             * @return double 
             */
            double v_Rrf(const double S) const;

            /**
             * @brief Returns the pressure in the right rarefaction region (between the head and the tail).
             * 
             * @param S 
             * @return double 
             */
            double P_Rrf(const double S) const;

            /**
             * @brief Returns the density for the right vacuum state Riemann problem.
             * 
             * @param S 
             * @return double 
             */
            double rho_L0(const double S) const;

            /**
             * @brief Returns the speed for the right vacuum state Riemann problem.
             * 
             * @param S 
             * @return double 
             */
            double v_L0(const double S) const;

            /**
             * @brief Returns the pressure for the right vacuum state Riemann problem.
             * 
             * @param S 
             * @return double 
             */
            double P_L0(const double S) const;

            /**
             * @brief Returns the density for the left vacuum state Riemann problem.
             * 
             * @param S 
             * @return double 
             */
            double rho_R0(const double S) const;

            /**
             * @brief Returns the speed for the left vacuum state Riemann problem.
             * 
             * @param S 
             * @return double 
             */
            double v_R0(const double S) const;

            /**
             * @brief Returns the pressure for the left vacuum state Riemann problem.
             * 
             * @param S 
             * @return double 
             */
            double P_R0(const double S) const;
        };
    }
}

#endif // EXACT_RIEMANN_SOLVER_HPP
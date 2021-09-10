#ifndef GAS_SOLVERS_HPP
#define GAS_SOLVERS_HPP

#include "../../Math/Vector.hpp"
#include "../../Mesh/Cell.hpp"

namespace Solvers
{
    namespace Gas
    {
        /**
         * @brief Tolerance used when the exact Riemman solver is called.
         * 
         */
        const double EXACT_RIEMANN_SOLVER_RELATIVE_TOLERANCE = 1e-6;

        /**
         * @brief Calculates the flux at the interface between cells L and R
         * using the exact Riemann Solver.
         * 
         * @warning Only valid when area is kept constant.
         * 
         * @param[in] L left cell.
         * @param[in] R right cell.
         * @param[out] F pointer to the numerical flux.
         * @param[out] S_max pointer to the maximum wave speed (in absolute value)
         * that intervenes in the Riemann problem.
         */
        void exact(const GasCell& L, const GasCell& R, Math::Vector<double>* F, double* S_max);

        /**
         * @brief Calculates the flux at the interface between cells L and R
         * using the HLL Solver.
         * 
         * @param[in] L left cell.
         * @param[in] R right cell.
         * @param[out] F pointer to the numerical flux.
         * @param[out] S_max pointer to the maximum wave speed (in absolute value)
         * that intervenes in the Riemann problem. 
         */
        void HLL(const GasCell& L, const GasCell& R, Math::Vector<double>* F, double* S_max);

        /**
         * @brief Calculates the flux at the interface between cells L and R
         * using the HLLC Solver.
         * 
         * @param[in] L left cell.
         * @param[in] R right cell.
         * @param[out] F pointer to the numerical flux.
         * @param[out] S_max pointer to the maximum wave speed (in absolute value)
         * that intervenes in the Riemann problem.
         */
        void HLLC(const GasCell& L, const GasCell& R, Math::Vector<double>* F, double* S_max);

        /**
         * @brief Calculates the flux at the interface between cells L and R
         * using the Roe Solver with Harten-Hyman Entropy Fix.
         * 
         * @param[in] L left cell.
         * @param[in] R right cell.
         * @param[out] F pointer to the numerical flux.
         * @param[out] S_max pointer to the maximum wave speed (in absolute value)
         * that intervenes in the Riemann problem.
         */
        void Roe(const GasCell&L, const GasCell& R, Math::Vector<double>* F, double* S_max);
    }
}

#endif // GAS_SOLVERS_HPP
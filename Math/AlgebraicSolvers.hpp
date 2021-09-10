/** 
 * @file AlgebraicSolvers.hpp
 * @brief This file contains functions that can be used to solve non-linear algebraic equations.
 * @author Andrés Laín Sanclemente
 * @version 0.2.0
 * @date 9th September 2021 
 * 
 */

#ifndef ALGEBRAIC_SOLVERS_HPP
#define ALGEBRAIC_SOLVERS_HPP

#include <functional>
#include "DualNumbers.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

namespace Math
{
    /**
     * @brief Contains solvers for algebraic equations.
     * 
     */
    namespace AlgebraicSolvers
    {
        /**
         * @brief Applies the bisection method to the function f in
         * the interval (a,b).
         * 
         * Returns the current estimation to the solution as soon as
         * one of the following things happens:
         *  -# The maximun number of iterations is reached.
         *  -# The absolute error is smaller then the tolerance given.
         *  -# The relative error is smaller then the tolerance given.
         * 
         * To disable either the absolute or relative error check, set
         * the corresponding tolerance to zero.
         * 
         * @param f the function.
         * @param a the left side of the interval.
         * @param b the right side of the interval.
         * @param iter_max the maximum number of iterations.
         * @param abs_tol the tolerance for the absolute error.
         * @param rel_tol the tolerance for the relative error.
         * @return double 
         */
        double bisection(std::function<double(const double x)> f,
            double a, double b,
            const unsigned int iter_max, const double abs_tol, const double rel_tol);

        /**
         * @brief Applies the secant method to the function f
         * with initial guesses x1 and x2.
         * 
         * Returns the current estimation to the solution as soon as
         * one of the following things happens:
         *  -# The maximun number of iterations is reached.
         *  -# The absolute error is smaller then the tolerance given.
         *  -# The relative error is smaller then the tolerance given.
         * 
         * To disable either the absolute or relative error check, set
         * the corresponding tolerance to zero.
         * 
         * @param f the function.
         * @param x1 the first initial estimation.
         * @param x2 the second initial estimation.
         * @param iter_max the maximum number of iterations.
         * @param abs_tol the tolerance for the absolute error.
         * @param rel_tol the tolerance for the relative error.
         * @return double 
         */
        double secant(std::function<double(const double x)> f,
            double x1, double x2,
            const unsigned int iter_max, const double abs_tol, const double rel_tol);

        /**
         * @brief Applies the Newton-Raphson method to the function f
         * with initial guess x0.
         * 
         * Returns the current estimation to the solution as soon as
         * one of the following things happens:
         *  -# The maximun number of iterations is reached.
         *  -# The absolute error is smaller then the tolerance given.
         *  -# The relative error is smaller then the tolerance given.
         * 
         * To disable either the absolute or relative error check, set
         * the corresponding tolerance to zero.
         * 
         * @param f the function.
         * @param dfdx the derivative of the function.
         * @param x0 the initial estimation.
         * @param iter_max the maximum number of iterations.
         * @param abs_tol the tolerance for the absolute error.
         * @param rel_tol the tolerance for the relative error.
         * @return double 
         */
        double newton_raphson(std::function<double(const double x)>f,
            std::function<double(const double x)>dfdx, const double x0,
            const unsigned int iter_max, const double abs_tol, const double rel_tol);

        /**
         * @brief Applies the Newton-Raphson method to the function f
         * with initial guess x0.
         * 
         * Returns the current estimation to the solution as soon as
         * one of the following things happens:
         *  -# The maximun number of iterations is reached.
         *  -# The absolute error is smaller then the tolerance given.
         *  -# The relative error is smaller then the tolerance given.
         * 
         * To disable either the absolute or relative error check, set
         * the corresponding tolerance to zero.
         * 
         * @param f the function (returning a dual number).
         * @param x0 the initial estimation.
         * @param iter_max the maximum number of iterations.
         * @param abs_tol the tolerance for the absolute error.
         * @param rel_tol the tolerance for the relative error.
         * @return double 
         */
        double newton_raphson(std::function<DualNumber<double>(const DualNumber<double> x)>f,
            const double x0,
            const unsigned int iter_max, const double abs_tol, const double rel_tol);


        /**
         * @brief Solves linear systems of the form Ax=b through LU decomposition with
         * partial pivoting.
         *
         * For efficiency reasons, the coefficient matrix must be supplied first, then
         * the LU decomposition is done and, afterwars, we may solve any system of the
         * form Ax=b without having to recompute the LU decomposition of A.
         * 
         */
        class LinearSystemSolver
        {
            public:
            /**
             * @brief Returns the solution of the system Ax=b.
             * 
             * @param b the independent vector.
             * @return Vector<double> 
             */
            Vector<double> solve(const Vector<double>& b) const;

            /**
             * @brief Construct a new Linear System Solver object
             * for the coefficient matrix A. It computes its LU decomposition
             * with partial pivoting.
             * 
             * @param A 
             */
            explicit LinearSystemSolver(const Matrix<double>& A);

            protected:

            /**
             * @brief The lower diagonal matrix of the LU decomposition of A (PA=LU).
             * 
             */
            Matrix<double> L;
            
            /**
             * @brief The upper diagonal matrix of the LU decomposition of A (PA=LU).
             * 
             */
            Matrix<double> U;

            /**
             * @brief The permutation vector. 
             * 
             * It is used to keep track of the permutations done when applying Gaussian.
             * The i-th component of the permutation vector stores which row of the matrix
             * A is now the i-th row of LU.
             * 
             */
            std::vector<unsigned int> P;
        };
    }
}

#endif // ALGEBRAIC_SOLVERS_HPP
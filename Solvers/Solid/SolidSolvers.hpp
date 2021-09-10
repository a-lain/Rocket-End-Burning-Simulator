#ifndef SOLID_SOLVERS_HPP
#define SOLID_SOLVERS_HPP

#include <functional>
#include "../../Mesh/Mesh.hpp"

/**
 * @brief This enum class is used to classify different possible boundary conditions.
 * 
 */
enum class SolidBoundaryConditionsType
{
    /**
     * @brief Temperature is known as a function of time.
     * 
     */
    FIXED_TEMPERATURE,

    /**
     * @brief The temperature gradient is known as a function of both temperature and time.
     * 
     */
    FIXED_GRADIENT,

    /**
     * @brief Special boundary condition used for the combustion simulation.
     * 
     */
    COMBUSTION_FRONT,
};

/**
 * @brief Object used to store left and right boundary conditions of the solid.
 * 
 */
class SolidBoundaryConditions
{
    public:
    /**
     * @brief Boundary condition at the left side of the mesh.
     * 
     */
    SolidBoundaryConditionsType left;

    /**
     * @brief Boundary condition at the right side of the mesh.
     * 
     */
    SolidBoundaryConditionsType right;

    /**
     * @brief Function of temperature and time that returns the value for the left boundary condition.
     * 
     */
    std::function<double(double T, double t)> left_condition;

    /**
     * @brief Function of temperature and time that provides the value for the right boundary condition.
     * 
     */
    std::function<double(double T, double t)> right_condition;

    /**
     * @brief Construct a new Solid Boundary Conditions object
     * 
     * @param left type of the left boundary condition.
     * @param right type of the right boundary condition.
     * @param left_condition function used to provide values for the left boundary condition.
     * @param right_condition function used to provide values for the right boundary condition.
     */
    SolidBoundaryConditions(const SolidBoundaryConditionsType left = SolidBoundaryConditionsType::FIXED_TEMPERATURE,
        const SolidBoundaryConditionsType right = SolidBoundaryConditionsType::FIXED_TEMPERATURE,
        std::function<double(double T, double t)> left_condition = [](double T, double t){return 0;},
        std::function<double(double T, double t)> right_condition = [](double T, double t){return 0;});
};

namespace Solvers
{
    namespace Solid
    {
        /**
         * @brief Applies Euler explicit method to advance a time step dt. If dt is greater than
         * the maximum stable time step, the time iteration is subcycled.
         * 
         * @param mesh solid mesh.
         * @param dt time step.
         * @param CFL Courant-Friedrichs-Levy number.
         * @param BC boundary conditions.
         * @param t time.
         */
        void euler_explicit(SolidMesh& mesh, double& dt, const double CFL, const SolidBoundaryConditions BC, const double t);

        /**
         * @brief Applies Euler implicit method to advance a time step dt.
         * 
         * @warning Robin boundary conditions only work when the gradient is a linear function of T.
         * 
         * @param mesh solid mesh.
         * @param dt time step.
         * @param CFL Courant-Friedrichs-Levy number.
         * @param BC boundary conditions.
         * @param t time.
         */
        void euler_implicit(SolidMesh& mesh, double& dt, const double CFL, const SolidBoundaryConditions BC, const double t);
    }
}

#endif // SOLID_SOLVERS_HPP
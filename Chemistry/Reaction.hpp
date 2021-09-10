/** 
 * @file Reaction.hpp
 * @brief This file contains objects used to store constants and to compute the
 * kinetics of chemical reactions.
 * @author Andrés Laín Sanclemente
 * @version 0.2.0
 * @date 9th September 2021 
 * 
 */

#ifndef REACTION_HPP
#define REACTION_HPP

#include <cstdio>
#include <string>

/**
 * @brief The objects of this library are used to store constants and to compute
 * the speed of Chemical Reactions.
 * 
 */
namespace Chemistry
{
    /**
     * @brief This object is used to store all properties of a
     * solid reactant and a gas product. It is also used to compute the
     * speed of the combustion front.
     * 
     */
    class SolidGasReaction
    {
        public:

        /**
         * @brief Solid density [kg/m^3].
         * 
         */
        double rho_s;

        /**
         * @brief Thermal conductivity of solid [W/(m·K)].
         * 
         */
        double k_s;

        /**
         * @brief Specific heat capacity at constant volume of solid [J/(kg·K)].
         * 
         */
        double cV_s;

        /**
         * @brief Thermal diffusivity of solid [m^2/s].
         * 
         */
        double alpha;
        
        /**
         * @brief Specific heat capacity at constant volume of gas [J/(kg·K)].
         * 
         */
        double cV_g;

        /**
         * @brief Gas constant of gas [J/(kg·K)].
         * 
         */
        double R_g;

        /**
         * @brief Adiabatic compression coefficent of gas [adimensional].
         * 
         * @return double 
         */
        double gamma;
        
        /**
         * @brief The prepotential factor of Vieille's law [m/s].
         * 
         */
        double a;
        
        /**
         * @brief The reference pressure for Vieille's law [Pa].
         * 
         */
        double P_ref;

        /**
         * @brief The exponent of Vieille's law [adimensional].
         * 
         */
        double n;

        /**
         * @brief The enthalpy increase of the chemical reaction [J/kg].
         * 
         */
        double delta_H;

        /**
         * @brief Computes the burning rate at the given conditions.
         * 
         * @param P the pressure of the gas.
         * @return double 
         */
        double v_q(const double P) const;

        /**
         * @brief Construct a new Solid Gas Reaction object
         * 
         * @param rho_s the solid density [kg/m^3].
         * @param k_s the solid thermal conductivity [W/(m·K)].
         * @param cV_s the solid heat capacity at constant volume [J/(kg·K)].
         * @param cV_g the gas heat capacity at constant volume [J/(kg·K)].
         * @param R_g the gas constant [J/(kg·K)].
         * @param a the prepotential factor of Vieille's law [m/s].
         * @param P_ref the reference pressure for Vieille's law [Pa].
         * @param n the exponent of Vieille's law [adimensional].
         * @param delta_H the enthalpy increase of the chemical reaction [J/kg].
         */
        SolidGasReaction(const double rho_s = 0, const double k_s = 0, const double cV_s = 0,
            const double cV_g = 0, const double R_g = 0, const double a = 0,
            const double P_ref = 0, const double n = 0, const double delta_H = 0);

        /**
         * @brief Construct a new Solid Gas Reaction object reading it from a file.
         * 
         * @param file 
         */
        explicit SolidGasReaction(FILE* file);

        /**
         * @brief Reads the object from file.
         * 
         * @param file 
         */
        void read_from_file(FILE* file);

        /**
         * @brief Writes the object to file.
         * 
         * @param file 
         */
        void write_to_file(FILE* file) const;

        /**
         * @brief Returns a string representation of the object.
         * 
         * @return std::string 
         */
        std::string to_string() const;
    };
}

#endif // REACTION_HPP
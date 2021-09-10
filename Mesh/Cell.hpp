/** 
 * @file Cell.hpp
 * @brief This files contains all declarations related to the basic objects of all meshes: the cells.
 * @author Andrés Laín Sanclemente
 * @version 0.2.0
 * @date 9th September 2021 
 * 
 */

#ifndef CELL_HPP
#define CELL_HPP

#include "../Math/Vector.hpp"
#include "../Chemistry/Reaction.hpp"
#include <functional>
#include <string>


/**
 * @brief The basic structure of the space discretization.
 * 
 */
class BaseCell
{
    public:
    /**
     * @brief Vector of cell conserved variables.
     * 
     */
    Math::Vector<double> U;

    /**
     * @brief Left limit of the cell.
     * 
     */
    double a;

    /**
     * @brief Right limit of the cell.
     * 
     */
    double b;

   /**
    * @brief Area of the cell.
    * 
    */
    double A;

    /**
     * @brief Pointer to the chemical reaction.
     * 
     */
    Chemistry::SolidGasReaction* QR;

    /**
     * @brief Construct a new Base Cell object.
     * 
     * @param U the vector of preserved variables.
     * @param a the left limit of the cell.
     * @param b the right limit of the cell.
     * @param A the area of the cell.
     * @param QR a pointer to the chemical reaction.
     */
    BaseCell(const Math::Vector<double>& U, const double a, const double b, const double A,
        Chemistry::SolidGasReaction* QR);

    /**
     * @brief Construct a new Cell object from file.
     * 
     * @param file 
     */
    BaseCell(FILE* file, Chemistry::SolidGasReaction* QR);

    /**
     * @brief Returns cell average position.
     * 
     * @return double 
     */
    double x() const;

    /**
     * @brief Returns cell length.
     * 
     * @return double 
     */
    double len() const;

    /**
     * @brief Reads the cell values from file.
     * 
     * @param file 
     */
    void read_from_file(FILE* file);

    /**
     * @brief Writes the cell values to file.
     * 
     * @param file 
     */
    void write_to_file(FILE* file) const;

    /**
     * @brief Returns a string representation of the cell.
     * 
     * @return std::string 
     */
    std::string to_string() const;
};


class GasCell: public BaseCell
{
    public:
    // Constructors.
    /**
     * @brief Construct a new Cell object.
     * 
     * @param U the vector of conserved variables.
     * @param a the left limit of the cell.
     * @param b the right limit of the cell.
     * @param A the area of the cell.
     * @param QR a reference to the quemical reaction.
     */
    GasCell(const Math::Vector<double>& U, double a, double b, const double A,
        Chemistry::SolidGasReaction* QR);

    /**
     * @brief Construct a new Cell object from file.
     * 
     * @param file 
     */
    GasCell(FILE* file, Chemistry::SolidGasReaction* QR);

    /**
     * @brief Returns the flux evaluated at the conserved variables of the cell,
     * that is, F(U).
     * 
     */
    Math::Vector<double> F() const;

    /**
     * @brief Right neighbour of the cell.
     * 
     */
    GasCell* right_neighbour;

    /**
     * @brief Left neighbour of the cell.
     * 
     */
    GasCell* left_neighbour;

    /**
     * @brief Gas density.
     * 
     * @return double 
     */
    double rho;

    /**
     * @brief Gas sound speed.
     * 
     * @return double 
     */
    double c;

    /**
     * @brief Gas speed.
     * 
     * @return double 
     */
    double v;

    /**
     * @brief Gas temperature.
     * 
     * @return double 
     */
    double T;

    /**
     * @brief Gas pressure.
     * 
     * @return double 
     */
    double P;

    /**
     * @brief Gas total enthalpy.
     * 
     */
    double H;

    /**
     * @brief Gas total energy.
     * 
     */
    double E;

    /**
     * @brief Mach number. 
     * 
     */
    double M;

    /**
     * @brief Calculates all cell variables (rho, T, P, ...) from the conserved quantities.
     * 
     */
    void update();

    /**
     * @brief Reads all cell values from file and calls update()
     * 
     * @param file 
     */
    void read_from_file(FILE* file);

    /**
     * @brief Returns a string representation of the GasCell.
     * 
     * @return std::string 
     */
    std::string to_string() const;
};

class SolidCell: public BaseCell
{
    public:
    /**
     * @brief Construct a new Solid Cell object.
     * 
     * @param U the vector of preserved variables.
     * @param a the left limit of the cell.
     * @param b the right limit of the cell.
     * @param A the area of the cell.
     * @param QR a pointer to the chemical reaction.
     */
    SolidCell(const Math::Vector<double>& U, double a, double b, const double A,
        Chemistry::SolidGasReaction* QR);

    /**
     * @brief Construct a new Cell object from file.
     * 
     * @param file 
     */
    SolidCell(FILE* file, Chemistry::SolidGasReaction* QR);

    /**
     * @brief Right neighbour of the cell.
     * 
     */
    SolidCell* right_neighbour;

    /**
     * @brief Left neighbour of the cell.
     * 
     */
    SolidCell* left_neighbour;

    /**
     * @brief Solid temperature.
     * 
     */
    double T;

    /**
     * @brief Updates the value of the temperature from the value of U.
     * 
     */
    void update();

    /**
     * @brief Reads the cell values from file and calls update()
     * 
     * @param file 
     */
    void read_from_file(FILE* file);

    /**
     * @brief Returns a string representation of the SolidCell.
     * 
     * @return std::string 
     */
    std::string to_string() const;
};

#endif // CELL_HPP
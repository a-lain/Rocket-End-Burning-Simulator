#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "Mesh/Mesh.hpp"
#include "Chemistry/Reaction.hpp"
#include "Utilities/FileArray.hpp"
#include "Math/Integration.hpp"
#include "CPGF/Plot2d/Graphic.hpp"
#include "Utilities/Progress.hpp"
#include "Solvers/Solid/SolidSolvers.hpp"
#include <functional>
#include <limits>
#include <vector>

/**
 * @brief This enum class is used to represent the different possible
 * boundary conditions.
 * 
 */
enum class GasBoundaryConditionsType
{
    /**
     * @brief Gas speed is set to zero. The gas cannot leave
     * the computational domain. It is a perfect adiabatic wall.
     * 
     */
    WALL,

    /**
     * @brief The gas can freely leave the computational domain.
     * 
     */
    FREE,

    /**
     * @brief The left neighbour of the first cell is the last cell
     * and the right neighbour of the last cell is the first cell.
     * 
     */
    PERIODIC,

    /**
     * @brief The pressure of the gas is fixed to a certain value.
     * This value can be a function of time.
     * 
     */
    FIXED_PRESSURE,

    /**
     * @brief The mass flow and the enthalpy are fixed to certain values.
     * These values can be functions of time.
     * 
     */
    FIXED_FLOW_ENTHALPY,

    /**
     * @brief Density, speed and pressure are fixed. Basically, everthing
     * is known as a function of time.
     * 
     */
    FIXED_DENSITY_SPEED_PRESSURE,

    /**
     * @brief This special boundary condition is used when boundary conditions
     * are handled differently, like at the solid-gas interface.
     * 
     */
    NONE,
};

/**
 * @brief This class is used to store the boundary conditions associated to
 * a simulation.
 * 
 */
class GasBoundaryConditions
{
    public:
    /**
     * @brief The type of boundary conditions to the left of the domain.
     * 
     */
    GasBoundaryConditionsType left;

    /**
     * @brief The type of boundary conditions to the right of the domain.
     * 
     */
    GasBoundaryConditionsType right;

    /**
     * @brief Function of time used to return the value of the first
     * left boundary condition.
     * 
     */
    std::function<double(double t)> left_condition_1;

    /**
     * @brief Function of time used to return the value of the second
     * left boundary condition.
     * 
     */
    std::function<double(double t)> left_condition_2;

    /**
     * @brief Function of time used to return the value of the third
     * left boundary condition.
     * 
     */
    std::function<double(double t)> left_condition_3;

    /**
     * @brief Function of time used to return the value of the first
     * right boundary condition.
     * 
     */
    std::function<double(double t)> right_condition_1;

    /**
     * @brief Function of time used to return the value of the second
     * right boundary condition.
     * 
     */
    std::function<double(double t)> right_condition_2;

    /**
     * @brief Function of time used to return the value of the third
     * right boundary condition.
     * 
     */
    std::function<double(double t)> right_condition_3;

    /**
     * @brief Construct a new Gas Boundary Conditions object.
     * 
     * @param left the type of the left boundary condition.
     * @param right the type of the right boundary condition.
     * @param left_condition_1 function of time used for the first left boundary condition.
     * @param left_condition_2 function of time used for the second left boundary condition.
     * @param left_condition_3 function of time used for the third left boundary condition.
     * @param right_condition_1 function of time used for the first right boundary condition.
     * @param right_condition_2 function of time used for the second right boundary condition.
     * @param right_condition_3 function of time used for the third right boundary condition.
     */
    GasBoundaryConditions(const GasBoundaryConditionsType& left = GasBoundaryConditionsType::WALL,
        const GasBoundaryConditionsType& right = GasBoundaryConditionsType::WALL,
        std::function<double(double t)> left_condition_1 = [](double t){return 0;},
        std::function<double(double t)> left_condition_2 = [](double t){return 0;},
        std::function<double(double t)> left_condition_3 = [](double t){return 0;},
        std::function<double(double t)> right_condition_1 = [](double t){return 0;},
        std::function<double(double t)> right_condition_2 = [](double t){return 0;},
        std::function<double(double t)> right_condition_3 = [](double t){return 0;});
};

/**
 * @brief This enum class is used to know what type of cells must be constructed
 * and how many solvers must be called.
 * 
 */
enum class SimulationType
{
    /**
     * @brief Used when the simulation only contains gas.
     * 
     */
    GAS,

    /**
     * @brief Used when the simulation only contains solid.
     * 
     */
    SOLID,

    /**
     * @brief Used when the simulation contains both. It is assumed that
     * the solid lies to the left of the interface and, consequently, that
     * the gas lies to the right.
     * 
     */
    BOTH,
};

/**
 * @brief This object is used to create some initial conditions and let the
 * domain evolve in time.
 * 
 */
class Simulation
{
    public:

    /**
     * @brief Gas Mesh at the current time.
     * 
     */
    GasMesh instant_gas_mesh;

    /**
     * @brief Solid Mesh at the current time.
     * 
     */
    SolidMesh instant_solid_mesh;

    /**
     * @brief The burning rate at the current time.
     * 
     */
    double instant_v_q;

    /**
     * @brief The position of the combustion front at the current time.
     * 
     */
    double instant_x_q;

    /**
     * @brief Current time.
     * 
     */
    double instant_t;

    /**
     * @brief The value of this variable varies from 0 to the adaptive_refinement_period - 1
     * and is used to determine when to refine the mesh. Its value is incremented
     * once per iteration.
     * 
     */
    unsigned int refine;

    /**
     * @brief The name of the simulation.
     * 
     */
    std::string name;

    /**
     * @brief Stores which type of simulation we have to solve.
     * 
     */
    SimulationType simulation_type;

    /**
     * @brief An array of gas meshes for each time instant t. It is
     * stored in a file.
     * 
     */
    Utilities::FileArray<GasMesh>* gas_mesh;   

    /**
     * @brief An array of solid meshes for each time instant t. It is
     * stored in a file.
     * 
     */
    Utilities::FileArray<SolidMesh>* solid_mesh;

    /**
     * @brief An array which contains the position of the combustion front
     * for each time t.
     * 
     */
    std::vector<double> x_q_array;

    /**
     * @brief An array which contains the burning rate for each time t.
     * 
     */
    std::vector<double> v_q_array;
     
    /**
     * @brief An array of all computed times.
     * 
     */
    std::vector<double> t_array;

    /**
     * @brief Gas boundary conditions.
     * 
     */
    GasBoundaryConditions gas_BC;

    /**
     * @brief Solid boundary conditions.
     * 
     */
    SolidBoundaryConditions solid_BC;

    /**
     * @brief Courant-Friedrichs-Lewy number.
     * 
     */
    double CFL;
    
    /**
     * @brief Determines how many times during the simulation
     * the current meshes and values are saved to the file.
     * 
     * The number of actual saves is N_saves + 2, because the
     * first and last state are saved by default. The saves are
     * evenly spread throghout the simulation time.
     * 
     */
    unsigned int N_saves;

    /**
     * @brief Determines whether to use adaptive refinement for the solid.
     * 
     */
    bool adaptive_refinement_solid;

    /**
     * @brief Determines whether to use adaptive refinement for the gas.
     * 
     */
    bool adaptive_refinement_gas;

    /**
     * @brief Controls how often the mesh is refined. In particular, it
     * represents how many iterations we wait after the mesh has been refined
     * to refine it again.
     * 
     * Setting this parameter to one means the mesh is refined in every iteration.
     * 
     */
    unsigned int adaptive_refinement_period;

    /**
     * @brief How many (CPU) tasks to use in other to execute the simulation.
     * 
     */
    unsigned int N_tasks;
    
    /**
     * @brief Function that computes the external forces acting on each gas cell.
     * 
     * @return std::function<double(const Cell& C, const double t)> 
     */
    std::function<double(const GasCell& C, const double t)> external_forces;

    /**
     * @brief Function that returns the numerical flux at the interface between gas cells A and B.
     * 
     */
    std::function<void(const GasCell& A, const GasCell& B, Math::Vector<double>* F, double* S_max)> convection_solver;

    /**
     * @brief Function that solves the solid part.
     * 
     */
    std::function<void(SolidMesh& mesh, double& dt, const double CFL, const SolidBoundaryConditions BC, const double t)> diffusion_solver;

    /**
     * @brief Variable used to store the progress of the simulation.
     * 
     */
    Utilities::Progress<double> progress;

    /**
     * @brief Constructs a Simulation object from initial conditions. It is assumed that the
     * simulation is of type BOTH.
     * 
     * @param name the name of the simulation.
     * @param QR the quemical reaction.
     * @param a the left limit of the mesh.
     * @param b the right limit of the mesh.
     * @param N_cells_solid the initial number of cells for the solid mesh.
     * @param N_cells_gas the initial number of cells for the gas mesh.
     * @param adaptive_refinement_solid whether to use adaptive refinement for the solid mesh.
     * @param adaptive_refinement_gas whether to use adaptive refinement for the gas mesh.
     * @param N_tasks the number of CPU tasks to use in order to compute the simulation.
     * @param CFL the Courant-Friedrichs-Lewy number.
     * @param N_saves how many times throughout the simulation is data stored in a file.
     * @param solid_BC solid boundary conditions for the simulation.
     * @param gas_BC gas boundary conditions for the simulation.
     * @param x_q the initial position of the solid-gas interface.
     * @param v the initial gas speed as a function of space.
     * @param P the initial gas pressure as a function of space
     * @param T the initial temperature as a function of space.
     * @param A the conduct area as a function of space.
     * @param convection_solver the solver to use in order to calculate intercell fluxes.
     * @param diffusion_solver the solver to use in order to solve the solid.
     * @param integrator the function to use in order to integrate the initial conditions.
     * @param external_forces a function that computes the external forces acting on a gas cell
     * for a specific moment in time.
     * @param adaptive_refinement_period determines how often the mesh is refined.
     * @param detail_subdivide_threshold if a cell detail is bigger than this number, we divide the
     * cell in two.
     * @param detail_merge_threshold if the details of two neighbouring cells are lower than this
     * number, they are merged.
     * @param max_length_factor a cell cannot be bigger than max_length_factor * total length of the mesh.
     * @param min_length_factor a cell cannot be small than min_length_factor * total length of the mesh.
     * @param boundary_cell_max_length_factor max_length_factor for cells near the boundary.
     * 
     */
    Simulation
    (
        const std::string& name,
        const Chemistry::SolidGasReaction& QR,
        const double a,
        const double b,
        const unsigned int N_cells_solid,
        const unsigned int N_cells_gas,
        const bool adaptive_refinement_solid,
        const bool adaptive_refinement_gas,
        const unsigned int N_tasks,
        const double CFL,
        const unsigned int N_saves,
        const SolidBoundaryConditions& solid_BC,
        const GasBoundaryConditions& gas_BC,
        const double x_q,
        std::function<double(double x)> v,
        std::function<double(double x)> P,
        std::function<double(double x)> T,
        std::function<double(double x)> A,
        std::function<void(const GasCell& A, const GasCell& B, Math::Vector<double>* F, double* S_max)> convection_solver,
        std::function<void(SolidMesh& mesh, double& dt, const double CFL,
            const SolidBoundaryConditions BC, const double t)> diffusion_solver,
        std::function<double(std::function<double(double x)> f, const double a, const double b)> integrator =
            [] (std::function<double(double x)> f, const double a, const double b)
                {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
        std::function<double(const GasCell& C, const double t)> external_forces = [] (const GasCell& C, const double t){return 0.;},
        const unsigned int adaptive_refinement_period = 1,
        const double detail_subdivide_threshold = 0.01,
        const double detail_merge_threshold = 0.001,
        const double max_length_factor = 1./50,
        const double min_length_factor = 1./2000,
        const double boundary_cell_max_length_factor = 1./1000
    );

    /**
     * @brief Constructs a Simulation object from initial conditions. It is assumed that the
     * simulation is of type GAS.
     * 
     * @param name the name of the simulation.
     * @param QR the quemical reaction object that stores the gas constants.
     * @param a the left limit of the mesh.
     * @param b the right limit of the mesh.
     * @param N_cells the initial number of cells.
     * @param adaptive_refinement whether to use adaptive refinement of the mesh.
     * @param N_tasks the number of CPU tasks to use in order to compute the simulation.
     * @param CFL the Courant-Friedrichs-Lewy number.
     * @param N_saves how many times throughout the simulation is data stored in a file.
     * @param BC boundary conditions for the simulation.
     * @param rho the initial gas density as a function of space.
     * @param v the initial gas speed as a function of space.
     * @param P the initial pressure as a function of space.
     * @param A the conduct area as a function of space.
     * @param convection_solver the solver to use in order to solve convection.
     * @param integrator the function to use in order to integrate the initial conditions.
     * @param external_forces a function that computes the external forces acting on a cell
     * for a specific moment in time.
     * @param adaptive_refinement_period determines how often the mesh is refined.
     * @param detail_subdivide_threshold if a cell detail is bigger than this number, we divide the
     * cell in two.
     * @param detail_merge_threshold if the details of two neighbouring cells are lower than this
     * number, they are merged.
     * @param max_length_factor a cell cannot be bigger than max_length_factor * total length of the mesh.
     * @param min_length_factor a cell cannot be small than min_length_factor * total length of the mesh.
     * @param boundary_cell_max_length_factor max_length_factor for cells near the boundary.
     * 
     */
    Simulation
    (
        const std::string& name,
        const Chemistry::SolidGasReaction& QR,
        const double a,
        const double b,
        const unsigned int N_cells,
        const bool adaptive_refinement,
        const unsigned int N_tasks,
        const double CFL,
        const unsigned int N_saves,
        const GasBoundaryConditions& BC,
        std::function<double(double x)> rho,
        std::function<double(double x)> v,
        std::function<double(double x)> P,
        std::function<double(double x)> A,
        std::function<void(const GasCell& A, const GasCell& B, Math::Vector<double>* F, double* S_max)> convection_solver,
        std::function<double(std::function<double(double x)> f, const double a, const double b)> integrator =
            [] (std::function<double(double x)> f, const double a, const double b)
                {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
        std::function<double(const GasCell& C, const double t)> external_forces = [] (const GasCell& C, const double t){return 0.;},
        const unsigned int adaptive_refinement_period = 1,
        const double detail_subdivide_threshold = 0.01,
        const double detail_merge_threshold = 0.001,
        const double max_length_factor = 1./50,
        const double min_length_factor = 1./2000,
        const double boundary_cell_max_length_factor = 1./1000
    );

    /**
     * @brief Constructs a Simulation object from initial conditions. It is assumed that the
     * simulation is of type SOLID.
     * 
     * @param name the name of the simulation.
     * @param QR the quemical reaction.
     * @param a the left limit of the mesh.
     * @param b the right limit of the mesh.
     * @param N_cells the initial number of cells.
     * @param adaptive_refinement whether to use adaptive refinement of the mesh.
     * @param N_tasks the number of CPU tasks to use in order to compute the simulation.
     * @param CFL the Courant-Friedrichs-Lewy number.
     * @param N_saves how many times throughout the simulation is data stored in a file.
     * @param BC boundary conditions for the simulation.
     * @param T the initial temperature as a function of space.
     * @param A the conduct area as a function of space.
     * @param diffusion_solver the solver to use in order to solve diffusion.
     * @param integrator the function to use in order to integrate the initial conditions.
     * @param adaptive_refinement_period determines how often the mesh is refined.
     * @param detail_subdivide_threshold if a cell detail is bigger than this number, we divide the
     * cell in two.
     * @param detail_merge_threshold if the details of two neighbouring cells are lower than this
     * number, they are merged.
     * @param max_length_factor a cell cannot be bigger than max_length_factor * total length of the mesh.
     * @param min_length_factor a cell cannot be small than min_length_factor * total length of the mesh.
     * @param boundary_cell_max_length_factor max_length_factor for cells near the boundary.
     * 
     */
    Simulation
    (
        const std::string& name,
        const Chemistry::SolidGasReaction& QR,
        const double a,
        const double b,
        const unsigned int N_cells,
        const bool adaptive_refinement,
        const unsigned int N_tasks,
        const double CFL,
        const unsigned int N_saves,
        const SolidBoundaryConditions& BC,
        std::function<double(double x)> T,
        std::function<double(double x)> A,
        std::function<void(SolidMesh& mesh, double& dt, const double CFL,
            const SolidBoundaryConditions BC, const double t)> diffusion_solver,
        std::function<double(std::function<double(double x)> f, const double a, const double b)> integrator =
            [] (std::function<double(double x)> f, const double a, const double b)
                {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
        const unsigned int adaptive_refinement_period = 1,
        const double detail_subdivide_threshold = 0.01,
        const double detail_merge_threshold = 0.001,
        const double max_length_factor = 1./50,
        const double min_length_factor = 1./2000,
        const double boundary_cell_max_length_factor = 1./1000
    );

    /**
     * @brief Construct a new Simulation object from a file.
     * 
     * @param file the name of the file where the simulation is stored.
     */
    explicit Simulation(const std::string& file);

    ~Simulation();

    /**
     * @brief Compute a new time step.
     * 
     * @param dt the maximum time step allowed.
     */
    void update(double dt = std::numeric_limits<double>::max());

    /**
     * @brief Runs the simulation until time t has been reached.
     * The last computed time is exactly t.
     * 
     * @param t 
     */
    void simulate_until(const double t);

    /**
     * @brief Writes the simulation into a file.
     * 
     * @param file the name of the file to write the simulation to.
     */
    void write_to_file(const std::string& file) const;

    /**
     * @brief Returns the position of the combustion at a specific
     * time t.
     * 
     * It does a linear interpolation using t_array and x_q_array.
     * 
     * @param t 
     * @return double 
     */
    double x_q(const double t) const;

    /**
     * @brief Returns the burning rate at a specific time t
     * 
     * It does a linear interpolation using t_array and v_q_array.
     * 
     * @param t 
     * @return double 
     */
    double v_q(const double t) const;

    /**
     * @brief Returns the function A(x) for a certain time t. 
     * 
     * A temporal and spatial linear interpolation is used.
     * 
     * @param t 
     * @return std::function<double(const double x)> 
     */
    std::function<double(const double x)> A(const double t) const;

    /**
     * @brief Returns the function rho(x) for a certain time t.
     * 
     * A temporal and spatial linear interpolation is used.
     * 
     * @param t 
     * @return std::function<double(const double x)> 
     */
    std::function<double(const double x)> rho(const double t) const;

    /**
     * @brief Returns the function c(x) for a certain time t.
     * 
     * A temporal and spatial linear interpolation is used.
     * 
     * @param t 
     * @return std::function<double(const double x)> 
     */
    std::function<double(const double x)> c(const double t) const;

    /**
     * @brief Returns the function v(x) for a certain time t.
     * 
     * A temporal and spatial linear interpolation is used.
     * 
     * @param t 
     * @return std::function<double(const double x)> 
     */
    std::function<double(const double x)> v(const double t) const;

    /**
     * @brief Returns the function T(x) for a certain time t.
     * 
     * A temporal and spatial linear interpolation is used.
     * 
     * @param t 
     * @return std::function<double(const double x)> 
     */
    std::function<double(const double x)> T(const double t) const;

    /**
     * @brief Returns the function P(x) for a certain time t.
     * 
     * A temporal and spatial linear interpolation is used.
     * 
     * @param t 
     * @return std::function<double(const double x)> 
     */
    std::function<double(const double x)> P(const double t) const;

    /**
     * @brief Returns the function M(x) for a certain time t.
     * 
     * A temporal and spatial linear interpolation is used.
     * 
     * @param t 
     * @return std::function<double(const double x)> 
     */
    std::function<double(const double x)> M(const double t) const;

    /**
     * @brief Returns a graphic where the cell boundaries are drawn
     * as a function of time.
     * 
     * @return CPGF::Plot2d::Graphic 
     */
    CPGF::Plot2d::Graphic* mesh_plot() const;

    protected:
    /**
     * @brief Coefficient that accompanies the term T_{n-1} in the energy
     * equation of the last solid cell.
     * 
     */
    double a;

    /**
     * @brief Coefficient that accompanies the term T_n in the energy equation
     * of the last solid cell.
     * 
     */
    double b;

    /**
     * @brief Independent term of the energy equation of the last solid cell.
     * 
     */
    double d;

    /**
     * @brief Density of the first gas cell.
     * 
     */
    double rho_g;

    /**
     * @brief Gas speed of the first gas cell.
     * 
     */
    double v_g;
};

#endif // SIMULATION_HPP
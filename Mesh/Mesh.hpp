/** 
 * @file Mesh.hpp
 * @brief This files contains all declarations related to the mesh.
 * @author Andrés Laín Sanclemente
 * @version 0.2.0
 * @date 9th September 2021 
 * 
 */

#ifndef MESH_HPP
#define MESH_HPP

#include "Cell.hpp"
#include <functional>
#include <vector>

/**
 * @brief This object represents a collection of cells. The number
 * of cells used can be changed. Methods for adaptive refinement
 * are included.
 * 
 */
template <class Cell>
class Mesh
{
    public:
    /**
     * @brief Pointer to the first cell of the mesh.
     * 
     */
    Cell* first_cell;

    /**
     * @brief Pointer to the last cell of the mesh.
     * 
     */
    Cell* last_cell;

    /**
     * @brief Chemical reaction object.
     * 
     */
    Chemistry::SolidGasReaction QR;

    /**
     * @brief Returns total number of cells.
     * 
     * @return unsigned int 
     */
    unsigned int N_cells() const;

    /**
     * @brief If a cell detail is bigger than this number, we divide the
     * cell in two.
     * 
     */
    double detail_subdivide_threshold;

    /**
     * @brief If the details of two neighbouring cells are lower than this
     * number, they are merged.
     * 
     */
    double detail_merge_threshold;

    /**
     * @brief Determines the maximum size of a cell. A cell cannot be
     * greater than max_length_factor times the length of the computational domain.
     * 
     */
    double max_length_factor;

    /**
     * @brief Determines the minimum size of a cell. A cell cannot be
     * smaller than min_length_factor times the length of the computational domain.
     * 
     */
    double min_length_factor;

    /**
     * @brief The max length factor for cells closer than 4 cells to the boundary.
     * 
     */
    double boundary_cell_max_length_factor;

    /**
     * @brief Function used to compute the area of the newly created cells when
     * using adaptive refinement.
     * 
     */
    std::function<double(const double x)> A_func;

    /**
     * @brief Iterator object to loop through all the cells of the mesh.
     * 
     */
    class Iterator
    {
        public:
        /**
         * @brief Pointer to the current cell.
         * 
         */
        Cell* cell;

        /**
         * @brief Construct a new Iterator object that points to nullptr.
         * 
         */
        Iterator();

        /**
         * @brief Construct a new Iterator object that points to @a cell.
         * 
         * @param cell 
         */
        Iterator(Cell* cell);
        
        /**
         * @brief Copy constructor.
         * 
         * @param J 
         */
        Iterator(const Iterator& J);

        /**
         * @brief Assignment operator.
         * 
         * @param J 
         * @return Iterator& 
         */
        Iterator& operator=(const Iterator& J);

        /**
         * @brief Returns the cell placed i positions to the right.
         * 
         * @param i 
         * @return Iterator& 
         */
        Iterator operator+(const int i);

        /**
         * @brief Returns the cell placed i positions to the left.
         * 
         * @param i 
         * @return Iterator& 
         */
        Iterator operator-(const int i);

        /**
         * @brief Advance one cell to the right in the mesh.
         * 
         * @return Iterator& 
         */
        Iterator operator++(int);

        /**
         * @brief Advance one cell to the right in the mesh.
         * 
         * @return Iterator& 
         */
        Iterator& operator++();

        /**
         * @brief Advance one cell to the left in the mesh.
         * 
         * @return Iterator& 
         */
        Iterator operator--(int);

        /**
         * @brief Advance one cell to the left in the mesh.
         * 
         * @return Iterator& 
         */
        Iterator& operator--();

        /**
         * @brief Dereference operator returns the cell it poins to.
         * 
         * @return Cell& 
         */
        Cell& operator*();

        /**
         * @brief Dereference operator returns the cell it poins to.
         * 
         * @return Cell* 
         */
        Cell* operator->();

        /**
         * @brief Implicit convertion into cell pointer.
         * 
         * @return Cell* 
         */
        operator Cell*() const;

        /**
         * @brief Returns true if both iterators point to the same cell and
         * false otherwise.
         * 
         * @param I 
         * @param J 
         * @return true 
         * @return false 
         */
        bool operator==(const Iterator& J);

        /**
         * @brief Returns false if both iterators point to the same cell and
         * true otherwise.
         * 
         * @param I 
         * @param J 
         * @return true 
         * @return false 
         */
        bool operator!=(const Iterator& J);
    };

    // Constructors.
    
    /**
     * @brief Construct a new Mesh object
     * 
     * @param first_cell pointer to the first cell.
     * @param last_cell pointer to the last cell.
     * @param QR chemical reaction.
     * @param A_func area function.
     * @param detail_subdivide_threshold determines when cells are subdivided.
     * @param detail_merge_threshold determines when cells are merged.
     * @param max_length_factor limits the maximum length of a created cell.
     * @param min_length_factor limits the minimum length of a created cell.
     * @param boundary_cell_max_length_factor limits the maxium length of cells near the boundary.
     */
    Mesh(Cell* first_cell = nullptr, Cell* last_cell = nullptr,
        const Chemistry::SolidGasReaction& QR = Chemistry::SolidGasReaction(),
        std::function<double(const double x)> A_func = [](double x){return 0;},
        const double detail_subdivide_threshold = 0.01,
        const double detail_merge_threshold = 0.001,
        const double max_length_factor = 1./50,
        const double min_length_factor = 1./2000,
        const double boundary_cell_max_length_factor = 1./1000);

    /**
     * @brief Perfomrs a deep copy of @a mesh .
     * 
     * @param mesh 
     */
    Mesh(const Mesh& mesh);

    /**
     * @brief Move constructor.
     * 
     * @param mesh 
     */
    Mesh(Mesh&& mesh);

    /**
     * @brief Assignment operator.
     * 
     * @param mesh 
     * @return Mesh& 
     */
    Mesh& operator=(const Mesh& mesh);

    /**
     * @brief Destroy the Mesh object.
     * 
     */
    ~Mesh();

    /**
     * @brief Frees memory associated to the mesh by deleting all cells.
     * 
     */
    void free();

    /**
     * @brief Returns an iterator to the first cell.
     * 
     * @return Mesh::Iterator 
     */
    Mesh::Iterator begin() const;

    /**
     * @brief Returns an iterator to the last cell.
     * 
     * @return Mesh::Iterator 
     */
    Mesh::Iterator rbegin() const;

    /**
     * @brief Returns an iterator to a virtual cell after the last one
     * (or before the first one).
     * 
     * @return Mesh::Iterator 
     */
    Mesh::Iterator end() const;

    /**
     * @brief Subdivides cell C. Returns a pointer to the new right cell.
     * 
     * @param C 
     */
    Cell* subdivide_at(Cell* C);

    /**
     * @brief Merges two cells into one. Returns a pointer to the new cell.
     * 
     * @param C1 
     * @param C2 
     * @return Mesh& 
     */
    Cell* merge_cells(Cell* L, Cell* R);

    /**
     * @brief Calculates the maximum and minimum value of the norm of the
     * vector of conserved quatities throughout the mesh. Then, substacts
     * the maximum value and the minimum value, saving the result
     * in the member ranges.
     * 
     */
    void calculate_variable_ranges();

    /**
     * @brief This is function to estimate how different are the values stored
     * in two neighbouring cells.
     * 
     * @param L 
     * @param R 
     * @return double 
     */
    double detail(Cell* L, Cell* R) const;

    /**
     * @brief Loops through the mesh and subdivides and merges
     * cells according to the thresholds.
     * 
     * @return Mesh& 
     */
    void optimize_mesh();

    /**
     * @brief Returns the average position of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> x() const;

    /**
     * @brief Returns the partition of the X axis given by the
     * boundaries of the cells.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> x_partition() const;
    
    /**
     * @brief Returns the area of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> A() const;

    /**
     * @brief Returns the vector of conserved variables of all cells in the mesh.
     * 
     * @return std::vector<Math::Vector<double>> 
     */
    std::vector<Math::Vector<double>> U() const;

    /**
     * @brief Writes the mesh to file.
     * 
     * @param file 
     */
    void read_from_file(FILE* file);

    /**
     * @brief Reads the mesh from file.
     * 
     * @param file 
     */
    void write_to_file(FILE* file) const;

    protected:
    /**
     * @brief The difference between the maximum and the minimum
     * norm of the vector of preserved quantities is stored here.
     * 
     */
    double ranges;
};

class GasMesh: public Mesh<GasCell>
{
    public:
    /**
     * @brief Construct a new Gas Mesh object
     * 
     * @param first_cell pointer to the first cell.
     * @param last_cell pointer to the last cell.
     * @param QR chemical reaction.
     * @param A_func area function.
     * @param detail_subdivide_threshold determines when cells are subdivided.
     * @param detail_merge_threshold determines when cells are merged.
     * @param max_length_factor limits the maximum length of a created cell.
     * @param min_length_factor limits the minimum length of a created cell.
     * @param boundary_cell_max_length_factor limits the maxium length of cells near the boundary.
     */
    GasMesh(GasCell* first_cell = nullptr, GasCell* last_cell = nullptr,
        const Chemistry::SolidGasReaction& QR = Chemistry::SolidGasReaction(),
        std::function<double(const double x)> A_func = [](double x){return 0;},
        const double detail_subdivide_threshold = 0.01,
        const double detail_merge_threshold = 0.0005,
        const double max_length_factor = 1./20,
        const double min_length_factor = 1./10000,
        const double boundary_cell_max_length_factor = 1./1000);

    /**
     * @brief Returns gas sound speed of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> c() const;

    /**
     * @brief Returns the density of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> rho() const;

    /**
     * @brief Returns gas speed of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> v() const;

    /**
     * @brief Returns the temperature of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> T() const;
    
    /**
     * @brief Returns the pressure of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> P() const;

    /**
     * @brief Returns the Mach number of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> M() const;
};

class SolidMesh: public Mesh<SolidCell>
{
    public:
    /**
     * @brief Construct a new Solid Mesh object
     * 
     * @param first_cell pointer to the first cell.
     * @param last_cell pointer to the last cell.
     * @param QR chemical reaction.
     * @param A_func area function.
     * @param detail_subdivide_threshold determines when cells are subdivided.
     * @param detail_merge_threshold determines when cells are merged.
     * @param max_length_factor limits the maximum length of a created cell.
     * @param min_length_factor limits the minimum length of a created cell.
     * @param boundary_cell_max_length_factor limits the maxium length of cells near the boundary.
     */
    SolidMesh(SolidCell* first_cell = nullptr, SolidCell* last_cell = nullptr,
        const Chemistry::SolidGasReaction& QR = Chemistry::SolidGasReaction(),
        std::function<double(const double x)> A_func = [](double x){return 0;},
        const double detail_subdivide_threshold = 0.01,
        const double detail_merge_threshold = 0.0005,
        const double max_length_factor = 1./20,
        const double min_length_factor = 1./2000,
        const double boundary_cell_max_length_factor = 1./1000);

    /**
     * @brief Returns the temperature of all cells in the mesh.
     * 
     * @return std::vector<double> 
     */
    std::vector<double> T() const;
};


#endif // MESH_HPP
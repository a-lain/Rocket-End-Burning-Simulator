#ifndef MESH_CPP
#define MESH_CPP

#include "Mesh.hpp"
#include "../Utilities/FileOperations.hpp"
#include "../Utilities/FileArray.hpp"
#include "../Utilities/FileArray.cpp"
#include "../Math/Integration.hpp"
#include <future>
#include <cmath>

using namespace Math;

template <class Cell>
unsigned int Mesh<Cell>::N_cells() const
{
    unsigned int i = 0;
    for(Mesh::Iterator C = begin(); C != end(); ++C)
    {
        i++;
    }
    return i;
}

template <class Cell>
Mesh<Cell>::Iterator::Iterator(): cell(nullptr)
{

}

template <class Cell>
Mesh<Cell>::Iterator::Iterator(Cell* cell):cell(cell)
{

}

template <class Cell>
Mesh<Cell>::Iterator::Iterator(const Mesh<Cell>::Iterator& J): cell(J.cell)
{

}

template <class Cell>
typename Mesh<Cell>::Iterator& Mesh<Cell>::Iterator::operator=(const Mesh<Cell>::Iterator& J)
{
    cell = J.cell;
    return *this;
}

template <class Cell>
typename Mesh<Cell>::Iterator Mesh<Cell>::begin() const
{
    return Mesh::Iterator(first_cell);
}

template <class Cell>
typename Mesh<Cell>::Iterator Mesh<Cell>::rbegin() const
{
    return Mesh::Iterator(last_cell);
}

template <class Cell>
typename Mesh<Cell>::Iterator Mesh<Cell>::end() const
{
    return Mesh::Iterator(nullptr);
}

template <class Cell>
typename Mesh<Cell>::Iterator& Mesh<Cell>::Iterator::operator++()
{
    cell = cell->right_neighbour;
    return *this;
}

template <class Cell>
typename Mesh<Cell>::Iterator Mesh<Cell>::Iterator::operator++(int)
{
    return Iterator(cell->right_neighbour);
}

template <class Cell>
typename Mesh<Cell>::Iterator& Mesh<Cell>::Iterator::operator--()
{
    cell = cell->left_neighbour;
    return *this;
}

template <class Cell>
typename Mesh<Cell>::Iterator Mesh<Cell>::Iterator::operator--(int)
{
    return Iterator(cell->left_neighbour);
}

template <class Cell>
typename Mesh<Cell>::Iterator Mesh<Cell>::Iterator::operator+(const int i)
{
    Mesh::Iterator I(*this);
    if (i >= 0)
    {
        for (unsigned int j = 0; j < (unsigned int)i; j++)
        {
            ++I;
        }
    }
    else
    {
        for (unsigned int j = 0; j < (unsigned int)(-i); j++)
        {
            --I;
        }
    }

    return I;
}

template <class Cell>
typename Mesh<Cell>::Iterator Mesh<Cell>::Iterator::operator-(const int i)
{
    return *this + (-i);
}

template <class Cell>
Cell& Mesh<Cell>::Iterator::operator*()
{
    return *cell;   
}

template <class Cell>
Cell* Mesh<Cell>::Iterator::operator->()
{
    return cell;
}

template <class Cell>
Mesh<Cell>::Iterator::operator Cell*() const
{
    return cell;
}

template <class Cell>
bool Mesh<Cell>::Iterator::operator==(const Mesh<Cell>::Iterator& J)
{
    if (this->cell != J.cell)
    {
        return false;
    }
    else
    {
        return true;
    }
}

template <class Cell>
bool Mesh<Cell>::Iterator::operator!=(const Mesh::Iterator& J)
{
    return !(*this == J);
}

template <class Cell>
Cell* Mesh<Cell>::subdivide_at(Cell* C)
{
    double A_L = Math::Integrators::Gauss_Konrad_G7_K15(A_func, C->a, C->x()) / (C->x() - C->a);
    double A_R = Math::Integrators::Gauss_Konrad_G7_K15(A_func, C->x(), C->b) / (C->b - C->x());
    if (std::isnan(A_L))
    {
        A_L = C->A;
    }
    if (std::isnan(A_R))
    {
        A_R = C->A;
    }

    Cell* L = new Cell(C->U, C->a, C->x(), A_L, &QR);
    Cell* R = new Cell(C->U, C->x(), C->b, A_R, &QR);
    L->right_neighbour = R;
    R->left_neighbour = L;
    if (C->left_neighbour != nullptr)
    {
        L->left_neighbour = C->left_neighbour;
        C->left_neighbour->right_neighbour = L;
    }
    else
    {
        L->left_neighbour = nullptr;
        first_cell = L;
    }
    if (C->right_neighbour != nullptr)
    {
        R->right_neighbour = C->right_neighbour;
        C->right_neighbour->left_neighbour = R;
    }
    else
    {
        R->right_neighbour = nullptr;
        last_cell = R;
    }

    delete C;

    return R;
}

template <class Cell>
Cell* Mesh<Cell>::merge_cells(Cell* L, Cell* R)
{
    double beta_L = L->len();
    double beta_R = R->len();
    double beta = beta_L + beta_R;
    beta_L /= beta;
    beta_R /= beta;
    double A_C = Math::Integrators::Gauss_Konrad_G7_K15(A_func, L->a, R->b) / (R->b - L->a);
    Cell* C = new Cell(beta_L*L->U + beta_R*R->U, L->a, R->b,
        A_C, &QR);
    if (L->left_neighbour != nullptr)
    {
        L->left_neighbour->right_neighbour = C;
        C->left_neighbour = L->left_neighbour;
    }
    else
    {
        C->left_neighbour = nullptr;
        first_cell = C;
    }
    if (R->right_neighbour != nullptr)
    {
        R->right_neighbour->left_neighbour = C;
        C->right_neighbour = R->right_neighbour;
    }
    else
    {
        C->right_neighbour = nullptr;
        last_cell = C;
    }

    delete L;
    delete R;

    return C;
}

template <class Cell>
void Mesh<Cell>::calculate_variable_ranges()
{
    double U_norm_max = 0;
    double U_norm_min = std::numeric_limits<double>::max();

    for (Mesh<Cell>::Iterator C = begin(); C != end(); ++C)
    {
        double temp = (C->U/C->A).norm_2();
        if (temp > U_norm_max)
        {
            U_norm_max = temp;
        }
        else if (temp < U_norm_min)
        {
            U_norm_min = temp;
        }
    }
    ranges = U_norm_max - U_norm_min;
}

template <class Cell>
double Mesh<Cell>::detail(Cell* L, Cell* R) const
{
    return (R->U/R->A - L->U/L->A).norm_2() / ranges;
}

template <class Cell>
void Mesh<Cell>::optimize_mesh()
{
    if (first_cell != nullptr && last_cell != nullptr)
    {
        double length = last_cell->b - first_cell->a;
        calculate_variable_ranges();

        bool has_just_merged_cells = false;

        for (Mesh::Iterator C = begin(); C != end(); ++C)
        {
            // We determine whether the cell is on the boundary.
            bool is_boundary_cell = false;
            Mesh::Iterator D1 = C;
            Mesh::Iterator D2 = C;
            unsigned int j = 0;
            while (j < 4)
            {
                if (D1.cell == nullptr)
                {
                    is_boundary_cell = true;
                    break;
                }
                if (D2.cell == nullptr)
                {
                    is_boundary_cell = true;
                    break;
                }
                ++j;
                --D1;
                ++D2;
            }

            if (C->left_neighbour == nullptr)
            {
                if (detail(C, C->right_neighbour) >= detail_subdivide_threshold)
                {
                    C.cell = subdivide_at(C);
                }
            }
            else if (C->right_neighbour == nullptr)
            {
                double left_detail = detail(C->left_neighbour, C);
                if (left_detail >= detail_subdivide_threshold)
                {
                    C.cell = subdivide_at(C);
                }
                else if (left_detail <= detail_merge_threshold
                    && C->left_neighbour->len() + C->len() < boundary_cell_max_length_factor*length)
                {
                    C.cell = merge_cells(C->left_neighbour, C);
                }
            }
            else
            {
                Cell* L = C->left_neighbour;
                Cell* R = C->right_neighbour;
                double left_detail = detail(L, C);
                double right_detail = detail(C, R);
                if (fmax(left_detail, right_detail) >= detail_subdivide_threshold
                    && C->len() / 2 > min_length_factor*length)
                {
                    C.cell = subdivide_at(C);
                    has_just_merged_cells = false;
                }
                else if (left_detail <= detail_merge_threshold && !has_just_merged_cells
                    && ((!is_boundary_cell && C->left_neighbour->len() + C->len() < max_length_factor*length)
                    || (is_boundary_cell && C->left_neighbour->len() + C->len() < boundary_cell_max_length_factor*length)))
                {
                    C.cell = merge_cells(C->left_neighbour, C);
                    has_just_merged_cells = true;
                }
                else
                {
                    has_just_merged_cells = false;
                }
            }
        }   
    }
}

template <class Cell>
Mesh<Cell>::Mesh(Cell* first_cell, Cell* last_cell,
    const Chemistry::SolidGasReaction& QR,
    std::function<double(const double x)> A_func,
    const double detail_subdivide_threshold, const double detail_merge_threshold,
    const double max_length_factor, const double min_length_factor,
    const double boundary_cell_max_length_factor):
    first_cell(first_cell), last_cell(last_cell), QR(QR),
    detail_subdivide_threshold(detail_subdivide_threshold),
    detail_merge_threshold(detail_merge_threshold),
    max_length_factor(max_length_factor),
    min_length_factor(min_length_factor),
    boundary_cell_max_length_factor(boundary_cell_max_length_factor),
    A_func(A_func), ranges(0)
{

}

template <class Cell>
Mesh<Cell>::Mesh(const Mesh& mesh)
{
    first_cell = new Cell(*mesh.first_cell);
    Cell* current_cell;
    Mesh::Iterator C_old = mesh.begin()++;
    Mesh::Iterator C_new = begin();
    while(C_old != mesh.end())
    {
        current_cell = new Cell(*C_old);
        C_new->right_neighbour = current_cell;
        current_cell->left_neighbour = C_new;
        ++C_old; ++C_new;
    }
    last_cell = C_new;

    A_func = mesh.A_func;
    detail_subdivide_threshold = mesh.detail_subdivide_threshold;
    detail_merge_threshold = mesh.detail_merge_threshold;
    max_length_factor = mesh.max_length_factor;
    min_length_factor = mesh.min_length_factor;
    boundary_cell_max_length_factor = mesh.boundary_cell_max_length_factor;
    QR = mesh.QR;
    ranges = mesh.ranges;
}

template <class Cell>
Mesh<Cell>& Mesh<Cell>::operator=(const Mesh& mesh)
{
    if (this != &mesh)
    {
        free();
    
        if (mesh.first_cell != nullptr && mesh.last_cell != nullptr)
        {
            first_cell = new Cell(*mesh.first_cell);
            Cell* current_cell;
            Mesh::Iterator C_old = mesh.begin()++;
            Mesh::Iterator C_new = begin();
            while(C_old != mesh.end())
            {
                current_cell = new Cell(*C_old);
                C_new->right_neighbour = current_cell;
                current_cell->left_neighbour = C_new;
                ++C_old; ++C_new;
            }
            last_cell = C_new;
        }

        A_func = mesh.A_func;
        detail_subdivide_threshold = mesh.detail_subdivide_threshold;
        detail_merge_threshold = mesh.detail_merge_threshold;
        max_length_factor = mesh.max_length_factor;
        min_length_factor = mesh.min_length_factor;
        boundary_cell_max_length_factor = mesh.boundary_cell_max_length_factor;
        QR = mesh.QR;
        ranges = mesh.ranges;
    }
    return *this;
}

template <class Cell>
void Mesh<Cell>::free()
{
    if (last_cell != nullptr && first_cell != nullptr)
    {
        for (Mesh::Iterator C = rbegin()--; C != end(); --C)
        {
            delete C->right_neighbour;
            C->right_neighbour = nullptr;
        }
        last_cell = nullptr;
        delete first_cell;
        first_cell = nullptr;
    }
}

template <class Cell>
Mesh<Cell>::~Mesh()
{
    free();
}

template <class Cell>
std::vector<double> Mesh<Cell>::x() const
{
    std::vector<double> res;
    for (Mesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->x());
    }
    return res;
}

template <class Cell>
std::vector<double> Mesh<Cell>::x_partition() const
{
    std::vector<double> res;
    for (Mesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->a);
    }
    res.push_back(last_cell->b);
    return res;
}

template <class Cell>
std::vector<double> Mesh<Cell>::A() const
{
    std::vector<double> res;
    for (Mesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->A);
    }
    return res;
}

template <class Cell>
std::vector<Vector<double>> Mesh<Cell>::U() const
{
    std::vector<Vector<double>> res;
    for (Mesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->U);
    }
    return res;
}

template <class Cell>
void Mesh<Cell>::read_from_file(FILE* file)
{
    free();
    using ::read_from_file;
    read_from_file(detail_merge_threshold, file);
    read_from_file(detail_subdivide_threshold, file);
    read_from_file(max_length_factor, file);
    read_from_file(min_length_factor, file);
    read_from_file(boundary_cell_max_length_factor, file);
    unsigned int N_cells;
    read_from_file(N_cells, file);
    if (N_cells != 0)
    {
        QR = Chemistry::SolidGasReaction(file);
        Cell* C = new Cell(file, &QR);
        Cell* L = C;
        first_cell = C;
        for (unsigned int i = 1; i < N_cells; i++)
        {
            C = new Cell(file, &QR);
            L->right_neighbour = C;
            C->left_neighbour = L;
            L = C;
        }
        last_cell = C;
    }
}

template <class Cell>
void Mesh<Cell>::write_to_file(FILE* file) const
{
    using ::write_to_file;
    write_to_file(detail_merge_threshold, file);
    write_to_file(detail_subdivide_threshold, file);
    write_to_file(max_length_factor, file);
    write_to_file(min_length_factor, file);
    write_to_file(boundary_cell_max_length_factor, file);
    write_to_file(N_cells(), file);
    if (first_cell != nullptr)
    {
        write_to_file(QR, file);
        for (Mesh::Iterator C = begin(); C != end(); ++C)
        {
            write_to_file(*C, file);
        }
    }
}

GasMesh::GasMesh(GasCell* first_cell, GasCell* last_cell,
    const Chemistry::SolidGasReaction& QR,
    std::function<double(const double x)> A_func,
    const double detail_subdivide_threshold, const double detail_merge_threshold,
    const double max_length_factor, const double min_length_factor,
    const double boundary_cell_max_length_factor):
        Mesh<GasCell>(first_cell, last_cell, QR, A_func,
            detail_subdivide_threshold, detail_merge_threshold,
            max_length_factor, min_length_factor,
            boundary_cell_max_length_factor)
{

}

std::vector<double> GasMesh::c() const
{
    std::vector<double> res;
    for (GasMesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->c);
    }
    return res;
}

std::vector<double> GasMesh::rho() const
{
    std::vector<double> res;
    for (GasMesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->rho);
    }
    return res;
}

std::vector<double> GasMesh::v() const
{
    std::vector<double> res;
    for (GasMesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->v);
    }
    return res;
}

std::vector<double> GasMesh::T() const
{
    std::vector<double> res;
    for (GasMesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->T);
    }
    return res;
}

std::vector<double> GasMesh::P() const
{
    std::vector<double> res;
    for (GasMesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->P);
    }
    return res;
}

std::vector<double> GasMesh::M() const
{
    std::vector<double> res;
    for (GasMesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->M);
    }
    return res;
}

SolidMesh::SolidMesh(SolidCell* first_cell, SolidCell* last_cell,
    const Chemistry::SolidGasReaction& QR,
    std::function<double(const double x)> A_func,
    const double detail_subdivide_threshold, const double detail_merge_threshold,
    const double max_length_factor, const double min_length_factor,
    const double boundary_cell_max_length_factor):
        Mesh<SolidCell>(first_cell, last_cell, QR, A_func,
            detail_subdivide_threshold, detail_merge_threshold,
            max_length_factor, min_length_factor,
            boundary_cell_max_length_factor)
{

}

std::vector<double> SolidMesh::T() const
{
    std::vector<double> res;
    for (SolidMesh::Iterator C = begin(); C != end(); ++C)
    {
        res.push_back(C->T);
    }
    return res;
}

template class Mesh<GasCell>;
template class Mesh<SolidCell>;

namespace Utilities
{
    template class FileArray<GasMesh>;
    template class FileArray<SolidMesh>;
}

#endif // MESH_CPP
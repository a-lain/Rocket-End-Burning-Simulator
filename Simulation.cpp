#ifndef SIMULATION_CPP
#define SIMULATION_CPP

#include "Simulation.hpp"
#include "Math/Interpolation.hpp"
#include "CPGF/Plot2d/LinePlot.hpp"
#include "Utilities/FileOperations.hpp"
#include <future>
#include <limits>
#include <cmath>
#include <iostream>
#include <mutex>

using namespace Math;
using namespace Utilities;
using namespace CPGF::Plot2d;
using namespace Math::Interpolation;

GasBoundaryConditions::GasBoundaryConditions(const GasBoundaryConditionsType& left,
    const GasBoundaryConditionsType& right,
    std::function<double(double t)> left_condition_1,
    std::function<double(double t)> left_condition_2,
    std::function<double(double t)> left_condition_3,
    std::function<double(double t)> right_condition_1,
    std::function<double(double t)> right_condition_2,
    std::function<double(double t)> right_condition_3):
        left(left), right(right), left_condition_1(left_condition_1),
        left_condition_2(left_condition_2), left_condition_3(left_condition_3),
        right_condition_1(right_condition_1), right_condition_2(right_condition_2),
        right_condition_3(right_condition_3)
{

}

Simulation::Simulation
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
        std::function<double(std::function<double(double x)> f, const double a, const double b)> integrator,
        std::function<double(const GasCell& C, const double t)> external_forces,
        const unsigned int adaptive_refinement_period,
        const double detail_subdivide_threshold,
        const double detail_merge_threshold,
        const double max_length_factor,
        const double min_length_factor,
        const double boundary_cell_max_length_factor
):
    name(name),
    simulation_type(SimulationType::BOTH),
    CFL(CFL), N_saves(N_saves), adaptive_refinement_solid(adaptive_refinement_solid),
    adaptive_refinement_gas(adaptive_refinement_gas),
    adaptive_refinement_period(adaptive_refinement_period), N_tasks(N_tasks), external_forces(external_forces),
    convection_solver(convection_solver), diffusion_solver(diffusion_solver)
{
    // Step I. We create the solid mesh.
    // We create the partition of the X axis.
    double* x_part = new double[N_cells_solid+1];
    double x = a;
    double delta_x = (x_q - a) / N_cells_solid;
    x_part[0] = a;
    for (unsigned int i = 1; i < N_cells_solid; i++)
    {
        x += delta_x;
        x_part[i] = x;
    }
    x_part[N_cells_solid] = x_q;
    // We create the cells.
    instant_solid_mesh = SolidMesh(nullptr, nullptr, QR, A, detail_subdivide_threshold, detail_merge_threshold,
        max_length_factor, min_length_factor, boundary_cell_max_length_factor);    
    SolidCell* L_s = nullptr;
    SolidCell* C_s = nullptr;
    SolidCell* first_s = nullptr;
    for (unsigned int i = 0; i < N_cells_solid; i++)
    {
        L_s = C_s;
        Vector<double> U = Vector<double>({integrator([&] (double x){return 
                T(x)*A(x);},
                x_part[i], x_part[i+1]) / delta_x});
        C_s = new SolidCell(U, x_part[i], x_part[i+1],
            integrator(A, x_part[i], x_part[i+1]) / delta_x, &instant_solid_mesh.QR);
        if (i > 0)
        {
            L_s->right_neighbour = C_s;
            C_s->left_neighbour = L_s;   
        }
        else
        {
            first_s = C_s;
        }
    }
    instant_solid_mesh.first_cell = first_s;
    instant_solid_mesh.last_cell = C_s;
    delete[] x_part;

    // Step II. We create the gas mesh.
    x_part = new double[N_cells_gas+1];
    x = x_q;
    delta_x = (b - x_q) / N_cells_gas;
    x_part[0] = x_q;
    for (unsigned int i = 1; i < N_cells_gas; i++)
    {
        x += delta_x;
        x_part[i] = x;
    }
    x_part[N_cells_gas] = b;

    // We create the cells.
    instant_gas_mesh = GasMesh(nullptr, nullptr, QR, A, detail_subdivide_threshold, detail_merge_threshold,
        max_length_factor, min_length_factor, boundary_cell_max_length_factor);
    GasCell* L_g = nullptr;
    GasCell* C_g = nullptr;
    GasCell* first_g = nullptr;
    for (unsigned int i = 0; i < N_cells_gas; i++)
    {
        L_g = C_g;
        std::function<double(const double x)> rho = [&](const double x){return P(x)/(QR.R_g*T(x));};
        Vector<double> U = Vector<double>({integrator([&] (double x){return 
                rho(x)*A(x);},
                x_part[i], x_part[i+1]) / delta_x,
            integrator([&] (double x){return 
                rho(x)*v(x)*A(x);},
                x_part[i], x_part[i+1]) / delta_x,
            integrator([&] (double x){return
                ((1./(QR.gamma - 1)*P(x) + 1./2*rho(x)*v(x)*v(x))*A(x));},
                x_part[i], x_part[i+1]) / delta_x});
        C_g = new GasCell(U, x_part[i], x_part[i+1],
            integrator(A, x_part[i], x_part[i+1]) / delta_x, &instant_gas_mesh.QR);
        if (i > 0)
        {
            L_g->right_neighbour = C_g;
            C_g->left_neighbour = L_g;   
        }
        else
        {
            first_g = C_g;
        }
    }
    instant_gas_mesh.first_cell = first_g;
    instant_gas_mesh.last_cell = C_g;
    delete[] x_part;

    // We prepare the simulation.
    instant_t = 0;
    refine = 0;
    instant_x_q = x_q;
    instant_v_q = QR.v_q(instant_gas_mesh.first_cell->P);
    this->solid_BC = solid_BC;
    this->solid_BC.right = SolidBoundaryConditionsType::COMBUSTION_FRONT;
    this->solid_BC.right_condition = [&a = this->a, &b = this->b, &d = d](const double T, const double t)
    {
        if (t < 0)
        {
            return a;
        }
        else if (t > 0 && t < 1)
        {
            return b;
        }
        else
        {
            return d;
        }
    };
    this->gas_BC = gas_BC;
    this->gas_BC.left = GasBoundaryConditionsType::NONE;

    // We create the mesh.
    gas_mesh = new FileArray<GasMesh>(name + "-gas.dat");
    solid_mesh = new FileArray<SolidMesh>(name + "-solid.dat");
    gas_mesh->push_back(instant_gas_mesh);
    solid_mesh->push_back(instant_solid_mesh);
    t_array.push_back(0);
    x_q_array.push_back(instant_x_q);
    v_q_array.push_back(instant_v_q);
}

Simulation::Simulation
(
    const std::string& name,
    const Chemistry::SolidGasReaction& QR,
    const double a,
    const double b,
    const unsigned int N_cells,
    const bool adaptive_refinement_gas,
    const unsigned int N_tasks,
    const double CFL,
    const unsigned int N_saves,
    const GasBoundaryConditions& BC,
    std::function<double(double x)> rho,
    std::function<double(double x)> v,
    std::function<double(double x)> P,
    std::function<double(double x)> A,
    std::function<void(const GasCell& A, const GasCell& B, Math::Vector<double>* F, double* S_max)> convection_solver,
    std::function<double(std::function<double(double x)> f, const double a, const double b)> integrator,
    std::function<double(const GasCell& C, const double t)> external_forces,
    const unsigned int adaptive_refinement_period,
    const double detail_subdivide_threshold,
    const double detail_merge_threshold,
    const double max_length_factor,
    const double min_length_factor,
    const double boundary_cell_max_length_factor
):
    name(name),
    simulation_type(SimulationType::GAS),
    gas_BC(BC), CFL(CFL), N_saves(N_saves), adaptive_refinement_solid(false),
    adaptive_refinement_gas(adaptive_refinement_gas),
    adaptive_refinement_period(adaptive_refinement_period),
    N_tasks(N_tasks), external_forces(external_forces), convection_solver(convection_solver)
{
    // We create the partition of the X axis.
    double* x_part = new double[N_cells+1];
    double x = a;
    double delta_x = (b - a) / N_cells;
    x_part[0] = a;
    for (unsigned int i = 1; i < N_cells; i++)
    {
        x += delta_x;
        x_part[i] = x;
    }
    x_part[N_cells] = b;

    // We create the cells.
    instant_gas_mesh = GasMesh(nullptr, nullptr, QR, A, detail_subdivide_threshold, detail_merge_threshold,
        max_length_factor, min_length_factor, boundary_cell_max_length_factor);
    GasCell* L = nullptr;
    GasCell* C = nullptr;
    GasCell* first = nullptr;
    for (unsigned int i = 0; i < N_cells; i++)
    {
        L = C;
        Vector<double> U = Vector<double>({integrator([&] (double x){return 
                rho(x)*A(x);},
                x_part[i], x_part[i+1]) / delta_x,
            integrator([&] (double x){return 
                rho(x)*v(x)*A(x);},
                x_part[i], x_part[i+1]) / delta_x,
            integrator([&] (double x){return
                ((1./(QR.gamma - 1)*P(x) + 1./2*rho(x)*v(x)*v(x))*A(x));},
                x_part[i], x_part[i+1]) / delta_x});
        C = new GasCell(U, x_part[i], x_part[i+1],
            integrator(A, x_part[i], x_part[i+1]) / delta_x, &instant_gas_mesh.QR);
        if (i > 0)
        {
            L->right_neighbour = C;
            C->left_neighbour = L;   
        }
        else
        {
            first = C;
        }
    }
    instant_gas_mesh.first_cell = first;
    instant_gas_mesh.last_cell = C;

    // We prepare the simulation.
    instant_solid_mesh = SolidMesh();
    instant_t = 0;
    instant_v_q = 0;
    instant_x_q = 0;
    refine = 0;

    // We create the mesh.
    gas_mesh = new FileArray<GasMesh>(name + "-gas.dat");
    solid_mesh = new FileArray<SolidMesh>(name + "-solid.dat");
    gas_mesh->push_back(instant_gas_mesh);
    solid_mesh->push_back(instant_solid_mesh);
    t_array.push_back(0);    

    // We delete temporary variables.
    delete[] x_part;
}

Simulation::Simulation
(
    const std::string& name,
    const Chemistry::SolidGasReaction& QR,
    const double a,
    const double b,
    const unsigned int N_cells,
    const bool adaptive_refinement_solid,
    const unsigned int N_tasks,
    const double CFL,
    const unsigned int N_saves,
    const SolidBoundaryConditions& BC,
    std::function<double(double x)> T,
    std::function<double(double x)> A,
    std::function<void(SolidMesh& mesh, double& dt, const double CFL,
        const SolidBoundaryConditions BC, const double t)> diffusion_solver,
    std::function<double(std::function<double(double x)> f, const double a, const double b)> integrator,
    const unsigned int adaptive_refinement_period,
    const double detail_subdivide_threshold,
    const double detail_merge_threshold,
    const double max_length_factor,
    const double min_length_factor,
    const double boundary_cell_max_length_factor
):
    name(name),
    simulation_type(SimulationType::SOLID),
    solid_BC(BC), CFL(CFL), N_saves(N_saves), adaptive_refinement_solid(adaptive_refinement_solid),
    adaptive_refinement_gas(false),
    adaptive_refinement_period(adaptive_refinement_period), N_tasks(N_tasks), diffusion_solver(diffusion_solver)
{
    // We create the partition of the X axis.
    double* x_part = new double[N_cells+1];
    double x = a;
    double delta_x = (b - a) / N_cells;
    x_part[0] = a;
    for (unsigned int i = 1; i < N_cells; i++)
    {
        x += delta_x;
        x_part[i] = x;
    }
    x_part[N_cells] = b;

    // We create the cells.
    instant_solid_mesh = SolidMesh(nullptr, nullptr, QR, A, detail_subdivide_threshold, detail_merge_threshold,
        max_length_factor, min_length_factor, boundary_cell_max_length_factor);
    SolidCell* L = nullptr;
    SolidCell* C = nullptr;
    SolidCell* first = nullptr;
    for (unsigned int i = 0; i < N_cells; i++)
    {
        L = C;
        Vector<double> U = Vector<double>({integrator([&] (double x){return
                T(x)*A(x);}, x_part[i], x_part[i+1]) / delta_x});
        C = new SolidCell(U, x_part[i], x_part[i+1],
            integrator(A, x_part[i], x_part[i+1]) / delta_x, &instant_solid_mesh.QR);
        if (i > 0)
        {
            L->right_neighbour = C;
            C->left_neighbour = L;   
        }
        else
        {
            first = C;
        }
    }
    instant_solid_mesh.first_cell = first;
    instant_solid_mesh.last_cell = C;

    // We prepare the simulation.
    instant_gas_mesh = GasMesh();
    instant_t = 0;
    instant_x_q = 0;
    instant_v_q = 0;
    refine = 0;

    // We create the mesh.
    solid_mesh = new FileArray<SolidMesh>(name + "-solid.dat");
    gas_mesh = new FileArray<GasMesh>(name + "-gas.dat");
    solid_mesh->push_back(instant_solid_mesh);
    gas_mesh->push_back(instant_gas_mesh);
    t_array.push_back(0);
    x_q_array.push_back(instant_x_q);
    v_q_array.push_back(instant_v_q);

    // We delete temporary variables.
    delete[] x_part;
}

Simulation::~Simulation()
{
    delete gas_mesh;
    delete solid_mesh;
}

void Simulation::update(double dt)
{
    Chemistry::SolidGasReaction& QR = instant_gas_mesh.QR;
    instant_v_q = (simulation_type == SimulationType::BOTH) ? QR.v_q(instant_gas_mesh.first_cell->P) : 0;

    // Step I. We solve the gas.
    if (simulation_type != SimulationType::SOLID)
    {
        // We do things for parallelization. We determine the intermediate cells.
        unsigned int N_cells = instant_gas_mesh.N_cells();
        Vector<double>* fluxes = new Vector<double>[N_cells - 1];
        double* S_max = new double[N_cells - 1];
        std::future<void>* threads = new std::future<void>[N_tasks];
        unsigned int step = N_cells / N_tasks;
        GasMesh::Iterator* cells = new GasMesh::Iterator[N_tasks+1];
        unsigned int* pos = new unsigned int[N_tasks];
        cells[0] = instant_gas_mesh.begin();
        pos[0] = 0;
        for (unsigned int i = 1; i < N_tasks; i++)
        {
            cells[i] = cells[i - 1] + step;
            pos[i] = pos[i - 1] + step;
        }
        cells[N_tasks] = instant_gas_mesh.rbegin();

        // Step II.1: We calculate the intercell fluxes and maximum wave speeds (except the ones on the boundary).
        // We parallelize the computation.

        // We define the task.
        auto f1 = [=] (unsigned int i, const GasMesh::Iterator start, const GasMesh::Iterator end)
        {
            GasMesh::Iterator L = start;
            GasMesh::Iterator R = L++;

            while(R.cell != (*end).right_neighbour)
            {
                convection_solver(*L, *R, fluxes + i, S_max + i);
                i++;
                L = R;
                ++R;
            }
        };

        // We launch the tasks.
        for (unsigned int i = 0; i < N_tasks; i++)
        {
            threads[i] = std::async(std::launch::async, f1, pos[i], cells[i], cells[i+1]);
        }

        // We wait for all tasks to be comleted.
        for (unsigned int i = 0; i < N_tasks; i++)
        {
            threads[i].wait();
        }

        // Step II.2: We determine the max dt.
        unsigned int k = 0;
        for (GasMesh::Iterator C = instant_gas_mesh.begin()++; C != instant_gas_mesh.rbegin(); ++C)
        {
            double temp = CFL * C->len() / fmax(S_max[k], S_max[k+1]);
            if (temp < dt)
            {
                dt = temp;
            }
            k++;
        }
        // We correct dt if necessary.
        if (simulation_type == SimulationType::BOTH)
        {
            GasMesh::Iterator C_g = instant_gas_mesh.begin();
            SolidMesh::Iterator C_s = instant_solid_mesh.rbegin();
            while (fabs(QR.rho_s*C_s->A*instant_v_q - fluxes[0][0])*dt / (C_g->U[0]*C_g->len()) > 0.1
                || fabs(C_g->P*C_g->A - fluxes[0][1])*dt / fmax(fabs(C_g->U[1])*C_g->len(), 0.1*C_g->U[0]*C_g->len()) > 0.05)
            {
                dt /= 1.5;
            }
        }

        // Step II.3: We apply wall and periodic boundary conditions.
        GasCell* C = instant_gas_mesh.first_cell;
        GasCell* ghost_cell;
        double S;
        Vector<double> F;
        switch (gas_BC.left)
        {
            case GasBoundaryConditionsType::WALL:
                ghost_cell = new GasCell(*C);
                ghost_cell->U[1] *= -1;
                ghost_cell->update();
                convection_solver(*ghost_cell, *C, &F, &S);
                C->U += dt / C->len() * (F - fluxes[0]);
                C->U[1] = 0;
                delete ghost_cell;
                break;   
            case GasBoundaryConditionsType::PERIODIC:
                ghost_cell = new GasCell(*instant_gas_mesh.last_cell);
                ghost_cell->update();
                convection_solver(*ghost_cell, *C, &F, &S);
                C->U += dt / C->len() * (F - fluxes[0]);
                delete ghost_cell;
                break;
            default:
                break;
        }
        C->update();
        C = instant_gas_mesh.last_cell;
        switch (gas_BC.right)
        {
            case GasBoundaryConditionsType::WALL:
                ghost_cell = new GasCell(*C);
                ghost_cell->U[1] *= -1;
                ghost_cell->update();
                convection_solver(*C, *ghost_cell, &F, &S);
                C->U += dt / C->len() * (fluxes[N_cells - 2] - F);
                C->U[1] = 0;
                delete ghost_cell;
                break;
            case GasBoundaryConditionsType::PERIODIC:
                ghost_cell = new GasCell(*C);
                ghost_cell->U[1] *= -1;
                ghost_cell->update();
                convection_solver(*C, *ghost_cell, &F, &S);
                C->U += dt / C->len() * (fluxes[N_cells - 2] - F);
                C->U[1] = 0;      
                delete ghost_cell;      
            default:
                break;
        }
        C->update();

        // Step II.4: We update all cells with the fluxes.
        // We parallelize this computation.
        cells[0] = instant_gas_mesh.begin();
        cells[N_tasks] = instant_gas_mesh.rbegin()--;

        // We define the task.
        auto f2 = [=] (unsigned int i, const GasMesh::Iterator start, const GasMesh::Iterator end)
        {
            GasMesh::Iterator C = start;
            while (C != end)
            {
                C->U += dt / C->len() * (fluxes[i] - fluxes[i+1]);
                C->update();
                i++;
                ++C;
            }
        };
        
        // We launch the tasks.
        for (unsigned int i = 0; i < N_tasks; i++)
        {
            threads[i] = std::async(std::launch::async, f2, pos[i], cells[i]++, cells[i+1]++);
        }

        // We wait for all tasks to be comleted.
        for (unsigned int i = 0; i < N_tasks; i++)
        {
            threads[i].wait();
        }

        // Step II.5: source terms.
        cells[0] = instant_gas_mesh.begin();
        cells[N_tasks] = instant_gas_mesh.end();
        // We define the task.
        auto f4 = [=] (const GasMesh::Iterator start, const GasMesh::Iterator end)
        {
            unsigned int i = 0;
            GasMesh::Iterator C = start;
            while (C != end)
            {
                double dAdx = 0;
                unsigned int divisor = 0;
                if (C->left_neighbour != nullptr)
                {
                    dAdx += (C->A - C->left_neighbour->A) / (C->x() - C->left_neighbour->x());
                    divisor++;
                }
                if (C->right_neighbour != nullptr)
                {
                    dAdx += (C->right_neighbour->A - C->A) / (C->right_neighbour->x() - C->x());
                    divisor++;
                }
                dAdx /= divisor;

                C->U += dt * Vector<double>({0,
                    C->rho*external_forces(*C, instant_t)*C->A +
                        C->P*dAdx,
                    0});
                C->update();
                ++C;
                i++;
            }
        };
        
        // We launch the tasks.
        for (unsigned int i = 0; i < N_tasks; i++)
        {
            threads[i] = std::async(std::launch::async, f4, cells[i], cells[i+1]);
        }

        // We wait for all tasks to be comleted.
        for (unsigned int i = 0; i < N_tasks; i++)
        {
            threads[i].wait();
        }

        // Step II.6: we apply the boundary conditions.
        GasCell* first_cell = instant_gas_mesh.begin();
        GasCell* last_cell = instant_gas_mesh.rbegin();

        double alpha;
        double U_0;
        double v;
        double Q;
        double rho;
        double P;

        alpha = (first_cell->x() - first_cell->right_neighbour->x())
            / (first_cell->right_neighbour->right_neighbour->x() - first_cell->right_neighbour->x());
        switch (gas_BC.left)
        {
            case GasBoundaryConditionsType::FREE:
                first_cell->U = alpha*first_cell->right_neighbour->right_neighbour->U + (1 - alpha)*first_cell->right_neighbour->U;
                break;
            case GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY:
                Q = gas_BC.left_condition_1(instant_t + dt);
                rho = alpha*first_cell->right_neighbour->right_neighbour->rho + (1 - alpha)*first_cell->right_neighbour->rho;
                U_0 = rho*first_cell->A;
                v = Q / U_0;
                first_cell->U = Vector<double>({U_0,
                    Q,
                    (gas_BC.left_condition_2(instant_t + dt)/QR.gamma + (QR.gamma-1)/(2*QR.gamma)*v*v)*U_0});
                break;     
            case GasBoundaryConditionsType::FIXED_PRESSURE:
                rho = alpha*first_cell->right_neighbour->right_neighbour->rho + (1 - alpha)*first_cell->right_neighbour->rho;
                v = alpha*first_cell->right_neighbour->right_neighbour->v + (1 - alpha)*first_cell->right_neighbour->v;
                U_0 = rho*first_cell->A;
                first_cell->U = Vector<double>({U_0,
                    U_0*v,
                    1./(QR.gamma - 1)*gas_BC.left_condition_1(instant_t + dt)*first_cell->A
                        + v*v/2.*U_0});
                break;
            case GasBoundaryConditionsType::FIXED_DENSITY_SPEED_PRESSURE:
                rho = gas_BC.left_condition_1(instant_t + dt);
                v = gas_BC.left_condition_2(instant_t + dt);
                P = gas_BC.left_condition_3(instant_t + dt);
                first_cell->U = Vector<double>({rho*first_cell->A,
                    rho*v*first_cell->A,
                    (1./(QR.gamma - 1)*P + rho*v*v/2.)*first_cell->A});
                break;
            default:
                break;
        }
        first_cell->update();

        alpha = (last_cell->x() - last_cell->left_neighbour->left_neighbour->x()) 
                    / (last_cell->left_neighbour->x() - last_cell->left_neighbour->left_neighbour->x());
        switch (gas_BC.right)
        {
            case GasBoundaryConditionsType::FREE:
                last_cell->U = alpha*last_cell->left_neighbour->U + (1 - alpha)*last_cell->left_neighbour->left_neighbour->U;
                break;
            case GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY:
                Q = gas_BC.right_condition_1(instant_t + dt);
                rho = alpha*last_cell->left_neighbour->rho + (1 - alpha)*last_cell->left_neighbour->left_neighbour->rho;
                U_0 = rho*last_cell->A;
                v = Q / U_0;
                last_cell->U = Vector<double>({U_0,
                    Q,
                    (gas_BC.right_condition_2(instant_t + dt)/QR.gamma + (QR.gamma-1)/(2*QR.gamma)*v*v)*U_0});
                break;  
            case GasBoundaryConditionsType::FIXED_PRESSURE:
                rho = alpha*last_cell->left_neighbour->rho + (1 - alpha)*last_cell->left_neighbour->left_neighbour->rho;
                v = alpha*last_cell->left_neighbour->v + (1 - alpha)*last_cell->left_neighbour->left_neighbour->v;
                U_0 = rho*last_cell->A;
                last_cell->U = Vector<double>({U_0,
                    U_0*v,
                    1./(QR.gamma - 1)*gas_BC.right_condition_1(instant_t + dt)*last_cell->A
                        + v*v/2.*U_0});
                break;
            case GasBoundaryConditionsType::FIXED_DENSITY_SPEED_PRESSURE:
                rho = gas_BC.right_condition_1(instant_t + dt);
                v = gas_BC.right_condition_2(instant_t + dt);
                P = gas_BC.right_condition_3(instant_t + dt);
                last_cell->U = Vector<double>({rho*last_cell->A,
                    rho*v*last_cell->A,
                    (1./(QR.gamma - 1)*P + rho*v*v/2.)*last_cell->A});
                break;
            default:
                break;
        }
        last_cell->update();

        // Step II.7: We update the values of the first gas cell (only when simulation_type == SimulationType::BOTH).
        if (simulation_type == SimulationType::BOTH)
        {
            // Step I: We determine the new position of the combustion front and join the necessary solid cells.
            double old_x_q = instant_x_q;
            instant_v_q = QR.v_q(instant_gas_mesh.first_cell->P);
            instant_x_q -= instant_v_q*dt;
            double delta_x = old_x_q - instant_x_q;

            // We move the solid.
            SolidMesh::Iterator C_s = instant_solid_mesh.rbegin();
            GasMesh::Iterator C_g = instant_gas_mesh.begin();
            while (C_s->a > instant_x_q)
            {
                double gradient = (C_s->T - C_s->left_neighbour->T) / (C_s->x() - C_s->left_neighbour->x());
                double I_T = C_s->left_neighbour->left_neighbour->T*C_s->left_neighbour->left_neighbour->len()
                    + C_s->left_neighbour->T*C_s->left_neighbour->len()
                    + C_s->T*C_s->len();
                C_s = instant_solid_mesh.merge_cells(C_s->left_neighbour, C_s);
                double T_n_left = (I_T - gradient*(C_s->x() - C_s->left_neighbour->x())*C_s->len()) / (C_s->len() + C_s->left_neighbour->len());
                double T_n = T_n_left + gradient*(C_s->x() - C_s->left_neighbour->x());
                C_s->U = T_n*C_s->A;
                C_s->left_neighbour->U = T_n_left*C_s->left_neighbour->A;
                C_s->update();
                C_s->left_neighbour->update();
            }

            C_g = instant_gas_mesh.begin();
            C_s = instant_solid_mesh.rbegin();
            double old_rho_g = C_g->rho;
            rho_g = (C_g->U[0]*C_g->len() + (QR.rho_s*C_s->A*instant_v_q - fluxes[0][0])*dt)
                / (C_g->A*(C_g->len() + delta_x));
            double old_v_g = C_g->v;
            v_g = (old_rho_g*old_v_g*C_g->A*C_g->len() + (old_rho_g*QR.R_g*C_g->T*C_g->A - fluxes[0][1])*dt)
                / (rho_g*C_g->A*(C_g->len() + 2*delta_x) - QR.rho_s*C_s->A*delta_x);
            a = -QR.k_s*C_s->A*dt / (C_s->x() - C_s->left_neighbour->x()) / C_s->left_neighbour->A;
            b = (QR.rho_s*QR.cV_s*C_s->A*(C_s->len() - delta_x) + rho_g*QR.cV_g*C_g->A*(C_g->len() + delta_x)
                +QR.k_s*C_s->A*dt / (C_s->x() - C_s->left_neighbour->x()) ) / C_s->A;
            d = QR.rho_s*QR.cV_s*C_s->A*C_s->len()*C_s->T + old_rho_g*QR.cV_g*C_g->A*C_g->len()*C_g->T
                + C_g->A/2 * (old_rho_g*old_v_g*old_v_g*C_g->len() - rho_g*v_g*v_g*(C_g->len() + delta_x))
                + (QR.rho_s*QR.delta_H*C_s->A*instant_v_q - fluxes[0][2])*dt;
        }

        delete[] fluxes;
        delete[] S_max;
        delete[] threads;
        delete[] pos;
        delete[] cells;
    }

    // Step III: Solid.
    if (simulation_type != SimulationType::GAS)
    {
        diffusion_solver(instant_solid_mesh, dt, CFL, solid_BC, instant_t);
    }

    // Step IV: We advance the combustion front.
    if (simulation_type == SimulationType::BOTH)
    {
        SolidMesh::Iterator C_s = instant_solid_mesh.rbegin();
        GasMesh::Iterator C_g = instant_gas_mesh.begin();
        C_g->U = Vector<double>({rho_g*C_g->A,
            rho_g*v_g*C_g->A,
            rho_g*C_g->A*(QR.cV_g*C_s->T + v_g*v_g/2)});
        C_g->update();

        if ((instant_x_q - instant_solid_mesh.first_cell->a) /
            (instant_gas_mesh.last_cell->b - instant_solid_mesh.first_cell->a) < 1./1000)
        {
            instant_solid_mesh = SolidMesh();
            simulation_type = SimulationType::GAS;
            gas_BC.left = GasBoundaryConditionsType::WALL;
        }
        else
        {
            C_s = instant_solid_mesh.rbegin();
            C_s->b = instant_x_q;

            // We move the left boundary of the first gas cell.
            C_g = instant_gas_mesh.begin();
            double dx = (instant_gas_mesh.last_cell->b - instant_gas_mesh.first_cell->a) / instant_gas_mesh.N_cells();
            while (C_g->b - instant_x_q > dx)
            {
                C_g->a = C_g->b - dx;
                double A_g = Math::Integrators::Gauss_Konrad_G7_K15(instant_gas_mesh.A_func, C_g->a - dx, C_g->a) / dx;
                GasCell* L_g = new GasCell(C_g->U, C_g->a - dx, C_g->a, A_g, C_g->QR);
                C_g->left_neighbour = L_g;
                L_g->right_neighbour = C_g;
                L_g->left_neighbour = nullptr;
                --C_g;
            }
            C_g->a = instant_x_q;
            instant_gas_mesh.first_cell = C_g;

            // std::cout << "t = " << instant_t <<  ", rho = " << C_g->rho << ", v = " << C_g->v << ", T = " << C_g->T
            //     << ", T_left = " << C_s->left_neighbour->T
            //     << ", P = " << C_g->P << ", v_q = " << instant_v_q << "\n";   
            // std::cout << "BC1: " << QR.rho_s*instant_v_q - C_g->rho*(C_g->v + instant_v_q) << "\n";
            // std::cout << "gradient = " << (C_s->T - C_s->left_neighbour->T) / (C_s->x() - C_s->left_neighbour->x()) << "\n";
            // std::cout << "BC2: " << QR.k_s*(C_s->T - C_s->left_neighbour->T) / (C_s->x() - C_s->left_neighbour->x())
            //     -QR.rho_s*instant_v_q*((QR.cV_s - QR.cV_g)*C_s->T + QR.delta_H - C_g->v*C_g->v/2) + C_g->P*C_g->v << "\n";
        }
    }

    // double mass = 0;
    // double E = 0;
    // for (SolidMesh::Iterator C = instant_solid_mesh.begin(); C != instant_solid_mesh.end(); ++C)
    // {
    //     mass += QR.rho_s*C->A*(C->b - C->a);
    //     E += QR.rho_s*(QR.cV_s*C->U[0] + QR.delta_H*C->A)*(C->b - C->a);
    // }
    // for (GasMesh::Iterator C = instant_gas_mesh.begin(); C != instant_gas_mesh.end(); ++C)
    // {
    //     mass += C->U[0]*(C->b - C->a);
    //     E += C->U[2]*(C->b - C->a);
    // }
    
    // Step V: we update time.
    instant_t += dt;

    // std::cout << "t = " << instant_t << ", m = " << mass << ", E = " << E << "\n";

    // Step VI: we refine the mesh.
    if (refine == adaptive_refinement_period - 1)
    {
        if (adaptive_refinement_gas)
        {
            instant_gas_mesh.optimize_mesh();
        }
        if (adaptive_refinement_solid)
        {
            instant_solid_mesh.optimize_mesh();
        }
        refine = 0;
    }
    else
    {
        refine++;
    }
}

void Simulation::simulate_until(const double t)
{
    this->t_array.clear();
    this->x_q_array.clear();
    this->v_q_array.clear();
    gas_mesh->clear();
    solid_mesh->clear();
    this->t_array.push_back(instant_t);
    this->x_q_array.push_back(instant_x_q);
    this->v_q_array.push_back(instant_v_q);
    gas_mesh->push_back(instant_gas_mesh);
    solid_mesh->push_back(instant_solid_mesh);
    progress = Progress<double>("Simulation", ProgressEstimation::LINEAR, instant_t, t);
    progress.start();
    progress.update_to_terminal();
    refine = 0;

    double delta_t = t / (N_saves + 1);
    double _t = delta_t;
    for (unsigned int i = 0; i < N_saves + 1; i++)
    {
        while(instant_t < _t)
        {
            update(_t - instant_t);
            progress = instant_t;
        }
        gas_mesh->push_back(instant_gas_mesh);
        solid_mesh->push_back(instant_solid_mesh);
        this->t_array.push_back(instant_t);
        this->x_q_array.push_back(instant_x_q);
        this->v_q_array.push_back(instant_v_q);
        _t += delta_t;
    }

    // We always save the last state.
    // gas_mesh->push_back(instant_gas_mesh);
    // solid_mesh->push_back(instant_solid_mesh);
    // this->t_array.push_back(instant_t);
    // this->x_q_array.push_back(instant_x_q);
    // this->v_q_array.push_back(instant_v_q);
    progress.finish();
}

double Simulation::x_q(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        return x_q_array[0]; 
    }
    else if (t > t_array[t_array.size() - 1])
    {
        return x_q_array[t_array.size() - 1];
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;
        return alpha * x_q_array[i] + beta * x_q_array[i+1];
    }
}

double Simulation::v_q(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        return v_q_array[0]; 
    }
    else if (t > t_array[t_array.size() - 1])
    {
        return v_q_array[t_array.size() - 1];
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;
        return alpha * v_q_array[i] + beta * v_q_array[i+1];
    }
}

std::function<double(const double x)> Simulation::A(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        SolidMesh solid_mesh_1 = (*solid_mesh)[0];
        GasMesh gas_mesh_1 = (*gas_mesh)[0];
        std::vector<double> A_1;
        std::vector<double> x_partition_1;
        for (SolidMesh::Iterator C = solid_mesh_1.begin(); C != solid_mesh_1.end(); ++C)
        {
            A_1.push_back(C->A);
            x_partition_1.push_back(C->a);
        }
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            A_1.push_back(C->A);
            x_partition_1.push_back(C->a);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            x_partition_1.push_back(solid_mesh_1.last_cell->b);
        }
        else
        {
            x_partition_1.push_back(gas_mesh_1.last_cell->b);
        }

        AverageLinearInterpolation I(A_1, x_partition_1);
        return [=](const double x)
        {
            return I(x);
        };     
    }
    else if (t > t_array[t_array.size() - 1])
    {
        SolidMesh solid_mesh_2 = (*solid_mesh)[t_array.size() - 1];
        GasMesh gas_mesh_2 = (*gas_mesh)[t_array.size() - 1];
        std::vector<double> A_2;
        std::vector<double> x_partition_2;
        for (SolidMesh::Iterator C = solid_mesh_2.begin(); C != solid_mesh_2.end(); ++C)
        {
            A_2.push_back(C->A);
            x_partition_2.push_back(C->a);
        }
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            A_2.push_back(C->A);
            x_partition_2.push_back(C->a);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            x_partition_2.push_back(solid_mesh_2.last_cell->b);
        }
        else
        {
            x_partition_2.push_back(gas_mesh_2.last_cell->b);
        }
        AverageLinearInterpolation I(A_2, x_partition_2);
        return [=](const double x)
        {
            return I(x);
        }; 
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        SolidMesh solid_mesh_1 = (*solid_mesh)[i];
        GasMesh gas_mesh_1 = (*gas_mesh)[i];
        std::vector<double> A_1;
        std::vector<double> x_partition_1;
        for (SolidMesh::Iterator C = solid_mesh_1.begin(); C != solid_mesh_1.end(); ++C)
        {
            A_1.push_back(C->A);
            x_partition_1.push_back(C->a);
        }
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            A_1.push_back(C->A);
            x_partition_1.push_back(C->a);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            x_partition_1.push_back(solid_mesh_1.last_cell->b);
        }
        else
        {
            x_partition_1.push_back(gas_mesh_1.last_cell->b);
        }
        AverageLinearInterpolation I1(A_1, x_partition_1);
        SolidMesh solid_mesh_2 = (*solid_mesh)[i+1];
        GasMesh gas_mesh_2 = (*gas_mesh)[i+1];
        std::vector<double> A_2;
        std::vector<double> x_partition_2;
        for (SolidMesh::Iterator C = solid_mesh_2.begin(); C != solid_mesh_2.end(); ++C)
        {
            A_2.push_back(C->A);
            x_partition_2.push_back(C->a);
        }
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            A_2.push_back(C->A);
            x_partition_2.push_back(C->a);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            x_partition_2.push_back(solid_mesh_2.last_cell->b);
        }
        else
        {
            x_partition_2.push_back(gas_mesh_2.last_cell->b);
        }
        AverageLinearInterpolation I2(A_2, x_partition_2);
        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;

        return [=](const double x)
        {
            return alpha * I1(x) + beta * I2(x);
        }; 
    }
}

std::function<double(const double x)> Simulation::rho(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        GasMesh gas_mesh_1 = (*gas_mesh)[0];
        std::vector<double> rho_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            rho_1.push_back(C->rho);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I(rho_1, x_partition_1);
        return [=](const double x)
        {
            return I(x);
        };     
    }
    else if (t > t_array[t_array.size() - 1])
    {
        GasMesh gas_mesh_2 = (*gas_mesh)[t_array.size() - 1];
        std::vector<double> rho_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            rho_2.push_back(C->rho);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I(rho_2, x_partition_2);
        return [=](const double x)
        {
            return I(x);
        }; 
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        GasMesh gas_mesh_1 = (*gas_mesh)[i];
        std::vector<double> rho_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            rho_1.push_back(C->rho);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I1(rho_1, x_partition_1);
        GasMesh gas_mesh_2 = (*gas_mesh)[i+1];
        std::vector<double> rho_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            rho_2.push_back(C->rho);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I2(rho_2, x_partition_2);
        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;

        return [=](const double x)
        {
            return alpha * I1(x) + beta * I2(x);
        }; 
    }
}

std::function<double(const double x)> Simulation::c(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        GasMesh gas_mesh_1 = (*gas_mesh)[0];
        std::vector<double> c_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            c_1.push_back(C->c);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I(c_1, x_partition_1);
        return [=](const double x)
        {
            return I(x);
        };     
    }
    else if (t > t_array[t_array.size() - 1])
    {
        GasMesh gas_mesh_2 = (*gas_mesh)[t_array.size() - 1];
        std::vector<double> c_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            c_2.push_back(C->c);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I(c_2, x_partition_2);
        return [=](const double x)
        {
            return I(x);
        }; 
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        GasMesh gas_mesh_1 = (*gas_mesh)[i];
        std::vector<double> c_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            c_1.push_back(C->c);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I1(c_1, x_partition_1);
        GasMesh gas_mesh_2 = (*gas_mesh)[i+1];
        std::vector<double> c_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            c_2.push_back(C->c);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I2(c_2, x_partition_2);
        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;

        return [=](const double x)
        {
            return alpha * I1(x) + beta * I2(x);
        }; 
    }
}

std::function<double(const double x)> Simulation::v(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        GasMesh gas_mesh_1 = (*gas_mesh)[0];
        std::vector<double> v_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            v_1.push_back(C->v);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I(v_1, x_partition_1);
        return [=](const double x)
        {
            return I(x);
        };     
    }
    else if (t > t_array[t_array.size() - 1])
    {
        GasMesh gas_mesh_2 = (*gas_mesh)[t_array.size() - 1];
        std::vector<double> v_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            v_2.push_back(C->v);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I(v_2, x_partition_2);
        return [=](const double x)
        {
            return I(x);
        }; 
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        GasMesh gas_mesh_1 = (*gas_mesh)[i];
        std::vector<double> v_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            v_1.push_back(C->v);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I1(v_1, x_partition_1);
        GasMesh gas_mesh_2 = (*gas_mesh)[i+1];
        std::vector<double> v_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            v_2.push_back(C->v);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I2(v_2, x_partition_2);
        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;

        return [=](const double x)
        {
            return alpha * I1(x) + beta * I2(x);
        }; 
    }
}

std::function<double(const double x)> Simulation::T(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        SolidMesh solid_mesh_1 = (*solid_mesh)[0];
        GasMesh gas_mesh_1 = (*gas_mesh)[0];
        std::vector<double> T_1;
        std::vector<double> x_partition_1;
        for (SolidMesh::Iterator C = solid_mesh_1.begin(); C != solid_mesh_1.end(); ++C)
        {
            T_1.push_back(C->T);
            x_partition_1.push_back(C->a);
        }
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            T_1.push_back(C->T);
            x_partition_1.push_back(C->a);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            x_partition_1.push_back(solid_mesh_1.last_cell->b);
        }
        else
        {
            x_partition_1.push_back(gas_mesh_1.last_cell->b);
        }

        AverageLinearInterpolation I(T_1, x_partition_1);
        return [=](const double x)
        {
            return I(x);
        };     
    }
    else if (t > t_array[t_array.size() - 1])
    {
        SolidMesh solid_mesh_2 = (*solid_mesh)[t_array.size() - 1];
        GasMesh gas_mesh_2 = (*gas_mesh)[t_array.size() - 1];
        std::vector<double> T_2;
        std::vector<double> x_partition_2;
        for (SolidMesh::Iterator C = solid_mesh_2.begin(); C != solid_mesh_2.end(); ++C)
        {
            T_2.push_back(C->T);
            x_partition_2.push_back(C->a);
        }
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            T_2.push_back(C->T);
            x_partition_2.push_back(C->a);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            x_partition_2.push_back(solid_mesh_2.last_cell->b);
        }
        else
        {
            x_partition_2.push_back(gas_mesh_2.last_cell->b);
        }
        AverageLinearInterpolation I(T_2, x_partition_2);
        return [=](const double x)
        {
            return I(x);
        }; 
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        SolidMesh solid_mesh_1 = (*solid_mesh)[i];
        GasMesh gas_mesh_1 = (*gas_mesh)[i];
        std::vector<double> T_1;
        std::vector<double> x_partition_1;
        for (SolidMesh::Iterator C = solid_mesh_1.begin(); C != solid_mesh_1.end(); ++C)
        {
            T_1.push_back(C->T);
            x_partition_1.push_back(C->a);
        }
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            T_1.push_back(C->T);
            x_partition_1.push_back(C->a);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            x_partition_1.push_back(solid_mesh_1.last_cell->b);
        }
        else
        {
            x_partition_1.push_back(gas_mesh_1.last_cell->b);
        }
        AverageLinearInterpolation I1(T_1, x_partition_1);
        SolidMesh solid_mesh_2 = (*solid_mesh)[i+1];
        GasMesh gas_mesh_2 = (*gas_mesh)[i+1];
        std::vector<double> T_2;
        std::vector<double> x_partition_2;
        for (SolidMesh::Iterator C = solid_mesh_2.begin(); C != solid_mesh_2.end(); ++C)
        {
            T_2.push_back(C->T);
            x_partition_2.push_back(C->a);
        }
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            T_2.push_back(C->T);
            x_partition_2.push_back(C->a);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            x_partition_2.push_back(solid_mesh_2.last_cell->b);
        }
        else
        {
            x_partition_2.push_back(gas_mesh_2.last_cell->b);
        }
        AverageLinearInterpolation I2(T_2, x_partition_2);
        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;

        return [=](const double x)
        {
            return alpha * I1(x) + beta * I2(x);
        }; 
    }
}

std::function<double(const double x)> Simulation::P(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        GasMesh gas_mesh_1 = (*gas_mesh)[0];
        std::vector<double> P_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            P_1.push_back(C->rho);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I(P_1, x_partition_1);
        return [=](const double x)
        {
            return I(x);
        };     
    }
    else if (t > t_array[t_array.size() - 1])
    {
        GasMesh gas_mesh_2 = (*gas_mesh)[t_array.size() - 1];
        std::vector<double> P_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            P_2.push_back(C->P);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I(P_2, x_partition_2);
        return [=](const double x)
        {
            return I(x);
        }; 
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        GasMesh gas_mesh_1 = (*gas_mesh)[i];
        std::vector<double> P_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            P_1.push_back(C->P);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I1(P_1, x_partition_1);
        SolidMesh solid_mesh_2 = (*solid_mesh)[i+1];
        GasMesh gas_mesh_2 = (*gas_mesh)[i+1];
        std::vector<double> P_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            P_2.push_back(C->P);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I2(P_2, x_partition_2);
        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;

        return [=](const double x)
        {
            return alpha * I1(x) + beta * I2(x);
        }; 
    }
}

std::function<double(const double x)> Simulation::M(const double t) const
{
    // We apply bisection for time.
    if (t < 0)
    {
        GasMesh gas_mesh_1 = (*gas_mesh)[0];
        std::vector<double> M_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            M_1.push_back(C->M);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I(M_1, x_partition_1);
        return [=](const double x)
        {
            return I(x);
        };     
    }
    else if (t > t_array[t_array.size() - 1])
    {
        GasMesh gas_mesh_2 = (*gas_mesh)[t_array.size() - 1];
        std::vector<double> M_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            M_2.push_back(C->M);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I(M_2, x_partition_2);
        return [=](const double x)
        {
            return I(x);
        }; 
    }
    else
    {
        unsigned int left = 0;
        unsigned int right = t_array.size();
        unsigned int i = (left + right) / 2;
        while (true)
        {
            if (t > t_array[i+1])
            {
                left = i;
            }
            else if (t < t_array[i])
            {
                right = i;
            }
            else
            {
                break;
            }

            i = (left + right) / 2;
        }

        GasMesh gas_mesh_1 = (*gas_mesh)[i];
        std::vector<double> M_1;
        std::vector<double> x_partition_1;
        for (GasMesh::Iterator C = gas_mesh_1.begin(); C != gas_mesh_1.end(); ++C)
        {
            M_1.push_back(C->M);
            x_partition_1.push_back(C->a);
        }
        x_partition_1.push_back(gas_mesh_1.last_cell->b);
        AverageLinearInterpolation I1(M_1, x_partition_1);
        GasMesh gas_mesh_2 = (*gas_mesh)[i+1];
        std::vector<double> M_2;
        std::vector<double> x_partition_2;
        for (GasMesh::Iterator C = gas_mesh_2.begin(); C != gas_mesh_2.end(); ++C)
        {
            M_2.push_back(C->M);
            x_partition_2.push_back(C->a);
        }
        x_partition_2.push_back(gas_mesh_2.last_cell->b);
        AverageLinearInterpolation I2(M_2, x_partition_2);
        double delta_t = t_array[i+1] - t_array[i];
        double alpha = (t - t_array[i]) / delta_t;
        double beta = (t_array[i+1] - t) / delta_t;

        return [=](const double x)
        {
            return alpha * I1(x) + beta * I2(x);
        }; 
    }
}

Graphic* Simulation::mesh_plot() const
{
    Axis* x_axis = new Axis(AxisType::HORIZONTAL, "$x$");
    Axis* t_axis = new Axis(AxisType::VERTICAL, "$t$");

    Graphic* G = new Graphic();
    GraphicObject* g;

    for (unsigned int i = 0; i < gas_mesh->size() - 1; i++)
    {
        // We copy the mesh from disk to RAM.
        SolidMesh tmp_solid_mesh = (*solid_mesh)[i];
        GasMesh tmp_gas_mesh = (*gas_mesh)[i];
        for (SolidMesh::Iterator C = tmp_solid_mesh.begin(); C != tmp_solid_mesh.end(); ++C)
        {
            g = new LinePlot(std::vector<double>({t_array[i], t_array[i+1]}),
                std::vector<double>({C->a, C->a}), CPGF::Color::BLUE, CPGF::LineWidth::ULTRA_THIN);
            G->add(g, t_axis, x_axis);
        }
        for (GasMesh::Iterator C = tmp_gas_mesh.begin(); C != tmp_gas_mesh.end(); ++C)
        {
            g = new LinePlot(std::vector<double>({t_array[i], t_array[i+1]}),
                std::vector<double>({C->a, C->a}), CPGF::Color::BLUE, CPGF::LineWidth::ULTRA_THIN);
            G->add(g, t_axis, x_axis);
        }
        if (simulation_type == SimulationType::SOLID)
        {
            g = new LinePlot(std::vector<double>({t_array[i], t_array[i+1]}),
            std::vector<double>({tmp_solid_mesh.last_cell->b, tmp_solid_mesh.last_cell->b}),
            CPGF::Color::BLUE, CPGF::LineWidth::ULTRA_THIN);
        }
        else
        {
            g = new LinePlot(std::vector<double>({t_array[i], t_array[i+1]}),
            std::vector<double>({tmp_gas_mesh.last_cell->b, tmp_gas_mesh.last_cell->b}),
            CPGF::Color::BLUE, CPGF::LineWidth::ULTRA_THIN);
        }
        G->add(g, t_axis, x_axis);
    }

    return G;
}

Simulation::Simulation(const std::string& file):
    refine(false),
    a(0), b(0), d(0),  rho_g(0), v_g(0)
{
    using ::read_from_file;
    FILE* F = fopen(file.c_str(), "rb");
    read_from_file(name, F);
    gas_mesh = new FileArray<GasMesh>(name + "-gas.dat");
    solid_mesh = new FileArray<SolidMesh>(name + "-solid.dat");
    unsigned int type;
    read_from_file(type, F);
    simulation_type = (SimulationType) type;
    // The next line is kept for downgrade compatiblity.
    Chemistry::SolidGasReaction QR = Chemistry::SolidGasReaction(F);
    read_from_file(x_q_array, F);
    read_from_file(v_q_array, F);
    read_from_file(t_array, F);
    read_from_file(CFL, F);
    read_from_file(N_saves, F);
    read_from_file(adaptive_refinement_solid, F);
    read_from_file(adaptive_refinement_gas, F);
    read_from_file(adaptive_refinement_period, F);
    read_from_file(N_tasks, F);
    fclose(F);

    unsigned int N = t_array.size() - 1;
    instant_t = t_array[N];
    instant_v_q = v_q_array[N];
    instant_x_q = x_q_array[N];
    instant_gas_mesh = (*gas_mesh)[N];
    instant_solid_mesh = (*solid_mesh)[N];
}

void Simulation::write_to_file(const std::string& file) const
{
    using ::write_to_file;
    FILE* F = fopen(file.c_str(), "wb");
    write_to_file(name, F);
    write_to_file((unsigned int) simulation_type, F);
    // This next line is kept for downward compatibility.
    write_to_file(Chemistry::SolidGasReaction(), F);
    write_to_file(x_q_array, F);
    write_to_file(v_q_array, F);
    write_to_file(t_array, F);
    write_to_file(CFL, F);
    write_to_file(N_saves, F);
    write_to_file(adaptive_refinement_solid, F);
    write_to_file(adaptive_refinement_gas, F);
    write_to_file(adaptive_refinement_period, F);
    write_to_file(N_tasks, F);
    fclose(F);
}

#endif // SIMULATION_CPP
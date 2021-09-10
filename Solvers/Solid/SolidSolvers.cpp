#ifndef SOLID_SOLVERS_CPP
#define SOLID_SOLVERS_CPP

#include "SolidSolvers.hpp"
#include <iostream>

using namespace Math;

SolidBoundaryConditions::SolidBoundaryConditions(const SolidBoundaryConditionsType left, const SolidBoundaryConditionsType right,
    std::function<double(double T, double t)> left_condition, std::function<double(double T, double t)> right_condition):
        left(left), right(right), left_condition(left_condition), right_condition(right_condition)
{

}

void Solvers::Solid::euler_explicit(SolidMesh& mesh, double& dt, const double CFL, const SolidBoundaryConditions BC, const double t)
{
    double _t = t;
    double _t_objective = t + dt;

    double alpha = mesh.first_cell->QR->alpha;

    while (_t < _t_objective)
    {
        // Step I. We determine dt.
        double _dt = _t_objective - _t;
        for (SolidMesh::Iterator C = mesh.begin(); C != mesh.end(); ++C)
        {
            double temp = CFL * C->len()*C->len() / (2*alpha);
            if (temp < _dt)
            {
                _dt = temp;
            }
        }

        // Step II. We update the cell values.
        double old_T = mesh.first_cell->T;
        if (BC.left == SolidBoundaryConditionsType::FIXED_TEMPERATURE)
        {
            SolidCell* C = mesh.first_cell;
            C->U[0] = BC.left_condition(0, _t) * C->A;
            C->update();
        }
        else
        {
            SolidCell* C = mesh.first_cell;
            double c = alpha*mesh.A_func(C->b) / (C->b - C->a) / (C->right_neighbour->x() - C->x());
            double b = - c;
            C->U[0] += _dt* (c*C->right_neighbour->T + b*C->T - mesh.A_func(C->a)*BC.left_condition(C->T, _t)/(C->b - C->a));
            C->update();
        }
        for (SolidMesh::Iterator C = mesh.begin()++; C != mesh.rbegin(); ++C)
        {
            double temp_T;
            double c = alpha*mesh.A_func(C->b) / (C->b - C->a) / (C->right_neighbour->x() - C->x());
            double a = alpha*mesh.A_func(C->a) / (C->b - C->a) / (C->x() - C->left_neighbour->x());
            double b = - (a + c);
            temp_T = C->T;
            C->U[0] += _dt * (c*C->right_neighbour->T + b*C->T + a*old_T);
            C->update();
            old_T = temp_T;
        }
        if (BC.right == SolidBoundaryConditionsType::FIXED_TEMPERATURE)
        {
            SolidCell* C = mesh.last_cell;
            C->U[0] = BC.right_condition(0, _t) * C->A;
            C->update();
        }
        else
        {
            SolidCell* C = mesh.last_cell;
            double a = alpha*mesh.A_func(C->a) / (C->b - C->a) / (C->x() - C->left_neighbour->x());
            double b = - a;
            C->U[0] += _dt* (b*C->T + a*old_T + mesh.A_func(C->b)*BC.right_condition(C->T, _t)/(C->b - C->a));
            C->update();
        }

        _t += _dt;
    }    
}

void Solvers::Solid::euler_implicit(SolidMesh& mesh, double& dt, const double CFL, const SolidBoundaryConditions BC, const double t)
{
    double alpha = mesh.first_cell->QR->alpha;

    // Step I. We construct the tridiagonal matrix.
    unsigned int N_cells = mesh.N_cells();
    double* a = new double[N_cells];
    double* b = new double[N_cells];
    double* c = new double[N_cells];
    double* d = new double[N_cells];

    if (BC.left == SolidBoundaryConditionsType::FIXED_TEMPERATURE)
    {
        SolidCell* C = mesh.first_cell;
        c[0] = 0;
        a[0] = 0;
        b[0] = 1;
        d[0] = BC.left_condition(0, t+dt)*C->A;
    }
    else
    {
        SolidCell* C = mesh.first_cell;
        double n = BC.left_condition(0, 0);
        double m = BC.left_condition(1, 0) - n;
        c[0] = -dt*alpha*mesh.A_func(C->b) / (C->b - C->a) / (C->right_neighbour->x() - C->x()) / C->right_neighbour->A;
        a[0] = 0;
        b[0] = dt*alpha/C->A*(mesh.A_func(C->b) / (C->b - C->a) / (C->right_neighbour->x() - C->x())
            + mesh.A_func(C->a)*m/(C->b - C->a)) + 1.;
        d[0] = C->U[0] - alpha*dt*mesh.A_func(C->a)*n / (C->b - C->a);
    }
    unsigned int i = 1;
    for (SolidMesh::Iterator C = mesh.begin()++; C != mesh.rbegin(); ++C)
    {
        c[i] = -dt*alpha*mesh.A_func(C->b) / (C->b - C->a) / (C->right_neighbour->x() - C->x()) / C->right_neighbour->A;
        a[i] = -dt*alpha*mesh.A_func(C->a) / (C->b - C->a) / (C->x() - C->left_neighbour->x()) / C->left_neighbour->A;
        b[i] = dt*alpha*(mesh.A_func(C->b) / (C->right_neighbour->x() - C->x())
            + mesh.A_func(C->a) / (C->x() - C->left_neighbour->x())) / (C->b - C->a) / C->A + 1.;
        d[i] = C->U[0];
        i++;
    }
    if (BC.right == SolidBoundaryConditionsType::FIXED_TEMPERATURE)
    {
        SolidCell* C = mesh.last_cell;
        c[N_cells - 1] = 0;
        a[N_cells - 1] = 0;
        b[N_cells - 1] = 1;
        d[N_cells - 1] = BC.right_condition(0, t+dt) * C->A;
    }
    else if (BC.right == SolidBoundaryConditionsType::COMBUSTION_FRONT)
    {
        c[N_cells - 1] = 0;
        b[N_cells - 1] = BC.right_condition(0, 0.5);
        a[N_cells - 1] = BC.right_condition(0, -1);
        d[N_cells - 1] = BC.right_condition(0, 2);
    }
    else
    {
        SolidCell* C = mesh.last_cell;
        double n = BC.right_condition(0, 0);
        double m = BC.right_condition(1, 0) - n;
        c[N_cells - 1] = 0;
        a[N_cells - 1] = -dt*alpha*mesh.A_func(C->a) / (C->b - C->a) / (C->x() - C->left_neighbour->x()) / C->left_neighbour->A;
        b[N_cells - 1] = dt*alpha/C->A*(mesh.A_func(C->a) / (C->b - C->a) / (C->x() - C->left_neighbour->x())
            - mesh.A_func(C->b)*m/(C->b - C->a)) + 1.;
        d[N_cells - 1] = C->U[0] + alpha*dt*mesh.A_func(C->b)*n / (C->b - C->a);
    }

    // Step II. We apply Gaussian Elimination (Tridiagonal Matrix algorithm).
    for (unsigned int i = 1; i < N_cells; i++)
    {
        double l = a[i] / b[i-1];
        b[i] = b[i] - l*c[i-1];
        d[i] = d[i] - l*d[i-1];
    }

    // Step III. We apply regressive substitution.
    mesh.last_cell->U[0] = d[N_cells - 1] / b[N_cells - 1];
    mesh.last_cell->update();
    i = N_cells - 2;
    for (SolidMesh::Iterator C = mesh.rbegin()--; C != mesh.end(); --C)
    {
        C->U[0] = (d[i] - c[i]*C->right_neighbour->U[0]) / b[i];
        C->update();
        i--;
    }

    // Step V. We free memory.
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
}

#endif // SOLID_SOLVERS_CPP
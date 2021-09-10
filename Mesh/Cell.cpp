#ifndef CELL_CPP
#define CELL_CPP

#include "Cell.hpp"
#include <iostream>
#include "../Utilities/FileOperations.hpp"

using namespace Math;

BaseCell::BaseCell(const Math::Vector<double>& U, const double a, const double b, const double A,
        Chemistry::SolidGasReaction* QR):
        U(U), a(a), b(b), A(A), QR(QR)
{

}

BaseCell::BaseCell(FILE* file, Chemistry::SolidGasReaction* QR):
    QR(QR)
{
    read_from_file(file);
}

double BaseCell::x() const
{
    return (b + a) / 2;
}

double BaseCell::len() const
{
    return b - a;
}

void BaseCell::read_from_file(FILE* file)
{
    using ::read_from_file;
    read_from_file(U, file);
    read_from_file(a, file);
    read_from_file(b, file);
    read_from_file(A, file);
}

void BaseCell::write_to_file(FILE* file) const
{
    using ::write_to_file;
    write_to_file(U, file);
    write_to_file(a, file);
    write_to_file(b, file);
    write_to_file(A, file);
}

std::string BaseCell::to_string() const
{
    char buffer[50];
    sprintf(buffer, "%p", (void*)this);
    std::string text = "Begin Cell\n\taddress = " + std::string(buffer);

    return text;
}

GasCell::GasCell(const Math::Vector<double>& U, double a, double b, const double A,
        Chemistry::SolidGasReaction* QR):
    BaseCell(U, a, b, A, QR), right_neighbour(nullptr), left_neighbour(nullptr),
        rho(0), c(0), v(0), T(0), P(0), H(0), E(0), M(0)
{
    update();
}

GasCell::GasCell(FILE* file, Chemistry::SolidGasReaction* QR):
    BaseCell(file, QR), right_neighbour(nullptr), left_neighbour(nullptr), rho(0), c(0), v(0), T(0), P(0), H(0), E(0), M(0)
{
    update();
}

Vector<double> GasCell::F() const
{
    double v2 = v*v;
    return Vector<double>({U[1], U[0]*(v2 + QR->R_g*T),
        U[1]*((QR->cV_g + QR->R_g)*T + v2/2)});
}

void GasCell::update()
{
    rho = U[0] / A;
    if (U[0] <= 1e-12)
    {
        rho = 1e-6;
        U[0] = 1e-6;
        U[1] = 0;
        v = 0;
        T = 0;
    }
    else
    {
        v = U[1] / U[0];
        T = (U[2] / U[0] - 1./2*v*v) / QR->cV_g;
    }
    P = rho * QR->R_g * T;
    c = sqrt(QR->gamma * QR->R_g * T);
    M = v / c;
    E = 1./2*rho*v*v + P/(QR->gamma - 1);
    H = (E + P) / rho;
}

void GasCell::read_from_file(FILE* file)
{
    BaseCell::read_from_file(file);
    update();
}

std::string GasCell::to_string() const
{
    std::string text = BaseCell::to_string();
    char buffer[50];
    sprintf(buffer, "%p", (void*)left_neighbour);
    text += "\n\tleft_neighbour = " + std::string(buffer);
    sprintf(buffer, "%p", (void*)right_neighbour);
    text += "\n\tright_neighbour = " + std::string(buffer);

    return text
        + "\n\tU = " + U.to_string()
        + "\n\ta = " + std::to_string(a)
        + "\n\tb = " + std::to_string(b)
        + "\n\tA = " + std::to_string(A)
        + "\n" + QR->to_string()
        + "\n\tρ = " + std::to_string(rho)
        + "\n\tc = " + std::to_string(c)
        + "\n\tv = " + std::to_string(v)
        + "\n\tT = " + std::to_string(T)
        + "\n\tP = " + std::to_string(P)
        + "\n\tH = " + std::to_string(H)
        + "\n\tE = " + std::to_string(E)
        + "\nEnd Cell";
}

SolidCell::SolidCell(const Math::Vector<double>& U, double a, double b, const double A,
        Chemistry::SolidGasReaction* QR):
    BaseCell(U, a, b, A, QR), right_neighbour(nullptr), left_neighbour(nullptr), T(0)
{
    update();
}

SolidCell::SolidCell(FILE* file, Chemistry::SolidGasReaction* QR):
    BaseCell(file, QR), right_neighbour(nullptr), left_neighbour(nullptr), T(0)
{
    update();
}

void SolidCell::update()
{
    T = U[0] / A;
}

void SolidCell::read_from_file(FILE* file)
{
    BaseCell::read_from_file(file);
    update();
}

std::string SolidCell::to_string() const
{
    std::string text = BaseCell::to_string();
    char buffer[50];
    sprintf(buffer, "%p", (void*)left_neighbour);
    text += "\n\tleft_neighbour = " + std::string(buffer);
    sprintf(buffer, "%p", (void*)right_neighbour);
    text += "\n\tright_neighbour = " + std::string(buffer);

    return text
        + "\n\tU = " + std::to_string(U[0])
        + "\n\ta = " + std::to_string(a)
        + "\n\tb = " + std::to_string(b)
        + "\n\tA = " + std::to_string(A)
        + "\n" + QR->to_string()
        + "\n\tT = " + std::to_string(T)
        + "\nEnd Cell";
}

#endif // CELL_CPP
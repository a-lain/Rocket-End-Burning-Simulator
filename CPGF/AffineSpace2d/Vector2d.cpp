#ifndef VECTOR2D_CPP
#define VECTOR2D_CPP

#include "Vector2d.hpp"
#include <cmath>

using namespace CPGF::AffineSpace;

namespace CPGF
{
    namespace AffineSpace
    {
        Vector2d operator+(const Vector2d& v, const Vector2d& w)
        {
            return Vector2d(v.x + w.x, v.y + w.y);
        }

        Vector2d operator-(const Vector2d& v, const Vector2d& w)
        {
            return Vector2d(v.x - w.x, v.y - w.y);
        }

        Vector2d operator*(const Vector2d& v, const Vector2d& w)
        {
            return Vector2d(v.x * w.x, v.y * w.y);
        }

        Vector2d operator/(const Vector2d& v, const Vector2d& w)
        {
            return Vector2d(v.x / w.x, v.y / w.y);
        }

        double operator|(const Vector2d& v, const Vector2d& w)
        {
            return v.x * w.x + v.y * w.y;
        }
    }
}

Vector2d& Vector2d::operator+()
{
    return *this;
}

Vector2d& Vector2d::operator-()
{
    x = -x;
    y = -y;
    return *this;
}

Vector2d& Vector2d::operator+=(const Vector2d& w)
{
    this->x += w.x;
    this->y += w.y;
    return *this;
}

Vector2d& Vector2d::operator-=(const Vector2d& w)
{
    this->x -= w.x;
    this->y -= w.y;
    return *this;
}

Vector2d& Vector2d::operator*=(const Vector2d& w)
{
    this->x *= w.x;
    this->y *= w.y;
    return *this;
}

Vector2d& Vector2d::operator/=(const Vector2d& w)
{
    this->x /= w.x;
    this->y /= w.y;
    return *this;
}

Vector2d Vector2d::perp() const
{
    return Vector2d(-y, x);
}

double Vector2d::norm() const
{
    return sqrt(*this | *this);
}

std::string Vector2d::to_string() const
{
    return "V(" + std::to_string(x) + ", " + std::to_string(y) + ")"; 
}

Vector2d::Vector2d() : Vector2d(0,0)
{

}

Vector2d::Vector2d(double x) : Vector2d(x,x)
{
    
}

Vector2d::Vector2d(double x, double y)
{
    this->x = x;
    this->y = y;
}

#endif // VECTOR2D_CPP
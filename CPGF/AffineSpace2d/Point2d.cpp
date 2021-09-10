#ifndef POINT2D_CPP
#define POINT2D_CPP

#include "Point2d.hpp"
#include <cmath>

using namespace CPGF::AffineSpace;

namespace CPGF
{
    namespace AffineSpace
    {
        Point2d operator+(const Point2d& P, const Vector2d& v)
        {
            return Point2d(P.x + v.x, P.y + v.y);
        }

        Vector2d operator-(const Point2d& P, const Point2d& Q)
        {
            return Vector2d(P.x - Q.x, P.y - Q.y);
        }

        bool operator==(const Point2d& P, const Point2d& Q)
        {
            return ((Q - P).norm() < 1e-10) ? true : false;
        }
    }
}

Point2d& Point2d::operator+=(const Vector2d& v)
{
    this->x += v.x;
    this->y += v.y;
    return *this;
}

double Point2d::angle_with_respect_to(const Point2d& Q)
{
    Vector2d v = *this - Q;
    double theta = acos(v.x / v.norm());
    return (v.y >= 0) ? theta : -theta + 2*M_PI;
}

Point2d& Point2d::rotate_with_respect_to(const Point2d& Q, const double theta)
{
    Vector2d v(Q - Point2d(0,0));
    *this += -v;
    *this = Point2d(cos(theta)*x + sin(theta)*y, -sin(theta)*x + cos(theta)*y);
    *this += v;
    return *this;
}    

Point2d& Point2d::rotate_to_with_respect_to(const Point2d& Q, const double theta)
{
    return rotate_with_respect_to(Q, theta - angle_with_respect_to(Q));
}

Point2d& Point2d::scale_with_respect_to(const Point2d& Q, const Vector2d& s)
{
    Vector2d v(Q - Point2d(0,0));
    *this += -v;
    x *= s.x;
    y *= s.y;
    *this += v;
    return *this;
}

std::string Point2d::to_string() const
{
    return "P(" + std::to_string(x) + ", " + std::to_string(y) + ")";
}

Point2d::Point2d() : Point2d(0,0)
{

}

Point2d::Point2d(double x, double y)
{
    this->x = x;
    this->y = y;
}

#endif // POINT2D_CPP
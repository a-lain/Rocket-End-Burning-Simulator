#ifndef OBJECT2D_CPP
#define OBJECT2D_CPP

#include "Object2d.hpp"

using namespace CPGF::AffineSpace;
using namespace CPGF::Objects2d;
using namespace CPGF::Basics;

Object2d::Object2d(): paths(std::vector<Path2d>())
{

}

Object2d::Object2d(const Basics::Path2d& path): paths(std::vector<Path2d>({path}))
{

}

Object2d::Object2d(const std::vector<Basics::Path2d>& paths): paths(paths)
{

}

unsigned int Object2d::size() const
{
    return paths.size();
}

Object2d& Object2d::translate(const Vector2d& v)
{
    for (unsigned int i = 0; i < paths.size(); i++)
    {
        paths[i].translate(v);
    }
    return *this;
}

Object2d& Object2d::rotate_with_respect_to(const Point2d& Q, const double theta)
{
    for (unsigned int i = 0; i < paths.size(); i++)
    {
        paths[i].rotate_with_respect_to(Q, theta);
    }
    return *this;
}

Object2d& Object2d::scale_with_respect_to(const Point2d& Q, const Vector2d& s)
{
    for (unsigned int i = 0; i < paths.size(); i++)
    {
        paths[i].scale_with_respect_to(Q, s);
    }
    return *this;
}

std::string Object2d::render_to_string() const
{
    std::string text = "";
    for (unsigned int i = 0; i < paths.size(); i++)
    {
        text += paths[i].render_to_string();
    }
    return text;
}

Object2d& Object2d::operator+=(const Object2d& B)
{
    for (unsigned int i = 0; i < B.size(); i++)
    {
        this->paths.push_back(B.paths[i]);
    }
    return *this;
}

Object2d& Object2d::operator+=(const Path2d& path)
{
    this->paths.push_back(path);
    return *this;
}

Object2d operator+(const Object2d& A, const Object2d& B)
{
    Object2d res(A);
    res += B;
    return res;
}

#endif // OBJECT2D_CPP

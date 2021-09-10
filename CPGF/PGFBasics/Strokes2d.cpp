#ifndef STROKES2D_CPP
#define STROKES2D_CPP

#include "Strokes2d.hpp"

using namespace CPGF::AffineSpace;
using namespace CPGF::Basics;

StraightStroke2d::StraightStroke2d()
{

}

StraightStroke2d::StraightStroke2d(const std::vector<Point2d> points): points(points)
{

}

StraightStroke2d* StraightStroke2d::clone() const
{
    return new StraightStroke2d(*this);
}

StraightStroke2d& StraightStroke2d::add_point(const Point2d& P)
{
    points.push_back(P);
    return *this;
}

StraightStroke2d& StraightStroke2d::operator+=(const Point2d& P)
{
    return this->add_point(P);
}

unsigned int StraightStroke2d::size() const
{
    return points.size();
}

Point2d& StraightStroke2d::start()
{
    return points[0];
}

Point2d StraightStroke2d::start() const
{
    return points[0];
}

Point2d& StraightStroke2d::end()
{
    return points[points.size() - 1];
}

Point2d StraightStroke2d::end() const
{
    return points[points.size() - 1];
}

StraightStroke2d& StraightStroke2d::translate(const Vector2d& v)
{
    for (unsigned int i = 0; i < points.size(); i++)
    {
        points[i] += v;
    }
    return *this;
}

StraightStroke2d& StraightStroke2d::rotate_with_respect_to(const Point2d& Q, const double theta)
{
    for (unsigned int i = 0; i < points.size(); i++)
    {
        points[i].rotate_with_respect_to(Q, theta);
    }
    return *this;
}

StraightStroke2d& StraightStroke2d::scale_with_respect_to(const AffineSpace::Point2d& Q, const AffineSpace::Vector2d& s)
{
    for (unsigned int i = 0; i < points.size(); i++)
    {
        points[i].scale_with_respect_to(Q, s);
    }
    return *this;
}

double StraightStroke2d::length() const
{
    if (size() <= 1)
    {
        return 0;
    }
    else
    {
        double res = 0;
        for (unsigned int i = 1; i < size(); i++)
        {
            res += (points[i] - points[i-1]).norm();
        }
        return res;
    }
}

double StraightStroke2d::area() const
{
    if (size() <= 1)
    {
        return 0;
    }
    else
    {
        double res = 0;
        for (unsigned int i = 1; i < size(); i++)
        {
            res += (points[i].y - points[i-1].y) * (points[i].x + points[i-1].x) / 2;
        }
        return res;
    }
}

std::vector<Point2d> StraightStroke2d::operator/(const SimpleStroke2d& B)
{
    // This may be an infinite loop...
    return *this / B;
}

std::vector<Point2d> StraightStroke2d::operator/(const StraightStroke2d& B)
{
    std::vector<Point2d> res;

    if (size() > 1 && B.size() > 1)
    {
        for (unsigned int i = 1; i < size(); i++)
        {
            // Director vector of the segment.
            Vector2d v = points[i] - points[i-1];
            // Perpendicular vector.
            Vector2d w = v.perp();

            for (unsigned int j = 1; j < B.size(); j++)
            {
                Vector2d u2 = B.points[j] - points[i-1];
                Vector2d u1 = B.points[j-1] - points[i-1];
                double y2 = u2 | w;
                double y1 = u1 | w;
                double x2 = u2 | v;
                double x1 = u1 | v;
                if (y1*y2 <= 0 && x1<=1 && x2<=1)
                {
                    res.push_back(points[i-1] + (x1 - y1*(x2 - x1)/(y2 - y1))*v);
                }
            }
        }
    }
    return res;
}

std::string StraightStroke2d::render_to_string() const
{
    if (points.size() == 0)
    {
        return "";
    }
    else
    {
        std::string text = "\\pgfpathmoveto{\\pgfpoint{"
            + std::to_string(points[0].x) + "cm}{" + std::to_string(points[0].y) + "cm}}\n";
        for (unsigned int i = 0; i < points.size() - 1; i++)
        {
            text += "\\pgfpathlineto{\\pgfpoint{" + std::to_string(points[i].x) + "cm}{"
                + std::to_string(points[i].y) + "cm}}\n";
        }

        // We decide whether we have to close the path.
        if (end() == start())
        {
            text += "\\pgfpathclose\n";
        }
        else
        {
            text += "\\pgfpathlineto{\\pgfpoint{" + std::to_string(points[size() - 1].x) + "cm}{"
                + std::to_string(points[size() - 1].y) + "cm}}\n";
        }
        return text;
    }
}

BezierStroke2d::BezierStroke2d(): BezierStroke2d(Point2d(0,0), Point2d(0,0), Point2d(0,0), Point2d(0,0))
{

}

BezierStroke2d::BezierStroke2d(const Point2d& P1, const Point2d& P2, const Point2d& Q1, const Point2d& Q2):
    P1(P1), P2(P2), Q1(Q1), Q2(Q2)
{

}

BezierStroke2d* BezierStroke2d::clone() const
{
    return new BezierStroke2d(*this);
}

Point2d& BezierStroke2d::start()
{
    return P1;
}

Point2d BezierStroke2d::start() const
{
    return P1;
}

Point2d& BezierStroke2d::end()
{
    return P2;
}

Point2d BezierStroke2d::end() const
{
    return P2;
}

BezierStroke2d& BezierStroke2d::translate(const Vector2d& v)
{
    P1 += v;
    P2 += v;
    Q1 += v;
    Q2 += v;
    return *this;
}

BezierStroke2d& BezierStroke2d::rotate_with_respect_to(const Point2d& Q, const double theta)
{
    P1.rotate_with_respect_to(Q, theta);
    P2.rotate_with_respect_to(Q, theta);
    Q1.rotate_with_respect_to(Q, theta);
    Q2.rotate_with_respect_to(Q, theta);
    return *this;
}

BezierStroke2d& BezierStroke2d::scale_with_respect_to(const Point2d& Q, const Vector2d& s)
{
    P1.scale_with_respect_to(Q, s);
    P2.scale_with_respect_to(Q, s);
    Q1.scale_with_respect_to(Q, s);
    Q2.scale_with_respect_to(Q, s);
    return *this;
}

double BezierStroke2d::length() const
{
    // IS YET TO BE IMPLEMENTED!
    return 0;
}

double BezierStroke2d::area() const
{
    // IS YET TO BE IMPLEMENTED!
    return 0;
}

std::vector<Point2d> BezierStroke2d::operator/(const SimpleStroke2d& B)
{
    // This may be an infinite loop...
    return *this / B;
}

std::string BezierStroke2d::render_to_string() const
{
    return "\\pgfpathmoveto{\\pgfpoint{" + std::to_string(P1.x) + "cm}{" + std::to_string(P1.y) + "cm}}\n" +
        "\\pgfpathcurveto{\\pgfpoint{" + std::to_string(Q1.x) + "cm}{" + std::to_string(Q1.y) + "cm}}{\\pgfpoint{" +
        std::to_string(Q2.x) + "cm}{" + std::to_string(Q2.y) + "cm}}{\\pgfpoint{" +
        std::to_string(P2.x) + "cm}{" + std::to_string(P2.y) + "cm}}\n";
}


#endif // STROKES2D_CPP
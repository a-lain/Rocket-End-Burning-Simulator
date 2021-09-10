#ifndef DATA_PLOT_CPP
#define DATA_PLOT_CPP

#include "DataPlot.hpp"
#include <stdexcept>
#include <limits>
#include <cmath>

using namespace CPGF::AffineSpace;
using namespace CPGF::Plot2d;
using namespace CPGF::Objects2d;
using namespace CPGF;

DataPlot::DataPlot(const std::vector<double>& Y, const std::vector<double>& X,
    const Color&color, const double size,
    const double opacity, std::function<Object2d(
        const Point2d pos,
        const Color& color, const double opacity, const double size)>
            shape,
    const std::string& legend)
{
    if (Y.size() != X.size())
    {
        throw std::runtime_error("Y and X do not have the same size!");
    }
    else
    {
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(X[i], Y[i]));
        }
        this->color = [=](unsigned int i){return color;};
        this->size = [=](unsigned int i){return size;};
        this->opacity = [=](unsigned int i){return opacity;};
        this->shape = [&](unsigned int i){return shape;};
        this->legend = legend;
    }
}

DataPlot::DataPlot(const std::vector<double>& Y, const std::vector<double>& X,
    std::function<Color(unsigned int)> color,
    std::function<double(unsigned int)> size,
    std::function<double(unsigned int)> opacity,
    std::function<std::function<Objects2d::Object2d(
        const AffineSpace::Point2d pos,
        const Color& color, const double opacity, const double size)>
        (unsigned int)> shape,
    const std::string& legend)
{
    if (Y.size() != X.size())
    {
        throw std::runtime_error("Y and X do not have the same size!");
    }
    else
    {
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(X[i], Y[i]));
        }
        this->color = color;
        this->size = size;
        this->opacity = opacity;
        this->shape = shape;
        this->legend = legend;
    }
}

DataPlot::DataPlot(const std::vector<double>& Y, const std::vector<double>& X,
    std::function<Color(AffineSpace::Point2d&)> color,
    std::function<double(AffineSpace::Point2d&)> size,
    std::function<double(AffineSpace::Point2d&)> opacity,
    std::function<std::function<Objects2d::Object2d(
        const AffineSpace::Point2d pos, 
        const Color& color, const double opacity, const double size)>
        (AffineSpace::Point2d&)> shape,
    const std::string& legend)
{
    if (Y.size() != X.size())
    {
        throw std::runtime_error("Y and X do not have the same size!");
    }
    else
    {
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(X[i], Y[i]));
        }
        this->color = [&](unsigned int i){return color(points[i]);};
        this->size = [&](unsigned int i){return size(points[i]);};
        this->opacity = [&](unsigned int i){return opacity(points[i]);};
        this->shape = [&](unsigned int i){return shape(points[i]);};
        this->legend = legend;
    }
}

DataPlot::DataPlot(const std::vector<AffineSpace::Point2d>& data,
    const Color&color, const double size,
    const double opacity, std::function<Objects2d::Object2d(
        const AffineSpace::Point2d pos,const Color& color,
        const double opacity, const double size)> shape,
    const std::string& legend):
        points(data), color([=](unsigned int i){return color;}),
        size([=](unsigned int i){return size;}),
        opacity([=](unsigned int i){return opacity;}),
        shape([=](unsigned int i){return shape;})
{
    this->legend = legend;
}

DataPlot::DataPlot(const std::vector<AffineSpace::Point2d>& data,
    std::function<Color(unsigned int)> color,
    std::function<double(unsigned int)> size,
    std::function<double(unsigned int)> opacity,
    std::function<std::function<Objects2d::Object2d(const AffineSpace::Point2d pos,
        const Color& color, const double opacity, const double size)>
        (unsigned int)> shape,
    const std::string& legend):
        points(data), color(color), size(size), opacity(opacity), shape(shape)
{
    this->legend = legend;
}

DataPlot::DataPlot(const std::vector<AffineSpace::Point2d>& data,
    std::function<Color(AffineSpace::Point2d&)> color,
    std::function<double(AffineSpace::Point2d&)> size,
    std::function<double(AffineSpace::Point2d&)> opacity,
    std::function<std::function<Objects2d::Object2d(const AffineSpace::Point2d pos,
        const Color& color, const double opacity, const double size)>
        (AffineSpace::Point2d&)> shape,
    const std::string& legend)
{
    this->points = data;
    this->color = [=](unsigned int i){return color(points[i]);};
    this->size = [=](unsigned int i){return size(points[i]);};
    this->opacity = [=](unsigned int i){return opacity(points[i]);};
    this->shape = [=](unsigned int i){return shape(points[i]);};
    this->legend = legend;
}

double DataPlot::x_min() const
{
    double res = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < points.size(); i++)
    {
        if (points[i].x < res)
        {
            res = points[i].x;
        }
    }
    return res;
}

double DataPlot::x_max() const
{
    double res = std::numeric_limits<double>::lowest();
    for (unsigned int i = 0; i < points.size(); i++)
    {
        if (points[i].x > res)
        {
            res = points[i].x;
        }
    }
    return res;
}

double DataPlot::y_min() const
{
    double res = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < points.size(); i++)
    {
        if (points[i].y < res)
        {
            res = points[i].y;
        }
    }
    return res;
}

double DataPlot::y_max() const
{
    double res = std::numeric_limits<double>::lowest();
    for (unsigned int i = 0; i < points.size(); i++)
    {
        if (points[i].y > res)
        {
            res = points[i].y;
        }
    }
    return res;
}

Object2d DataPlot::miniature(const AffineSpace::Point2d& pos) const
{
    return Object2d(shape(0)(pos, color(0), opacity(0), size(0)));
}

Scene2d DataPlot::render_to_scene(
    std::function<AffineSpace::Point2d(const AffineSpace::Point2d& P)> transform,
    const double x_min, const double x_max, const double y_min, const double y_max) const
{
    Scene2d scene;
    for (unsigned int i = 0; i < points.size(); i++)
    {
        if (points[i].x >= x_min && points[i].x <= x_max && points[i].y >= y_min && points[i].y <= y_max
            && std::isfinite(points[i].x) && std::isfinite(points[i].y))
        {
            Object2d* object = new Object2d(shape(i)(transform(points[i]), color(i), opacity(i), size(i)));
            scene += *object;
        }
    }

    return scene;
}

Object2d Shapes::Circle(const AffineSpace::Point2d& pos,
    const Color& color, const double opacity,
    const double size)
{
    return Objects2d::Circle(pos, size/10, false, true, color, color, opacity,
        LineWidth::ULTRA_THIN, DashPatterns::SOLID);
}

#endif // DATA_PLOT_CPP
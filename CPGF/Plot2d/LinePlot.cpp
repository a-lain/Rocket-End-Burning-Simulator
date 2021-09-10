#ifndef LINE_PLOT_CPP
#define LINE_PLOT_CPP

#include "LinePlot.hpp"
#include <limits>
#include <stdexcept>
#include <cmath>
#include <iostream>

using namespace CPGF;
using namespace CPGF::Objects2d;
using namespace CPGF::AffineSpace;
using namespace CPGF::Plot2d;

LinePlot::LinePlot(const std::vector<double>& Y, const std::vector<double>& X,
    const Color& color, const double line_width,
    const double opacity, const std::vector<double>& dash_pattern,
    const std::string& legend)
{
    if (Y.size() != X.size())
    {
        throw std::runtime_error("Y does not have the same size as X!");
    }
    else
    {
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(X[i], Y[i]));
        }
        this->color = [=](unsigned int i){return color;};
        this->line_width = [=](unsigned int i){return line_width;};
        this->opacity = [=](unsigned int i){return opacity;};
        this->dash_pattern = [=](unsigned int i){return dash_pattern;};
        this->legend = legend;
        this->const_parameters = true;
    }
}

LinePlot::LinePlot(const std::vector<double>& Y, const std::vector<double>& X,
    std::function<Color(unsigned int)> color, std::function<double(unsigned int)> line_width,
    std::function<double(unsigned int)> opacity,
    std::function<std::vector<double>(unsigned int)> dash_pattern,
    const std::string& legend)
{
    if (Y.size() != X.size())
    {
        throw std::runtime_error("Y does not have the same size as X!");
    }
    else
    {
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(X[i], Y[i]));
        }
        this->color = color;
        this->line_width = line_width;
        this->opacity = opacity;
        this->dash_pattern = dash_pattern;
        this->legend = legend;
        this->const_parameters = false;
    }
}

LinePlot::LinePlot(const std::vector<double>& Y, const std::vector<double>& X,
    std::function<Color(AffineSpace::Point2d&)> color, std::function<double(AffineSpace::Point2d&)> line_width,
    std::function<double(AffineSpace::Point2d&)> opacity,
    std::function<std::vector<double>(AffineSpace::Point2d&)> dash_pattern,
    const std::string& legend)
{
    if (Y.size() != X.size())
    {
        throw std::runtime_error("Y does not have the same size as X!");
    }
    else
    {
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(X[i], Y[i]));
        }
        this->color = [&](unsigned int i){return color(points[i]);};
        this->line_width = [&](unsigned int i){return line_width(points[i]);};
        this->opacity = [&](unsigned int i){return opacity(points[i]);};
        this->dash_pattern = [&](unsigned int i){return dash_pattern(points[i]);};
        this->legend = legend;
        this->const_parameters = false;
    }
}

LinePlot::LinePlot(const std::vector<AffineSpace::Point2d>& data,
    const Color& color, const double line_width,
    const double opacity, const std::vector<double>& dash_pattern,
    const std::string& legend):
        points(data), color([=](unsigned int i){return color;}),
        line_width([=](unsigned int i){return line_width;}),
        opacity([=](unsigned int i){return opacity;}),
        dash_pattern([=](unsigned int i){return dash_pattern;})
{
    this->legend = legend;
    this->const_parameters = false;
}

LinePlot::LinePlot(const std::vector<AffineSpace::Point2d>& data,
    std::function<Color(unsigned int)> color, std::function<double(unsigned int)> line_width,
    std::function<double(unsigned int)> opacity,
    std::function<std::vector<double>(unsigned int)> dash_pattern,
    const std::string& legend):
        points(data), color(color), line_width(line_width), opacity(opacity), dash_pattern(dash_pattern)
{
    this->legend = legend;
    this->const_parameters = false;
}

LinePlot::LinePlot(const std::vector<AffineSpace::Point2d>& data,
    std::function<Color(AffineSpace::Point2d&)> color, std::function<double(AffineSpace::Point2d&)> line_width,
    std::function<double(AffineSpace::Point2d&)> opacity,
    std::function<std::vector<double>(AffineSpace::Point2d&)> dash_pattern,
    const std::string& legend):
        points(data)
{
    this->color = [=](unsigned int i){return color(points[i]);};
    this->line_width = [=](unsigned int i){return line_width(points[i]);};
    this->opacity = [=](unsigned int i){return opacity(points[i]);};
    this->dash_pattern = [=](unsigned int i){return dash_pattern(points[i]);};
    this->legend = legend;
    this->const_parameters = false;
}

double LinePlot::x_min() const
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

double LinePlot::x_max() const
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

double LinePlot::y_min() const
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

double LinePlot::y_max() const
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

Object2d LinePlot::miniature(const AffineSpace::Point2d& pos) const
{
    return Line(pos + Vector2d(-MINIATURE_HALF_WIDTH*0.75, 0), pos + Vector2d(MINIATURE_HALF_WIDTH*0.75, 0),
        color(0), opacity(0), line_width(0), dash_pattern(0));
}

Scene2d LinePlot::render_to_scene(
    std::function<AffineSpace::Point2d(const AffineSpace::Point2d& P)> transform,
    const double x_min, const double x_max, const double y_min, const double y_max) const
{
    Scene2d scene;
    if (!const_parameters)
    {
        for (unsigned int i = 0; i < points.size() - 1; i++)
        {
            // For the moment, we only trace completely contained segments.
            if (points[i].x >= x_min && points[i].x <= x_max && points[i].y >= y_min && points[i].y <= y_max
                && points[i+1].x >= x_min && points[i+1].x <= x_max && points[i+1].y >= y_min && points[i+1].y <= y_max
                && std::isfinite(points[i].x) && std::isfinite(points[i].y)
                && std::isfinite(points[i+1].x) && std::isfinite(points[i+1].y))
            {
                Object2d* object = new Object2d(Line(transform(points[i]), transform(points[i+1]), color(i),
                    opacity(i), line_width(i), dash_pattern(i)));
                scene += *object;
            }
        }   
    }
    else
    {
        std::vector<Point2d> transformed_points;
        for (unsigned int i = 0; i < points.size(); i++)
        {
            if (std::isfinite(points[i].x) && std::isfinite(points[i].y))
            {
                transformed_points.push_back(transform(points[i]));
            }
        }
        Object2d* object = new Object2d(Line(transformed_points, color(0), opacity(0), line_width(0), dash_pattern(0)));
        scene += *object;
    }

    return scene;
}

AveragePlot::AveragePlot(const std::vector<double>& Y, const std::vector<double>& partition,
    const Color& color, const double line_width,
    const double opacity, const std::vector<double>& dash_pattern,
    const std::string& legend) : 
        LinePlot(builder(Y, partition, color, line_width, opacity, dash_pattern, legend))
{

}

AveragePlot::AveragePlot(const std::vector<double>& Y, const std::vector<double>& partition,
    std::function<Color(unsigned int)> color, std::function<double(unsigned int)> line_width,
    std::function<double(unsigned int)> opacity,
    std::function<std::vector<double>(unsigned int)> dash_pattern,
    const std::string& legend):
        LinePlot(builder(Y, partition, color, line_width, opacity, dash_pattern, legend))
{

}

AveragePlot::AveragePlot(const std::vector<double>& Y, const std::vector<double>& partition,
    std::function<Color(AffineSpace::Point2d&)> color, std::function<double(AffineSpace::Point2d&)> line_width,
    std::function<double(AffineSpace::Point2d&)> opacity,
    std::function<std::vector<double>(AffineSpace::Point2d&)> dash_pattern,
    const std::string& legend):
        LinePlot(builder(Y, partition, color, line_width, opacity, dash_pattern, legend))
{

}

LinePlot AveragePlot::builder(const std::vector<double>& Y, const std::vector<double>& partition,
    const Color& color, const double line_width,
    const double opacity, const std::vector<double>& dash_pattern,
    const std::string& legend)
{
    if (partition.size() != Y.size() + 1)
    {
        throw std::runtime_error("The partition must have one more point than the average values.");
    }
    else
    {
        std::vector<Point2d> points;
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(partition[i], Y[i]));
            points.push_back(Point2d(partition[i+1], Y[i]));
        }
        return LinePlot(points, color, line_width, opacity, dash_pattern, legend);
    }
}

LinePlot AveragePlot::builder(const std::vector<double>& Y, const std::vector<double>& partition,
    std::function<Color(unsigned int)> color, std::function<double(unsigned int)> line_width,
    std::function<double(unsigned int)> opacity,
    std::function<std::vector<double>(unsigned int)> dash_pattern,
    const std::string& legend)
{
    if (partition.size() != Y.size() + 1)
    {
        throw std::runtime_error("The partition must have one more point than the average values.");
    }
    else
    {
        std::vector<Point2d> points;
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(partition[i], Y[i]));
            points.push_back(Point2d(partition[i+1], Y[i]));
        }
        return LinePlot(points, color, line_width, opacity, dash_pattern, legend);
    }
}

LinePlot AveragePlot::builder(const std::vector<double>& Y, const std::vector<double>& partition,
    std::function<Color(AffineSpace::Point2d&)> color, std::function<double(AffineSpace::Point2d&)> line_width,
    std::function<double(AffineSpace::Point2d&)> opacity,
    std::function<std::vector<double>(AffineSpace::Point2d&)> dash_pattern,
    const std::string& legend)
{
    if (partition.size() != Y.size() + 1)
    {
        throw std::runtime_error("The partition must have one more point than the average values.");
    }
    else
    {
        std::vector<Point2d> points;
        for (unsigned int i = 0; i < Y.size(); i++)
        {
            points.push_back(Point2d(partition[i], Y[i]));
            points.push_back(Point2d(partition[i+1], Y[i]));
        }
        return LinePlot(points, color, line_width, opacity, dash_pattern, legend);
    }
}

#endif // LINE_PLOT_CPP
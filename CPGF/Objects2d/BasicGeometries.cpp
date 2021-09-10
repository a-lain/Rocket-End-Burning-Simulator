#ifndef BASICGEOMETRIES_CPP
#define BASICGEOMETRIES_CPP

#include "BasicGeometries.hpp"

using namespace CPGF::Objects2d;
using namespace CPGF::AffineSpace;
using namespace CPGF::Basics;

Line::Line(const Point2d& A, const Point2d& B,
    const Color& color, const double opacity,
    const double line_width,
    const std::vector<double>& dash_pattern):
        Line(std::vector<Point2d>({A, B}),
        color, opacity, line_width, dash_pattern)
{

}

Line::Line(const std::vector<Point2d>& points,
    const Color& color, const double opacity,
    const double line_width, const std::vector<double>& dash_pattern):
        Object2d(builder(points, color, opacity, line_width,
        dash_pattern))
{

}

Object2d Line::builder(const std::vector<Point2d>& points,
        const Color& color, const double opacity,
        const double line_width,
        const std::vector<double>& dash_pattern)
{
    PGFConf conf(DrawType::DRAW, color, opacity, line_width,
        dash_pattern, 0, LineCap::BUTT, LineJoin::BEVEL);
    StraightStroke2d stroke(points);
    return Object2d(Path2d(stroke, conf));
}

Point2d& Line::start()
{
    return paths[0].start();
}

Point2d& Line::end()
{
    return paths[0].end();
}

Circle::Circle(const Point2d& pos, const double radius,
    const bool draw, const bool fill, const Color& draw_color,
    const Color& fill_color, const double opacity, const double line_width,
    const std::vector<double>& dash_pattern):
        Object2d(builder(pos, radius, draw, fill, draw_color,
            fill_color, opacity, line_width, dash_pattern))
{

}

Object2d Circle::builder(const Point2d& pos, const double radius,
    const bool draw, const bool fill, const Color& draw_color,
    const Color& fill_color, const double opacity, const double line_width,
    const std::vector<double>& dash_pattern)
{
    double c = 0.551915024494 * radius;
    Point2d* P = new Point2d[4];
    Point2d* Q = new Point2d[8];
    P[0] = pos + Vector2d(radius, 0);
    P[1] = pos + Vector2d(0, radius);
    P[2] = pos + Vector2d(-radius, 0);
    P[3] = pos + Vector2d(0, -radius);
    Q[0] = P[0] + Vector2d(0, c);
    Q[1] = P[1] + Vector2d(c, 0);
    Q[2] = P[1] + Vector2d(-c, 0);
    Q[3] = P[2] + Vector2d(0, c);
    Q[4] = P[2] + Vector2d(0, -c);
    Q[5] = P[3] + Vector2d(-c, 0);
    Q[6] = P[3] + Vector2d(c, 0);
    Q[7] = P[0] + Vector2d(0, -c);

    BezierStroke2d strokes[4];
    strokes[0] = BezierStroke2d(P[0], P[1], Q[0], Q[1]);
    strokes[1] = BezierStroke2d(P[1], P[2], Q[2], Q[3]);
    strokes[2] = BezierStroke2d(P[2], P[3], Q[4], Q[5]);
    strokes[3] = BezierStroke2d(P[3], P[0], Q[6], Q[7]);

    Object2d circle;
    if (fill)
    {
        double epsilon = radius / 1000; 
        PGFConf conf_fill(DrawType::FILL, fill_color, opacity, line_width);
        StraightStroke2d square_stroke(std::vector<Point2d>({P[0] + Vector2d(epsilon, 0),
            P[1] + Vector2d(0, epsilon), P[2] + Vector2d(-epsilon, 0), P[3] + Vector2d(0, -epsilon),
            P[0] + Vector2d(epsilon, 0)}));
        circle += Path2d(std::vector<SimpleStroke2d*>({strokes, strokes + 1, strokes + 2, strokes + 3}), conf_fill);
        circle += Path2d(square_stroke, conf_fill);
    }
    if (draw)
    {
        PGFConf conf_draw(DrawType::DRAW, draw_color, opacity, line_width, dash_pattern, 0,
            LineCap::BUTT, LineJoin::BEVEL);
        circle += Path2d(std::vector<SimpleStroke2d*>({strokes, strokes + 1, strokes + 2, strokes + 3}), conf_draw);
    }
    
    delete[] P;
    delete[] Q;

    return circle;
}

Arrow::Arrow(const AffineSpace::Point2d& start,
    const AffineSpace::Point2d& end,
    const double arrow_head_length, const double arrow_head_width,
    const Color& color, const double opacity, const double line_width,
    const std::vector<double>& dash_pattern):
        Object2d(builder(start, end, arrow_head_length, arrow_head_width,
        color, opacity, line_width, dash_pattern))
{

}

Object2d Arrow::builder(const AffineSpace::Point2d& start,
    const AffineSpace::Point2d& end,
    const double arrow_head_length, const double arrow_head_width,
    const Color& color, const double opacity, const double line_width,
    const std::vector<double>& dash_pattern)
{
    PGFConf conf(DrawType::DRAW, color, opacity, line_width,
        dash_pattern, 0, LineCap::BUTT, LineJoin::MITER);
    StraightStroke2d stroke1(std::vector<Point2d>({start, end}));
    Vector2d v = end - start;
    Vector2d w = v.perp();
    StraightStroke2d stroke2(std::vector<Point2d>({
        end + (-arrow_head_length)*v + arrow_head_width/2*w,
        end, end + (-arrow_head_length)*v + (-arrow_head_width/2)*w}));
    return Object2d(Path2d(std::vector<SimpleStroke2d*>({&stroke1, &stroke2}), conf));
}

#endif // BASICGEOMETRIES_CPP
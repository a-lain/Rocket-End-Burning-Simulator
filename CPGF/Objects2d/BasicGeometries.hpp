#ifndef BASICGEOMETRIES_HPP
#define BASICGEOMETRIES_HPP

#include "Object2d.hpp"
#include "../PGFBasics/Path2d.hpp"

namespace CPGF
{
    namespace Objects2d
    {
        class Line: public Object2d
        {
            public:
            Line(const AffineSpace::Point2d& A, const AffineSpace::Point2d& B,
                const Color& color = Color::BLACK, const double opacity = 1,
                const double line_width = LineWidth::SEMITHICK,
                const std::vector<double>& dash_pattern = DashPatterns::SOLID);
            Line(const std::vector<AffineSpace::Point2d>& points,
                const Color& color = Color::BLACK, const double opacity = 1,
                const double line_width = LineWidth::SEMITHICK,
                const std::vector<double>& dash_pattern = DashPatterns::SOLID);

            AffineSpace::Point2d& start();
            AffineSpace::Point2d& end();

            protected:
            static Object2d builder(const std::vector<AffineSpace::Point2d>& points,
                const Color& color, const double opacity, const double line_width,
                const std::vector<double>& dash_pattern);
        };

        class Circle: public Object2d
        {
            public:
            Circle(const AffineSpace::Point2d& pos, const double radius,
                const bool draw = true, const bool fill = false,
                const Color& draw_color = Color::BLACK,
                const Color& fill_color = Color::WHITE,
                const double opacity = 1,
                const double line_width = LineWidth::SEMITHICK,                
                const std::vector<double>& dash_pattern = DashPatterns::SOLID);

            AffineSpace::Point2d center() const;
            double radius() const;

            protected:
            static Object2d builder(const AffineSpace::Point2d& pos, const double radius,
                const bool draw, const bool fill, const Color& draw_color,
                const Color& fill_color, const double opacity, const double line_width,
                const std::vector<double>& dash_pattern);
        };

        class Arrow: public Object2d
        {
            public:
            Arrow(const AffineSpace::Point2d& start,
                const AffineSpace::Point2d& end,
                const double arrow_head_length = 0.15, const double arrow_head_width = 0.3,
                const Color& color = Color::BLACK, const double opacity = 1,
                const double line_width = LineWidth::SEMITHICK,
                const std::vector<double>& dash_pattern = DashPatterns::SOLID);

            protected:
            static Object2d builder(const AffineSpace::Point2d& start,
                const AffineSpace::Point2d& end,
                const double arrow_head_length, const double arrow_head_width,
                const Color& color, const double opacity, const double line_width,
                const std::vector<double>& dash_pattern);
        };
    }
}

#endif // BASICGEOMETRIES_HPP
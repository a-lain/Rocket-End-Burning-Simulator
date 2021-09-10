#ifndef GRAPHICOBJECT_HPP
#define GRAPHICOBJECT_HPP

#include "../Scene2d.hpp"
#include "../AffineSpace2d/Point2d.hpp"
#include <functional>

namespace CPGF
{
    namespace Plot2d
    {
        class GraphicObject
        {
            public:
            virtual double x_min() const = 0;
            virtual double x_max() const = 0;
            virtual double y_min() const = 0;
            virtual double y_max() const = 0;

            virtual Objects2d::Object2d miniature(const AffineSpace::Point2d& pos) const = 0;

            virtual Scene2d render_to_scene(
                std::function<AffineSpace::Point2d(const AffineSpace::Point2d& P)> transform,
                const double x_min, const double x_max, const double y_min, const double y_max) const = 0;
            std::string legend;

            constexpr static const double MINIATURE_HALF_WIDTH = 0.5;
        };
    }
}

#endif // GRAPHICOBJECT_HPP
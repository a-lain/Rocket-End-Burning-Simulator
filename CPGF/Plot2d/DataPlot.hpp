#ifndef DATA_PLOT_HPP
#define DATA_PLOT_HPP

#include "../AffineSpace2d/Point2d.hpp"
#include "../Objects2d/BasicGeometries.hpp"
#include "GraphicObject.hpp"
#include <functional>

namespace CPGF
{
    namespace Plot2d
    {
        class Shapes
        {
            public:
            static Objects2d::Object2d Circle(const AffineSpace::Point2d& pos,
                const Color& color, const double opacity,
                const double size);
            static Objects2d::Object2d Square(const AffineSpace::Point2d& pos,
                const Color& color, const double opacity,
                const double size);
        };

        class DataPlot: public GraphicObject
        {
            public:
            std::vector<AffineSpace::Point2d> points;
            std::function<Color(unsigned int)> color;
            std::function<double(unsigned int)> size;
            std::function<double(unsigned int)> opacity;
            std::function<std::function<Objects2d::Object2d(const AffineSpace::Point2d pos,
                const Color& color, const double opacity, const double size)>
                (unsigned int)> shape;


            DataPlot(const std::vector<double>& Y, const std::vector<double>& X,
                const Color&color = Color::BLUE, const double size = 1,
                const double opacity = 1, std::function<Objects2d::Object2d(
                    const AffineSpace::Point2d pos,
                    const Color& color, const double opacity, const double size)>
                        shape = Shapes::Circle,
                const std::string& legend = "");

            DataPlot(const std::vector<double>& Y, const std::vector<double>& X,
                std::function<Color(unsigned int)> color,
                std::function<double(unsigned int)> size,
                std::function<double(unsigned int)> opacity,
                std::function<std::function<Objects2d::Object2d(
                    const AffineSpace::Point2d pos,
                    const Color& color, const double opacity, const double size)>
                    (unsigned int)> shape,
                const std::string& legend = "");

            DataPlot(const std::vector<double>& Y, const std::vector<double>& X,
                std::function<Color(AffineSpace::Point2d&)> color,
                std::function<double(AffineSpace::Point2d&)> size,
                std::function<double(AffineSpace::Point2d&)> opacity,
                std::function<std::function<Objects2d::Object2d(
                    const AffineSpace::Point2d pos, 
                    const Color& color, const double opacity, const double size)>
                    (AffineSpace::Point2d&)> shape,
                const std::string& legend = "");

            DataPlot(const std::vector<AffineSpace::Point2d>& data,
                const Color&color = Color::BLUE, const double size = 1,
                const double opacity = 1, std::function<Objects2d::Object2d(
                    const AffineSpace::Point2d pos,const Color& color,
                    const double opacity, const double size)> shape =
                    Shapes::Circle,
                const std::string& legend = "");

            DataPlot(const std::vector<AffineSpace::Point2d>& data,
                std::function<Color(unsigned int)> color,
                std::function<double(unsigned int)> size,
                std::function<double(unsigned int)> opacity,
                std::function<std::function<Objects2d::Object2d(const AffineSpace::Point2d pos,
                    const Color& color, const double opacity, const double size)>
                    (unsigned int)> shape,
                const std::string& legend = "");

            DataPlot(const std::vector<AffineSpace::Point2d>& data,
                std::function<Color(AffineSpace::Point2d&)> color,
                std::function<double(AffineSpace::Point2d&)> size,
                std::function<double(AffineSpace::Point2d&)> opacity,
                std::function<std::function<Objects2d::Object2d(const AffineSpace::Point2d pos,
                    const Color& color, const double opacity, const double size)>
                    (AffineSpace::Point2d&)> shape,
                const std::string& legend = "");

            double x_min() const override;
            double x_max() const override;
            double y_min() const override;
            double y_max() const override;

            Objects2d::Object2d miniature(const AffineSpace::Point2d& pos) const override;

            Scene2d render_to_scene(
                std::function<AffineSpace::Point2d(const AffineSpace::Point2d& P)> transform,
                const double x_min, const double x_max, const double y_min, const double y_max) const override;
        };
    }
}

#endif // DATA_PLOT_HPP
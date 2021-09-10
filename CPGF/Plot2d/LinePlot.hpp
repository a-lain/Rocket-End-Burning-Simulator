#ifndef LINE_PLOT_HPP
#define LINE_PLOT_HPP

#include "../AffineSpace2d/Point2d.hpp"
#include "../Objects2d/BasicGeometries.hpp"
#include "GraphicObject.hpp"
#include <functional>

namespace CPGF
{
    namespace Plot2d
    {
        class LinePlot: public GraphicObject
        {
            public:
            std::vector<AffineSpace::Point2d> points;
            std::function<Color(unsigned int)> color;
            std::function<double(unsigned int)> line_width;
            std::function<double(unsigned int)> opacity;
            std::function<std::vector<double>(unsigned int)> dash_pattern;

            LinePlot(const std::vector<double>& Y, const std::vector<double>& X,
                const Color& color = Color::BLUE, const double line_width = LineWidth::THIN,
                const double opacity = 1, const std::vector<double>& dash_pattern = DashPatterns::SOLID,
                const std::string& legend = "");

            LinePlot(const std::vector<double>& Y, const std::vector<double>& X,
                std::function<Color(unsigned int)> color, std::function<double(unsigned int)> line_width,
                std::function<double(unsigned int)> opacity,
                std::function<std::vector<double>(unsigned int)> dash_pattern,
                const std::string& legend = "");

            LinePlot(const std::vector<double>& Y, const std::vector<double>& X,
                std::function<Color(AffineSpace::Point2d&)> color, std::function<double(AffineSpace::Point2d&)> line_width,
                std::function<double(AffineSpace::Point2d&)> opacity,
                std::function<std::vector<double>(AffineSpace::Point2d&)> dash_pattern,
                const std::string& legend = "");

            LinePlot(const std::vector<AffineSpace::Point2d>& data,
                const Color& color = Color::BLUE, const double line_width = LineWidth::THIN,
                const double opacity = 1, const std::vector<double>& dash_pattern = DashPatterns::SOLID,
                const std::string& legend = "");

            LinePlot(const std::vector<AffineSpace::Point2d>& data,
                std::function<Color(unsigned int)> color, std::function<double(unsigned int)> line_width,
                std::function<double(unsigned int)> opacity,
                std::function<std::vector<double>(unsigned int)> dash_pattern,
                const std::string& legend = "");

            LinePlot(const std::vector<AffineSpace::Point2d>& data,
                std::function<Color(AffineSpace::Point2d&)> color, std::function<double(AffineSpace::Point2d&)> line_width,
                std::function<double(AffineSpace::Point2d&)> opacity,
                std::function<std::vector<double>(AffineSpace::Point2d&)> dash_pattern,
                const std::string& legend = "");

            double x_min() const override;
            double x_max() const override;
            double y_min() const override;
            double y_max() const override;

            Objects2d::Object2d miniature(const AffineSpace::Point2d& pos) const override;

            Scene2d render_to_scene(
                std::function<AffineSpace::Point2d(const AffineSpace::Point2d& P)> transform,
                const double x_min, const double x_max, const double y_min, const double y_max) const override;

            protected:
            bool const_parameters;
        };

        class AveragePlot: public LinePlot
        {
            public:
            /**
             * @brief Allows to average values over a partition of the real line. Returns a LinePlot.
             * 
             * @param Y 
             * @param partition 
             * @param color 
             * @param line_width 
             * @param opacity 
             * @param dash_pattern 
             * @param legend 
             * @return LinePlot 
             */
            AveragePlot(const std::vector<double>& Y, const std::vector<double>& partition,
                const Color& color = Color::BLUE, const double line_width = LineWidth::THIN,
                const double opacity = 1, const std::vector<double>& dash_pattern = DashPatterns::SOLID,
                const std::string& legend = "");

            AveragePlot(const std::vector<double>& Y, const std::vector<double>& partition,
                std::function<Color(unsigned int)> color, std::function<double(unsigned int)> line_width,
                std::function<double(unsigned int)> opacity,
                std::function<std::vector<double>(unsigned int)> dash_pattern,
                const std::string& legend = "");

            AveragePlot(const std::vector<double>& Y, const std::vector<double>& partition,
                std::function<Color(AffineSpace::Point2d&)> color, std::function<double(AffineSpace::Point2d&)> line_width,
                std::function<double(AffineSpace::Point2d&)> opacity,
                std::function<std::vector<double>(AffineSpace::Point2d&)> dash_pattern,
                const std::string& legend = "");

            protected:
            static LinePlot builder(const std::vector<double>& Y, const std::vector<double>& partition,
                const Color& color, const double line_width,
                const double opacity, const std::vector<double>& dash_pattern,
                const std::string& legend);

            static LinePlot builder(const std::vector<double>& Y, const std::vector<double>& partition,
                std::function<Color(unsigned int)> color, std::function<double(unsigned int)> line_width,
                std::function<double(unsigned int)> opacity,
                std::function<std::vector<double>(unsigned int)> dash_pattern,
                const std::string& legend);
                
            static LinePlot builder(const std::vector<double>& Y, const std::vector<double>& partition,
                std::function<Color(AffineSpace::Point2d&)> color, std::function<double(AffineSpace::Point2d&)> line_width,
                std::function<double(AffineSpace::Point2d&)> opacity,
                std::function<std::vector<double>(AffineSpace::Point2d&)> dash_pattern,
                const std::string& legend);
        };
    } 
}

#endif // LINE_PLOT_HPP
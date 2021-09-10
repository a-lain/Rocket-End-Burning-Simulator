#ifndef AXIS_HPP
#define AXIS_HPP

#include "../PGFBasics/PGFConf.hpp"
#include "../Scene2d.hpp"
#include <functional>
#include "../../Utilities/FormatNumber.hpp"
#include "../Objects2d/BasicGeometries.hpp"
#include <cmath>

namespace CPGF
{
    namespace Plot2d
    {
        enum class AxisType
        {
            HORIZONTAL,
            VERTICAL,
        };

        enum class AxisScale
        {
            LINEAR,
            LOG,
        };

        /**
         * @brief This class is used to determine where the numbers are displayed
         * on the axis.
         * 
         * The position is given by where the numbers stand with respect to the axis
         * when we make the axis growth direction be the positive vertical.
         * 
         */
        enum class NumberPosition
        {
            LEFT,
            RIGHT,
        };

        class Axis
        {
            public:

            /**
             * @brief Determines whether the axis is horizontal or vertical.
             * 
             */
            AxisType axis_type;

            /**
             * @brief Represents the scale of the axis.
             * 
             */
            AxisScale axis_scale;

            /**
             * @brief Determines the relative size between the graph and the labels.
             * Use scale=1 for full-paged graphics and scale=0.5 for half-paged graphics.
             * 
             */
            double scale;

            /**
             * @brief Determines whether the axis is rendered or not.
             * 
             */
            bool visible;

            /**
             * @brief Controlls whether the axis is inverted, a.k.a. points in the opposite direction.
             * 
             */
            bool inverted;

            /**
             * @brief The axis label.
             * 
             */
            std::string label;

            /**
             * @brief A real number between zero and one which represents the position of the axis on the graph.
             * 
             */
            double position;

            /**
             * @brief Color of the axis.
             * 
             */
            Color color;

            /**
             * @brief Axis line width.
             * 
             */
            double line_width;

            /**
             * @brief Axis opacity.
             * 
             */
            double opacity;

            /**
             * @brief Axis dash pattern.
             * 
             */
            std::vector<double> dash_pattern;

            /**
             * @brief Determines the length of the axis if it is vertical. To obtain
             * the axis length, we have to multiply the aspect_ratio by the default length
             * of all X axes.
             * 
             * Default value is 1/sqrt(2).
             * 
             */
            double aspect_ratio;

            /**
             * @brief Arrow head length.
             * 
             */
            double arrow_head_length;

            /**
             * @brief Arrow width.
             * 
             */
            double arrow_head_width;

            /**
             * @brief Determines how much the arrow sticks out of the axis.
             * 
             */
            double arrow_length;

            /**
             * @brief Function that creates the arrow head.
             * 
             */
            std::function<Objects2d::Object2d(const AffineSpace::Point2d& start,
                const AffineSpace::Point2d& end,
                const double arrow_head_length, const double arrow_head_width,
                const Color& color, const double opacity, const double line_width,
                const std::vector<double> dash_pattern)> arrow;

            /**
             * @brief Number of major ticks.
             * 
             * Set to 0 to disable.
             * 
             */
            unsigned int N_major_ticks;

            /**
             * @brief Controlls how large the major ticks are in the perpendicular
             * direction to the axis.
             * 
             */
            double major_tick_deviation;

            /**
             * @brief The width of the major ticks will be the axis witdh divided by
             * this factor.
             * 
             */
            double major_tick_line_width_divisor;

            /**
             * @brief Controlls whether medium ticks are displayed.
             * 
             */
            bool show_medium_ticks;

            /**
             * @brief Controlls how large the medium ticks are in the perpendicular
             * direction to the axis.
             * 
             */
            double medium_tick_deviation;

            /**
             * @brief The width of the medium ticks will be the axis witdh divided by
             * this factor.
             * 
             */
            double medium_tick_line_width_divisor;

            /**
             * @brief Controlls whether small ticks are displayed.
             * 
             */
            bool show_small_ticks;

            /**
             * @brief Controlls how large the small ticks are in the perpendicular
             * direction to the axis.
             * 
             */
            double small_tick_deviation;

            /**
             * @brief The width of the small ticks will be the axis witdh divided by
             * this factor.
             * 
             */
            double small_tick_line_width_divisor;

            /**
             * @brief Controlls whether the number scale is displayed.
             * 
             */
            bool show_numbers;

            /**
             * @brief Determines the position of the numbers.
             * 
             */
            NumberPosition number_position;

            /**
             * @brief Determines whether the major grid lines are drawn.
             * 
             */
            bool show_major_grid_lines;

            /**
             * @brief The dash pattern of the major grid lines.
             * 
             */
            std::vector<double> major_grid_lines_dash_pattern;

            /**
             * @brief The width of the major grid lines will be the axis
             * width divided by this factor.
             * 
             */
            double major_grid_line_line_width_divisor;

            /**
             * @brief The opacity of the major grid lines.
             * 
             */
            double major_grid_line_opacity;

            /**
             * @brief Determines whether the medium grid lines are drawn.
             * 
             */
            bool show_medium_grid_lines;

            /**
             * @brief The dash pattern of the medium grid lines.
             * 
             */
            std::vector<double> medium_grid_lines_dash_pattern;

            /**
             * @brief The width of the medium grid lines will be the axis
             * width divided by this factor.
             * 
             */
            double medium_grid_line_line_width_divisor;

            /**
             * @brief The opacity of the medium grid lines.
             * 
             */
            double medium_grid_line_opacity;

            /**
             * @brief Determines whether the small grid lines are drawn.
             * 
             */
            bool show_small_grid_lines;

            /**
             * @brief The dash pattern of the medium grid lines.
             * 
             */
            std::vector<double> small_grid_lines_dash_pattern;

            /**
             * @brief The width of the small grid lines will be the axis
             * width divided by this factor.
             * 
             */
            double small_grid_line_line_width_divisor;

            /**
             * @brief The opacity of the small grid lines.
             * 
             */
            double small_grid_line_opacity;

            void set_max_value(const double value);
            void reset_max_value();
            void set_min_value(const double value);
            void reset_min_value();

            void update_max_min(const double max_value, const double min_value);

            void calculate_transformations();

            double get_min_value() const;
            double get_max_value() const;

            /**
             * @brief Applied to a coordiante, it returns the coordinate on the final drawing.
             * 
             * @return std::function<double(double)> 
             */
            std::function<double(double)> axis_transform() const;

            std::function<double(double)> axis_transform_inverse() const;

            Scene2d render_to_scene() const;


            Axis(const AxisType axis_type, const std::string& label = "",
                const AxisScale axis_scale = AxisScale::LINEAR,
                const double scale = 1,
                const bool visible = true, const bool inverted = false,
                const double position = 0, const Color& color = Color::BLACK,
                const double line_width = LineWidth::THIN,
                const double opacity = 1, const std::vector<double>& dash_pattern = DashPatterns::SOLID,
                const double aspect_ratio = 0.6,
                const double arrow_head_length = 0.3,
                const double arrow_head_width = 0.75,
                const double arrow_length = 0.5,
                std::function<Objects2d::Object2d(const AffineSpace::Point2d& start,
                const AffineSpace::Point2d& end,
                const double arrow_head_length, const double arrow_head_width,
                const Color& color, const double opacity, const double line_width,
                const std::vector<double>& dash_pattern)> arrow = 
                    [] (const AffineSpace::Point2d& start,
                        const AffineSpace::Point2d& end,
                        const double arrow_head_length, const double arrow_head_width,
                        const Color& color, const double opacity, const double line_width,
                        const std::vector<double>& dash_pattern)
                        {return Objects2d::Arrow(start, end, arrow_head_length, arrow_head_width, color,
                        opacity, line_width, dash_pattern);},
                const unsigned int N_major_ticks = 9,
                const double major_tick_deviation = 0.25,
                const double major_tick_line_width_divisor = 3,
                const bool show_medium_ticks = true,
                const double medium_tick_deviation = 0.20,
                const double medium_tick_line_width_divisor = 4,
                const bool show_small_ticks = true,
                const double small_tick_deviation = 0.15,
                const double small_tick_line_width_divisor = 5,
                const bool show_numbers = true,
                const bool show_major_grid_lines = false,
                const std::vector<double>& major_grid_lines_dash_pattern = DashPatterns::SOLID,
                const double major_grid_line_line_width_divisor = 2,
                const double major_grid_line_opacity = 0.75,
                const bool show_medium_grid_lines = false,
                const std::vector<double>& medium_grid_lines_dash_pattern = DashPatterns::SOLID,
                const double medium_grid_line_line_width_divisor = 4,
                const double medium_grid_line_opacity = 0.5,
                const bool show_small_grid_lines = false,
                const std::vector<double>& small_grid_lines_dash_pattern = DashPatterns::SOLID,
                const double small_grid_line_line_width_divisor = 8,
                const double small_grid_line_opacity = 0.25);

            protected:
            bool user_defined_max_value;
            double max_value;
            bool user_defined_min_value;
            double min_value;

            double x_max;
            
            // Transformation parameters.
            double _a;
            double _b;

            int digit_max;
            unsigned int precision;

            public:
            // Predefined constants.

            /**
             * @brief End position of the axis (if it is horizontal) on the final drawing.
             * 
             */
            constexpr static const double X_MAX = 15;

            /**
             * @brief Start position of the axis (if it is horizontal) on the final drawing.
             * 
             */
            constexpr static const double X_MIN = 0;  

            /**
             * @brief Start position of the axis (if it is vertical) on the final drawing.
             * 
             */
            constexpr static const double Y_MIN = 0;

            /**
             * @brief Determines how much the numbers are displaced with respect to the axis.
             * 
             */
            constexpr static const double NUMBER_DISPLACEMENT = 0.5;

            /**
             * @brief Determines how much the label is displaced with respect to the axis.
             * 
             */
            constexpr static const double LABEL_DISPLACEMENT_HORIZONTAL = 1.2;

            constexpr static const double LABEL_DISPLACEMENT_VERTICAL = 2;
        };
    }
}
#endif // AXIS_HPP
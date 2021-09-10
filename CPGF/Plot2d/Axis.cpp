#ifndef AXIS_CPP
#define AXIS_CPP

#include "Axis.hpp"
#include "../Objects2d/BasicGeometries.hpp"
#include "../Text/Text.hpp"
#include <cmath>
#include <limits>

#include <iostream>

using namespace CPGF::Plot2d;
using namespace CPGF::Objects2d;
using namespace CPGF::AffineSpace;
using namespace CPGF;

Axis::Axis(const AxisType axis_type, const std::string& label,
    const AxisScale axis_scale,
    const double scale,
    const bool visible, const bool inverted,
    const double position, const Color& color,
    const double line_width,
    const double opacity, const std::vector<double>& dash_pattern,
    const double aspect_ratio,
    const double arrow_head_length,
    const double arrow_head_width,
    const double arrow_length,
    std::function<Objects2d::Object2d(const AffineSpace::Point2d& start,
    const AffineSpace::Point2d& end,
    const double arrow_head_length, const double arrow_head_width,
    const Color& color, const double opacity, const double line_width,
    const std::vector<double>& dash_pattern)> arrow,
    const unsigned int N_major_ticks,
    const double major_tick_deviation,
    const double major_tick_line_width_divisor,
    const bool show_medium_ticks,
    const double medium_tick_deviation,
    const double medium_tick_line_width_divisor,
    const bool show_small_ticks,
    const double small_tick_deviation,
    const double small_tick_line_width_divisor,
    const bool show_numbers,
    const bool show_major_grid_lines,
    const std::vector<double>& major_grid_lines_dash_pattern,
    const double major_grid_line_line_width_divisor,
    const double major_grid_line_opacity,
    const bool show_medium_grid_lines,
    const std::vector<double>& medium_grid_lines_dash_pattern,
    const double medium_grid_line_line_width_divisor,
    const double medium_grid_line_opacity,
    const bool show_small_grid_lines,
    const std::vector<double>& small_grid_lines_dash_pattern,
    const double small_grid_line_line_width_divisor,
    const double small_grid_line_opacity):
        axis_type(axis_type), axis_scale(axis_scale),
        scale(scale),
        visible(visible), inverted(inverted), label(label),
        position(position), color(color), line_width(line_width),
        opacity(opacity), dash_pattern(dash_pattern), aspect_ratio(aspect_ratio),
        arrow_head_length(arrow_head_length), arrow_head_width(arrow_head_width),
        arrow_length(arrow_length), arrow(arrow), N_major_ticks(N_major_ticks),
        major_tick_deviation(major_tick_deviation),
        major_tick_line_width_divisor(major_tick_line_width_divisor),
        show_medium_ticks(show_medium_ticks),
        medium_tick_deviation(medium_tick_deviation),
        medium_tick_line_width_divisor(medium_tick_line_width_divisor),
        show_small_ticks(show_small_ticks),
        small_tick_deviation(small_tick_deviation),
        small_tick_line_width_divisor(small_tick_line_width_divisor),
        show_numbers(show_numbers),
        show_major_grid_lines(show_major_grid_lines),
        major_grid_lines_dash_pattern(major_grid_lines_dash_pattern),
        major_grid_line_line_width_divisor(major_grid_line_line_width_divisor),
        major_grid_line_opacity(major_grid_line_opacity),
        show_medium_grid_lines(show_medium_grid_lines),
        medium_grid_lines_dash_pattern(medium_grid_lines_dash_pattern),
        medium_grid_line_line_width_divisor(medium_grid_line_line_width_divisor),
        medium_grid_line_opacity(medium_grid_line_opacity),
        show_small_grid_lines(show_small_grid_lines),
        small_grid_lines_dash_pattern(small_grid_lines_dash_pattern),
        small_grid_line_line_width_divisor(small_grid_line_line_width_divisor),
        small_grid_line_opacity(small_grid_line_opacity),
        user_defined_max_value(false),
        max_value(std::numeric_limits<double>::lowest()),
        user_defined_min_value(false),
        min_value(std::numeric_limits<double>::max())
{
    if (axis_type == AxisType::VERTICAL)
    {
        this->number_position = NumberPosition::LEFT;
    }
    else
    {
        this->number_position = NumberPosition::RIGHT;
    }
}

void Axis::set_max_value(const double value)
{
    user_defined_max_value = true;
    max_value = value;
}

void Axis::reset_max_value()
{
    user_defined_max_value = false;
}

void Axis::set_min_value(const double value)
{
    user_defined_min_value = true;
    min_value = value;
}

void Axis::reset_min_value()
{
    user_defined_min_value = false;
}

void Axis::update_max_min(const double max_value, const double min_value)
{
    if (!user_defined_max_value)
    {
        if (max_value > this->max_value)
        {
            this->max_value = max_value;
        }
    }
    if (!user_defined_min_value)
    {
        if (min_value < this->min_value)
        {
            this->min_value = min_value;
        }
    }
}

double Axis::get_max_value() const
{
    return max_value;
}

double Axis::get_min_value() const
{
    return min_value;
}

void Axis::calculate_transformations()
{
    x_max = X_MAX*scale;
    if (axis_scale == AxisScale::LINEAR)
    {
        double delta_val = (max_value - min_value) / (N_major_ticks - 1);
        digit_max = Utilities::position_of_most_significant_digit(fmax(fabs(min_value), fabs(max_value)));
        precision = (unsigned int) (digit_max
            - Utilities::position_of_most_significant_digit(delta_val) + 2);
        max_value = Utilities::ceil_to_precision(max_value,
            std::max(0, Utilities::position_of_most_significant_digit(max_value) - digit_max + (int)precision));
        min_value = Utilities::floor_to_precision(min_value,
            std::max(0, Utilities::position_of_most_significant_digit(min_value) - digit_max + (int)precision));

        if (axis_type == AxisType::HORIZONTAL)
        {
            if (!inverted)
            {
                _a = (x_max - X_MIN) / (max_value - min_value);
                _b = X_MIN - _a*min_value;
            }
            else
            {
                _a = (X_MIN - x_max) / (max_value - min_value);
                _b = X_MIN - _a*max_value;
            }
        }
        else // axis_type = AxisType::VERTICAL
        {
            if (!inverted)
            {
                _a = (x_max*aspect_ratio - Y_MIN) / (max_value - min_value);
                _b = Y_MIN - _a*min_value;
            }
            else
            {
                _a = (Y_MIN - x_max*aspect_ratio) / (max_value - min_value);
                _b = Y_MIN - _a*max_value;
            }
        }
    }
    else // axis_scale = AxisScale::LOG
    {
        unsigned int power_max = (unsigned int)ceil(log10(max_value/min_value));
        if (axis_type == AxisType::HORIZONTAL)
        {
            if (!inverted)
            {
                _a = (x_max - X_MIN) / power_max;
                _b = X_MIN;
            }
            else
            {
                _a = (X_MIN - x_max) / power_max;
                _b = x_max;
            }
        }
        else // axis_type = AxisType::VERTICAL
        {
            if (!inverted)
            {
                _a = (x_max*aspect_ratio - Y_MIN) / power_max;
                _b = Y_MIN;
            }
            else
            {
                _a = (Y_MIN - x_max*aspect_ratio) / power_max;
                _b = x_max*aspect_ratio;
            }
        }
    }
}

std::function<double(double)> Axis::axis_transform() const
{
    if (axis_scale == AxisScale::LINEAR)
    {
        return [=](double x){return _a*x + _b;};
    }
    else // axis_scale = AxisScale::LOG
    {
        return [=](double x){return _a*log10(x / min_value) + _b;};
    }
}

std::function<double(double)> Axis::axis_transform_inverse() const
{
    if (axis_scale == AxisScale::LINEAR)
    {
        return [=](double y){return (y - _b)/_a;};
    }
    else // axis_scale = AxisScale::LOG
    {
        return [=](double y){return min_value*pow(10, (y-_b)/_a);};
    }
}

Scene2d Axis::render_to_scene() const
{
    Scene2d scene;
    Object2d* object;
    Text* text;
    std::function<double(double)> transform = axis_transform();
    std::function<double(double)> transform_inverse = axis_transform_inverse();

    double y_max = (x_max - X_MIN)*aspect_ratio + Y_MIN;
 
    if (visible)
    {
        if (axis_type == AxisType::HORIZONTAL)
        {
            // The axis vertical position.
            double y = (y_max - Y_MIN)*position + Y_MIN;

            // We draw the axis.
            object = new Line(Point2d(X_MIN, y), Point2d(x_max, y),
                color, opacity, line_width, dash_pattern);
            scene += *object;

            // We draw the arrow
            if (!inverted)
            {
                object = new Object2d(arrow(Point2d(x_max, y), Point2d(x_max + arrow_length, y),
                    arrow_head_length, arrow_head_width,
                    color, opacity, line_width, dash_pattern));
                scene += *object;
            }
            else
            {
                object = new Object2d(arrow(Point2d(X_MIN, y), Point2d(X_MIN - arrow_length, y),
                    arrow_head_length, arrow_head_width,
                    color, opacity, line_width, dash_pattern));
                scene += *object;
            }
            
            // We draw the ticks.
            if (N_major_ticks != 0)
            {
                std::vector<double> major_ticks;
                std::vector<double> medium_ticks;
                std::vector<double> small_ticks;

                // We calculate all ticks.
                if (axis_scale == AxisScale::LINEAR)
                {
                    double x = (inverted) ? x_max : X_MIN;
                    double delta_x = (x_max - X_MIN) / (N_major_ticks - 1)
                        * ((inverted) ? -1 : 1);

                    for (unsigned int i = 0; i < N_major_ticks - 1; i++)
                    {
                        major_ticks.push_back(x);
                        medium_ticks.push_back(x + delta_x/2);
                        for (int j = 1; j <= 4; j++)
                        {
                            small_ticks.push_back(x + j*delta_x/10);
                        }
                        for (int j = 6; j <= 9; j++)
                        {
                            small_ticks.push_back(x + j*delta_x/10);
                        }
                        x += delta_x;
                    }

                    major_ticks.push_back(x);
                }
                else // axis_scale = AxisScale::LOG
                {
                    unsigned int power_max = (unsigned int)ceil(log10(max_value/min_value));
                    unsigned int delta_alpha = power_max / N_major_ticks + 1;

                    for (unsigned int alpha = 0; alpha <= power_max; alpha += delta_alpha)
                    {
                        major_ticks.push_back(_a*alpha + _b);
                        if (alpha + delta_alpha - 1 + log10(5) <= power_max)
                        {
                            medium_ticks.push_back(_a*(alpha + delta_alpha - 1 + log10(5)) + _b);
                        }                        
                        for (int j = 2; j <= 4; j++)
                        {
                            if (alpha + delta_alpha - 1 + log10(j) <= power_max)
                            {
                                small_ticks.push_back(_a*(alpha + delta_alpha - 1 + log10(j)) + _b);
                            }
                        }
                        for (int j = 6; j <= 9; j++)
                        {
                            if (alpha + delta_alpha - 1 + log10(j) <= power_max)
                            {
                                small_ticks.push_back(_a*(alpha + delta_alpha - 1 + log10(j)) + _b);
                            }
                        }
                    }
                }

                // We draw the major ticks.
                for (unsigned int i = 0; i < major_ticks.size(); i++)
                {
                    object = new Line(Point2d(major_ticks[i], -major_tick_deviation),
                        Point2d(major_ticks[i], major_tick_deviation), color, opacity,
                        line_width / major_tick_line_width_divisor,
                        DashPatterns::SOLID);
                    scene += *object;
                }

                // We draw the medium ticks.
                if (show_medium_ticks)
                {
                    for (unsigned int i = 0; i < medium_ticks.size(); i++)
                    {
                        object = new Line(Point2d(medium_ticks[i], -medium_tick_deviation),
                            Point2d(medium_ticks[i], medium_tick_deviation), color, opacity,
                            line_width / medium_tick_line_width_divisor,
                            DashPatterns::SOLID);
                        scene += *object;
                    }
                }

                // We draw the small ticks.
                if (show_small_ticks)
                {
                    for (unsigned int i = 0; i < small_ticks.size(); i++)
                    {
                        object = new Line(Point2d(small_ticks[i], -small_tick_deviation),
                            Point2d(small_ticks[i], small_tick_deviation), color, opacity,
                            line_width / small_tick_line_width_divisor,
                            DashPatterns::SOLID);
                        scene += *object;
                    }
                }

                // We draw the major grid lines.
                if (show_major_grid_lines)
                {
                    for (unsigned int i = 0; i < major_ticks.size(); i++)
                    {
                        object = new Line(Point2d(major_ticks[i], Y_MIN),
                            Point2d(major_ticks[i], y_max), color,
                            major_grid_line_opacity,
                            line_width / major_grid_line_line_width_divisor,
                            major_grid_lines_dash_pattern);
                        scene += *object;
                    }
                }

                // We draw the medium grid lines.
                if (show_medium_grid_lines)
                {
                    for (unsigned int i = 0; i < medium_ticks.size(); i++)
                    {
                        object = new Line(Point2d(medium_ticks[i], Y_MIN),
                            Point2d(medium_ticks[i], y_max), color,
                            medium_grid_line_opacity,
                            line_width / medium_grid_line_line_width_divisor,
                            medium_grid_lines_dash_pattern);
                        scene += *object;
                    }
                }

                // We draw the small grid lines.
                if (show_small_grid_lines)
                {
                    for (unsigned int i = 0; i < small_ticks.size(); i++)
                    {
                        object = new Line(Point2d(small_ticks[i], Y_MIN),
                            Point2d(small_ticks[i], y_max), color,
                            small_grid_line_opacity,
                            line_width / small_grid_line_line_width_divisor,
                            small_grid_lines_dash_pattern);
                        scene += *object;
                    }
                }

                // We place the text corresponding to each tick.
                if (show_numbers)
                {
                    double y_text = y + ((number_position == NumberPosition::RIGHT) ?
                        -NUMBER_DISPLACEMENT : NUMBER_DISPLACEMENT);
                    TextAlignment text_alignment = (number_position == NumberPosition::RIGHT) ?
                        TextAlignment::TOP : TextAlignment::BOTTOM;
                    
                    for (unsigned int i = 0; i < major_ticks.size(); i++)
                    {
                        text = new Text(Point2d(major_ticks[i], y_text),
                            "$" + Utilities::format_number(transform_inverse(major_ticks[i]), true,
                                (axis_scale == AxisScale::LINEAR) ?
                                std::max(0, Utilities::position_of_most_significant_digit(transform_inverse(major_ticks[i])) - digit_max + (int)precision)
                                : 3) + "$", color, text_alignment, 0);
                        scene += *text;
                    }
                }

                // We write the axis label.
                text = new Text(Point2d((x_max - X_MIN) / 2 + X_MIN, y + ((number_position == NumberPosition::RIGHT) ?
                    -LABEL_DISPLACEMENT_HORIZONTAL : LABEL_DISPLACEMENT_HORIZONTAL)),
                    "\\Large" + label, color, TextAlignment::CENTER, 0);
                scene += *text;
            }
        }
        else // axis_type == AxisType::VERTICAL
        {
            double x = X_MIN + (x_max - X_MIN)*position;

            // We draw the axis.
            object = new Line(Point2d(x, Y_MIN), Point2d(x, y_max),
                color, opacity, line_width, dash_pattern);
            scene += *object;

            // We draw the arrow
            if (!inverted)
            {
                object = new Object2d(arrow(Point2d(x, y_max), Point2d(x, y_max + arrow_length),
                    arrow_head_length, arrow_head_width,
                    color, opacity, line_width, dash_pattern));
                scene += *object;
            }
            else
            {
                object = new Object2d(arrow(Point2d(x, Y_MIN), Point2d(x, Y_MIN - arrow_length),
                    arrow_head_length, arrow_head_width,
                    color, opacity, line_width, dash_pattern));
                scene += *object;
            }
            
            // We draw the ticks.
            if (N_major_ticks != 0)
            {
                std::vector<double> major_ticks;
                std::vector<double> medium_ticks;
                std::vector<double> small_ticks;

                // We calculate all ticks.
                if (axis_scale == AxisScale::LINEAR)
                {
                    double y = (inverted) ? y_max : Y_MIN;
                    double delta_y = (y_max - Y_MIN) / (N_major_ticks - 1)
                        * ((inverted) ? -1 : 1);

                    for (unsigned int i = 0; i < N_major_ticks - 1; i++)
                    {
                        major_ticks.push_back(y);
                        medium_ticks.push_back(y + delta_y/2);
                        for (int j = 1; j <= 4; j++)
                        {
                            small_ticks.push_back(y + j*delta_y/10);
                        }
                        for (int j = 6; j <= 9; j++)
                        {
                            small_ticks.push_back(y + j*delta_y/10);
                        }
                        y += delta_y;
                    }

                    major_ticks.push_back(y);
                }
                else // axis_scale = AxisScale::LOG
                {
                    unsigned int power_max = (unsigned int)ceil(log10(max_value/min_value));
                    unsigned int delta_alpha = power_max / N_major_ticks + 1;

                    for (unsigned int alpha = 0; alpha <= power_max; alpha += delta_alpha)
                    {
                        major_ticks.push_back(_a*alpha + _b);
                        if (alpha + delta_alpha - 1 + log10(5) <= power_max)
                        {
                            medium_ticks.push_back(_a*(alpha + delta_alpha - 1 + log10(5)) + _b);
                        }                        
                        for (int j = 2; j <= 4; j++)
                        {
                            if (alpha + delta_alpha - 1 + log10(j) <= power_max)
                            {
                                small_ticks.push_back(_a*(alpha + delta_alpha - 1 + log10(j)) + _b);
                            }
                        }
                        for (int j = 6; j <= 9; j++)
                        {
                            if (alpha + delta_alpha - 1 + log10(j) <= power_max)
                            {
                                small_ticks.push_back(_a*(alpha + delta_alpha - 1 + log10(j)) + _b);
                            }
                        }
                    }
                }

                // We draw the major ticks.
                for (unsigned int i = 0; i < major_ticks.size(); i++)
                {
                    object = new Line(Point2d(-major_tick_deviation, major_ticks[i]),
                        Point2d(major_tick_deviation, major_ticks[i]), color, opacity,
                        line_width / major_tick_line_width_divisor,
                        DashPatterns::SOLID);
                    scene += *object;
                }

                // We draw the medium ticks.
                if (show_medium_ticks)
                {
                    for (unsigned int i = 0; i < medium_ticks.size(); i++)
                    {
                        object = new Line(Point2d(-medium_tick_deviation, medium_ticks[i]),
                            Point2d(medium_tick_deviation, medium_ticks[i]), color, opacity,
                            line_width / medium_tick_line_width_divisor,
                            DashPatterns::SOLID);
                        scene += *object;
                    }
                }

                // We draw the small ticks.
                if (show_small_ticks)
                {
                    for (unsigned int i = 0; i < small_ticks.size(); i++)
                    {
                        object = new Line(Point2d(-small_tick_deviation, small_ticks[i]),
                            Point2d(small_tick_deviation, small_ticks[i]), color, opacity,
                            line_width / small_tick_line_width_divisor,
                            DashPatterns::SOLID);
                        scene += *object;
                    }
                }

                // We draw the major grid lines.
                if (show_major_grid_lines)
                {
                    for (unsigned int i = 0; i < major_ticks.size(); i++)
                    {
                        object = new Line(Point2d(X_MIN, major_ticks[i]),
                            Point2d(x_max, major_ticks[i]), color,
                            major_grid_line_opacity,
                            line_width / major_grid_line_line_width_divisor,
                            major_grid_lines_dash_pattern);
                        scene += *object;
                    }
                }

                // We draw the medium grid lines.
                if (show_medium_grid_lines)
                {
                    for (unsigned int i = 0; i < medium_ticks.size(); i++)
                    {
                        object = new Line(Point2d(X_MIN, medium_ticks[i]),
                            Point2d(x_max, medium_ticks[i]), color,
                            medium_grid_line_opacity,
                            line_width / medium_grid_line_line_width_divisor,
                            medium_grid_lines_dash_pattern);
                        scene += *object;
                    }
                }

                // We draw the small grid lines.
                if (show_small_grid_lines)
                {
                    for (unsigned int i = 0; i < small_ticks.size(); i++)
                    {
                        object = new Line(Point2d(X_MIN, small_ticks[i]),
                            Point2d(x_max, small_ticks[i]), color,
                            small_grid_line_opacity,
                            line_width / small_grid_line_line_width_divisor,
                            small_grid_lines_dash_pattern);
                        scene += *object;
                    }
                }

                // We place the text corresponding to each tick.
                if (show_numbers)
                {
                    double x_text = x + ((number_position == NumberPosition::RIGHT) ?
                        NUMBER_DISPLACEMENT : -NUMBER_DISPLACEMENT);
                    TextAlignment text_alignment = (number_position == NumberPosition::RIGHT) ?
                        TextAlignment::LEFT : TextAlignment::RIGHT;

                    for (unsigned int i = 0; i < major_ticks.size(); i++)
                    {
                        text = new Text(Point2d(x_text, major_ticks[i]),
                            "$" + Utilities::format_number(transform_inverse(major_ticks[i]), true,
                            (axis_scale == AxisScale::LINEAR) ?
                            std::max(0, Utilities::position_of_most_significant_digit(transform_inverse(major_ticks[i])) - digit_max + (int)precision)
                            : 3) + "$", color, text_alignment, 0);
                        scene += *text;
                    }
                }
            }

            // We write the axis label.
            double rot = (number_position == NumberPosition::LEFT) ? 90 : -90;         
            text = new Text(Point2d(x + ((number_position == NumberPosition::RIGHT) ?
                LABEL_DISPLACEMENT_VERTICAL : -LABEL_DISPLACEMENT_VERTICAL), (y_max - Y_MIN) / 2 + Y_MIN),
                "\\Large" + label, color,
                TextAlignment::CENTER, rot);
            scene += *text;
        }
    }

    return scene;
}

#endif // AXIS_CPP
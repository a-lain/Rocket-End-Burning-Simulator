#ifndef GRAPHIC_CPP
#define GRAPHIC_CPP

#include "Graphic.hpp"

using namespace CPGF::Plot2d;
using namespace CPGF::AffineSpace;
using namespace CPGF::Objects2d;
using namespace CPGF;

Graphic::Graphic(const bool show_legend):
    show_legend(show_legend), legend_position(LegendPosition::ABOVE)
{

}

void Graphic::add(GraphicObject* object, Axis* Y, Axis* X)
{
    // We add the graphic object.
    graphic_objects.push_back(std::tuple<GraphicObject*, Axis*, Axis*>(object, Y, X));
    
    // We search the list of axes and look for Y and X.
    bool Y_included = false;
    bool X_included = false;
    for (unsigned int i = 0; i < axes.size(); i++)
    {
        if (axes[i] == Y)
        {
            Y_included = true;
        }
        else if (axes[i] == X)
        {
            X_included = true;
        }
    }
    if (!Y_included)
    {
        axes.push_back(Y);
    }
    if (!X_included)
    {
        axes.push_back(X);
    }
}

Scene2d Graphic::render_to_scene() const
{
    Scene2d scene;

    // We update the max and min value of each axis.
    GraphicObject* object;
    Axis* X_axis;
    Axis* Y_axis;
    for (unsigned int i = 0; i < graphic_objects.size(); i++)
    {
        std::tie(object, Y_axis, X_axis) = graphic_objects[i];
        Y_axis->update_max_min(object->y_max(), object->y_min());
        X_axis->update_max_min(object->x_max(), object->x_min());
    }

    // We render the axes.
    for (unsigned int i = 0; i < axes.size(); i++)
    {
        axes[i]->calculate_transformations();
        scene += axes[i]->render_to_scene();
    }

    // We render every object.
    for (unsigned int i = 0; i < graphic_objects.size(); i++)
    {
        std::tie(object, Y_axis, X_axis) = graphic_objects[i];
        std::function<Point2d(const Point2d&)> transform =
            [&] (const Point2d& P){return Point2d(X_axis->axis_transform()(P.x), Y_axis->axis_transform()(P.y));};
        scene += object->render_to_scene(transform, X_axis->get_min_value(), X_axis->get_max_value(),
            Y_axis->get_min_value(), Y_axis->get_max_value());
    }

    // We render the legend.
    if (show_legend)
    {
        unsigned int N_legends = 0;
        // We count how many legends we have to show.
        for (unsigned int i = 0; i < graphic_objects.size(); i++)
        {
            std::tie(object, Y_axis, X_axis) = graphic_objects[i];
            if (object->legend != "")
            {
                N_legends++;
            }
        }

        if (N_legends != 0)
        {
            double pos_max = 0;
            double pos_min = 1;
            double aspect_ratio_max = 0;
            double scale_max = 0;
            double x_miniatures;
            double y_miniatures;
            double x_med;
            unsigned int j;
            unsigned int k;
            Text* text;
            Object2d* miniature;
            switch (legend_position)
            {
                case LegendPosition::RIGHT:
                    for (unsigned int i = 0; i < axes.size(); i++)
                    {
                        if (axes[i]->axis_type == AxisType::VERTICAL
                            && axes[i]->number_position == NumberPosition::RIGHT
                            && axes[i]->position > pos_max)
                        {
                            pos_max = axes[i]->position;
                        }
                        if (axes[i]->aspect_ratio > aspect_ratio_max)
                        {
                            aspect_ratio_max = axes[i]->aspect_ratio;
                        }
                        if (axes[i]->scale > scale_max)
                        {
                            scale_max = axes[i]->scale;
                        }
                    }
                    x_miniatures = std::max(pos_max*(scale_max*Axis::X_MAX - Axis::X_MIN) + Axis::X_MIN
                        + Axis::LABEL_DISPLACEMENT_VERTICAL, scale_max*Axis::X_MAX) + LEGEND_MARGIN;
                    y_miniatures = aspect_ratio_max*(scale_max*Axis::X_MAX - Axis::X_MIN)/2 + Axis::Y_MIN
                        + N_legends*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE;
                    for (unsigned int i = 0; i < graphic_objects.size(); i++)
                    {
                        std::tie(object, Y_axis, X_axis) = graphic_objects[i];
                        if (object->legend != "")
                        {
                            miniature = new Object2d(object->miniature(Point2d(x_miniatures, y_miniatures)));
                            text = new Text(Point2d(x_miniatures + GraphicObject::MINIATURE_HALF_WIDTH, y_miniatures),
                                object->legend, Color::BLACK, TextAlignment::LEFT, 0);
                            scene += *miniature;
                            scene += *text;
                            y_miniatures -= 2*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE;
                        }
                    }
                    break;
                case LegendPosition::LEFT:
                    for (unsigned int i = 0; i < axes.size(); i++)
                    {
                        if (axes[i]->axis_type == AxisType::VERTICAL 
                            && axes[i]->number_position == NumberPosition::LEFT
                            && axes[i]->position < pos_min)
                        {
                            pos_min = axes[i]->position;
                        }
                        if (axes[i]->aspect_ratio > aspect_ratio_max)
                        {
                            aspect_ratio_max = axes[i]->aspect_ratio;
                        }
                        if (axes[i]->scale > scale_max)
                        {
                            scale_max = axes[i]->scale;
                        }
                    }
                    x_miniatures = std::min(pos_min*(scale_max*Axis::X_MAX - Axis::X_MIN) + Axis::X_MIN
                        - Axis::LABEL_DISPLACEMENT_VERTICAL, Axis::X_MIN) - LEGEND_MARGIN;
                    y_miniatures = aspect_ratio_max*(scale_max*Axis::X_MAX - Axis::X_MIN)/2 + Axis::Y_MIN
                        + N_legends*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE;
                    for (unsigned int i = 0; i < graphic_objects.size(); i++)
                    {
                        std::tie(object, Y_axis, X_axis) = graphic_objects[i];
                        if (object->legend != "")
                        {
                            miniature = new Object2d(object->miniature(Point2d(x_miniatures, y_miniatures)));
                            text = new Text(Point2d(x_miniatures - GraphicObject::MINIATURE_HALF_WIDTH, y_miniatures),
                                object->legend, Color::BLACK, TextAlignment::RIGHT, 0);
                            scene += *miniature;
                            scene += *text;
                            y_miniatures -= 2*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE;
                        }
                    }
                    break;
                case LegendPosition::ABOVE:
                    for (unsigned int i = 0; i < axes.size(); i++)
                    {
                        if (axes[i]->axis_type == AxisType::HORIZONTAL
                            && axes[i]->number_position == NumberPosition::LEFT
                            && axes[i]->position > pos_max)
                        {
                            pos_max = axes[i]->position;
                        }
                        if (axes[i]->aspect_ratio > aspect_ratio_max)
                        {
                            aspect_ratio_max = axes[i]->aspect_ratio;
                        }
                        if (axes[i]->scale > scale_max)
                        {
                            scale_max = axes[i]->scale;
                        }
                    }
                    y_miniatures = std::max(pos_max*aspect_ratio_max*(scale_max*Axis::X_MAX - Axis::X_MIN) + Axis::Y_MIN
                        + Axis::LABEL_DISPLACEMENT_HORIZONTAL, aspect_ratio_max*(scale_max*Axis::X_MAX - Axis::X_MIN)+ Axis::Y_MIN)
                        + LEGEND_MARGIN + 2*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE * (N_legends / 2);
                    x_med = (scale_max*Axis::X_MAX + Axis::X_MIN) / 2;
                    
                    j = 0;
                    if (N_legends % 2 != 0)
                    {
                        while (std::get<0>(graphic_objects[j])->legend == "")
                        {
                            j++;
                        }
                        std::tie(object, Y_axis, X_axis) = graphic_objects[j];
                        miniature = new Object2d(object->miniature(
                            Point2d(x_med,
                                y_miniatures + 2*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE)));
                        text = new Text(Point2d(x_med + GraphicObject::MINIATURE_HALF_WIDTH, y_miniatures + 2*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE),
                            object->legend, Color::BLACK, TextAlignment::LEFT);
                        scene += *miniature;
                        scene += *text;
                        j++;                     
                    }
                    while (j < graphic_objects.size() && std::get<0>(graphic_objects[j])->legend == "")
                    {
                        j++;
                    }
                    while (j < graphic_objects.size())
                    {
                        std::tie(object, Y_axis, X_axis) = graphic_objects[j];
                        miniature = new Object2d(object->miniature(
                            Point2d(0, y_miniatures)));
                        text = new Text(Point2d(GraphicObject::MINIATURE_HALF_WIDTH, y_miniatures), object->legend, Color::BLACK,
                            TextAlignment::LEFT);
                        scene += *miniature;
                        scene += *text;
                        j++;
                        while (std::get<0>(graphic_objects[j])->legend == "")
                        {
                            j++;
                        } 
                        std::tie(object, Y_axis, X_axis) = graphic_objects[j];
                        miniature = new Object2d(object->miniature(
                            Point2d(x_med, y_miniatures)));
                        text = new Text(Point2d(x_med + GraphicObject::MINIATURE_HALF_WIDTH, y_miniatures), object->legend, Color::BLACK,
                            TextAlignment::LEFT);
                        scene += *miniature;
                        scene += *text;   
                        y_miniatures -= 2*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE;
                        j++;
                        while (j < graphic_objects.size() && std::get<0>(graphic_objects[j])->legend == "")
                        {
                            j++;
                        }
                    }
                    break;
                case LegendPosition::BELOW:
                    for (unsigned int i = 0; i < axes.size(); i++)
                    {
                        if (axes[i]->axis_type == AxisType::HORIZONTAL
                            && axes[i]->number_position == NumberPosition::RIGHT
                            && axes[i]->position < pos_min)
                        {
                            pos_min = axes[i]->position;
                        }
                        if (axes[i]->aspect_ratio > aspect_ratio_max)
                        {
                            aspect_ratio_max = axes[i]->aspect_ratio;
                        }
                        if (axes[i]->scale > scale_max)
                        {
                            scale_max = axes[i]->scale;
                        }
                    }
                    y_miniatures = std::min(pos_min*aspect_ratio_max*(scale_max*Axis::X_MAX - Axis::X_MIN) + Axis::Y_MIN
                        - Axis::LABEL_DISPLACEMENT_HORIZONTAL, Axis::X_MIN)
                        - LEGEND_MARGIN;
                    x_med = (scale_max*Axis::X_MAX + Axis::X_MIN) / 2;

                    j = 0;
                    k = 0;
                    while (j < graphic_objects.size() && k < 2*(N_legends/2) && std::get<0>(graphic_objects[j])->legend == "")
                    {
                        j++;
                    }
                    while (j < graphic_objects.size() && k < 2*(N_legends/2))
                    {
                        std::tie(object, Y_axis, X_axis) = graphic_objects[j];
                        miniature = new Object2d(object->miniature(
                            Point2d(0, y_miniatures)));
                        text = new Text(Point2d(GraphicObject::MINIATURE_HALF_WIDTH, y_miniatures), object->legend, Color::BLACK,
                            TextAlignment::LEFT);
                        scene += *miniature;
                        scene += *text;  
                        j++; 
                        while (std::get<0>(graphic_objects[j])->legend == "")
                        {
                            j++;
                        } 
                        std::tie(object, Y_axis, X_axis) = graphic_objects[j];
                        miniature = new Object2d(object->miniature(
                            Point2d(x_med, y_miniatures)));
                        text = new Text(Point2d(x_med + GraphicObject::MINIATURE_HALF_WIDTH, y_miniatures), object->legend, Color::BLACK,
                            TextAlignment::LEFT);
                        scene += *miniature;
                        scene += *text;   
                        y_miniatures -= 2*LEGEND_VERTICAL_DISPLACEMENT_PER_LINE;
                        j++;
                        k += 2;
                        while (j < graphic_objects.size() && k < 2*(N_legends/2) && std::get<0>(graphic_objects[j])->legend == "")
                        {
                            j++;
                        }
                    }
                    if (N_legends %= 2)
                    {
                        while (std::get<0>(graphic_objects[j])->legend == "")
                        {
                            j++;
                        }
                        std::tie(object, Y_axis, X_axis) = graphic_objects[j];
                        miniature = new Object2d(object->miniature(
                            Point2d(x_med,
                                y_miniatures)));
                        text = new Text(Point2d(x_med + GraphicObject::MINIATURE_HALF_WIDTH, y_miniatures),
                            object->legend, Color::BLACK,
                            TextAlignment::LEFT);
                        scene += *miniature;
                        scene += *text;
                    }
                    break;
            }        
        }
    }

    return scene;
}

#endif // GRAPHIC_CPP
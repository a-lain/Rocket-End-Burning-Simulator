#ifndef GRAPHIC_HPP
#define GRAPHIC_HPP

#include <vector>
#include "GraphicObject.hpp"
#include "../Scene2d.hpp"
#include "Axis.hpp"

namespace CPGF
{
    namespace Plot2d
    {
        enum class LegendPosition
        {
            LEFT,
            ABOVE,
            RIGHT,
            BELOW,
        };

        class Graphic
        {
            public:
            std::vector<std::tuple<GraphicObject*, Axis*, Axis*>> graphic_objects;

            bool show_legend;

            LegendPosition legend_position;

            void add(GraphicObject* object, Axis* Y, Axis* X);

            Scene2d render_to_scene() const;

            Graphic(const bool show_legend = false);

            std::vector<Axis*> axes;

            constexpr static const double LEGEND_MARGIN = 0.25;
            constexpr static const double LEGEND_VERTICAL_DISPLACEMENT_PER_LINE = 0.25;
            constexpr static const double CHARACTER_WIDTH = 0.1;
        };
    }
}

#endif // GRAPHIC_HPP
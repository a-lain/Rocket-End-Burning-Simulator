#ifndef PGFCONF_HPP
#define PGFCONF_HPP

#include <vector>
#include <string>

namespace CPGF
{
    /**
     * @brief Determines whether a path has to be drawn (only the contour
     * is drawn and the interior remains white) or to be filled (the
     * contour is not drawn but the interior is).
     * 
     */
    enum class DrawType
    {
        DRAW,
        FILL,
    };

    /**
     * @brief This class provides a handful of useful predefined constants
     * to express line width.
     * 
     */
    class LineWidth
    {
        public:
        static const double ULTRA_THIN;
        static const double VERY_THIN;
        static const double THIN;
        static const double SEMITHICK;
        static const double THICK;
        static const double VERY_THICK;
        static const double ULTRA_THICK;
    };

    /**
     * @brief This determines how the extremes of an open path are drawn.
     * 
     */
    enum class LineCap
    {
        ROUND,
        RECT,
        BUTT
    };

    /**
     * @brief This determines how the junction between two non-parallel
     * lines is drawn.
     * 
     */
    enum class LineJoin
    {
        ROUND,
        BEVEL,
        MITER
    };

    // Remember miter limit. Currently not available.

    /**
     * @brief An object used to represent an rgb color. r, g and b
     * must all be real numbers between 0 and 1.
     * 
     */
    class Color
    {
        public:
        /**
         * @brief The amount of red of the color. \f$0\le r\le1\f$.
         * 
         */
        double r;

        /**
         * @brief The amount of green of the color. \f$0\le g\le1\f$.
         * 
         */
        double g;

        /**
         * @brief The amount of blue of the color. \f$0\le b\le1\f$.
         * 
         */
        double b;

        /**
         * @brief Returns a string representation of the color.
         * 
         * @return std::string 
         */
        std::string to_string();

        /**
         * @brief Returns the red color (1,0,0).
         * 
         */
        Color();

        /**
         * @brief Returns the (r,g,b) color.
         * 
         * @param r a real number between zero and one.
         * @param g a real number between zero and one.
         * @param b a real number between zero and one.
         */
        Color(const double r, const double g, const double b);

        /**
         * @brief Returns a mix of colors A and B in proportions
         * alpha and 1 - alpha.
         * 
         * @param A 
         * @param B 
         * @param alpha 
         * @return Color 
         */
        static Color mix(const Color& A, const Color& B, const double alpha);

        static Color from_RGB(const unsigned char R, const unsigned char G, const unsigned char B);

        // Color list taken from https://www.rapidtables.com/web/color/RGB_Color.html
        static const Color MAROON;
        static const Color DARK_RED;
        static const Color BROWN;
        static const Color FIREBRICK;
        static const Color CRIMSON;
        static const Color RED;
        static const Color TOMATO;
        static const Color CORAL;
        static const Color INDIAN_RED;
        static const Color LIGHT_CORAL;
        static const Color DARK_SALMON;
        static const Color SALMON;
        static const Color LIGHT_SALMON;
        static const Color ORANGE_RED;
        static const Color DARK_ORANGE;
        static const Color ORANGE;
        static const Color GOLD;
        static const Color DARK_GOLDEN_ROD;
        static const Color GOLDEN_ROD;
        static const Color PALE_GOLDEN_ROD;
        static const Color DARK_KHAKI;
        static const Color KHAKI;
        static const Color OLIVE;
        static const Color YELLOW;
        static const Color YELLOW_GREEN;
        static const Color DARK_OLIVE_GREEN;
        static const Color OLIVE_DRAB;
        static const Color LAWN_GREEN;
        static const Color CHART_REUSE;
        static const Color GREEN_YELLOW;
        static const Color DARK_GREEN;
        static const Color GREEN;
        static const Color FOREST_GREEN;
        static const Color LIME;
        static const Color LIME_GREEN;
        static const Color LIGHT_GREEN;
        static const Color PALE_GREEN;
        static const Color DARK_SEA_GREEN;
        static const Color MEDIUM_SPRING_GREEN;
        static const Color SPRING_GREEN;
        static const Color SEA_GREEN;
        static const Color MEDIUM_AQUA_MARINE;
        static const Color MEDIUM_SEA_GREEN;
        static const Color LIGHT_SEA_GREEN;
        static const Color DARK_SLATE_GRAY;
        static const Color TEAL;
        static const Color DARK_CYAN;
        static const Color CYAN;
        static const Color LIGHT_CYAN;
        static const Color DARK_TURQUOISE;
        static const Color TURQUOISE;
        static const Color MEDIUM_TURQUOISE;
        static const Color PALE_TORQUOISE;
        static const Color AQUA_MARINE;
        static const Color POWDER_BLUE;
        static const Color CADET_BLUE;
        static const Color STEEL_BLUE;
        static const Color CORN_FLOWER_BLUE;
        static const Color DEEP_SKY_BLUE;
        static const Color DODGER_BLUE;
        static const Color LIGHT_BLUE;
        static const Color SKY_BLUE;
        static const Color LIGHT_SKY_BLUE;
        static const Color MIDNIGHT_BLUE;
        static const Color NAVY;
        static const Color DARK_BLUE;
        static const Color MEDIUM_BLUE;
        static const Color BLUE;
        static const Color ROYAL_BLUE;
        static const Color BLUE_VIOLET;
        static const Color INDIGO;
        static const Color DARK_SLATE_BLUE;
        static const Color SLATE_BLUE;
        static const Color MEDIUM_SLATE_BLUE;
        static const Color MEDIUM_PURPLE;
        static const Color DARK_MAGENTA;
        static const Color DARK_VIOLET;
        static const Color DARK_ORCHID;
        static const Color MEDIUM_ORCHID;
        static const Color PURPLE;
        static const Color THISTLE;
        static const Color PLUM;
        static const Color VIOLET;
        static const Color MAGENTA;
        static const Color ORCHID;
        static const Color MEDIUM_VIOLET_RED;
        static const Color PALE_VIOLET_RED;
        static const Color DEEP_PINK;
        static const Color HOT_PINK;
        static const Color LIGHT_PINK;
        static const Color PINK;
        static const Color ANTIQUE_WHITE;
        static const Color BEIGE;
        static const Color BISQUE;
        static const Color BLANCHED_ALMOND;
        static const Color WHEAT;
        static const Color CORN_SILK;
        static const Color LEMON_CHIFFON;
        static const Color LIGHT_GOLDEN_ROD_YELLOW;
        static const Color LIGHT_YELLOW;
        static const Color SADDLE_BROWN;
        static const Color SIENNA;
        static const Color CHOCOLATE;
        static const Color PERU;
        static const Color SANDY_BROWN;
        static const Color BURLY_WOOD;
        static const Color TAN;
        static const Color ROSY_BROWN;
        static const Color MOCCASIN;
        static const Color NAVAJO_WHITE;
        static const Color PEACH_STUFF;
        static const Color MISTY_ROSE;
        static const Color LAVENDER_BLUSH;
        static const Color LINEN;
        static const Color OLD_LACE;
        static const Color PAPAYA_WHIP;
        static const Color SEA_SHELL;
        static const Color MINT_CREAM;
        static const Color SLATE_GRAY;
        static const Color LIGHT_SLATE_GRAY;
        static const Color LIGHT_STEEL_BLUE;
        static const Color LAVENDER;
        static const Color FLORAL_WHITE;
        static const Color ALICE_BLUE;
        static const Color GHOST_WHITE;
        static const Color HONEYDEW;
        static const Color IVORY;
        static const Color AZURE;
        static const Color SNOW;
        static const Color BLACK;
        static const Color DIM_GRAY;
        static const Color GRAY;
        static const Color DARK_GRAY;
        static const Color SILVER;
        static const Color LIGHT_GRAY;
        static const Color GAINSBORO;
        static const Color WHITE_SMOKE;
        static const Color WHITE;
    };

    /**
     * @brief This class provides a handful of useful predefined constants
     * to express DashPatterns.
     * 
     */
    class DashPatterns
    {
        public:
        static const std::vector<double> SOLID;
        static const std::vector<double> DASHED;
        static const std::vector<double> DOTTED;
    };

    /**
     * @brief This object is used to store all information needed to draw
     * a path.
     * 
     */
    class PGFConf
    {
        public:
        DrawType draw_type;
        LineCap line_cap;
        LineJoin line_join;
        Color color;

        /**
         * @brief This double array is used to decide whether the line
         * drawn is continuous or discontinuous (i.e. it has gaps).
         * 
         * \todo Explain this better.
         * 
         */
        std::vector<double> dash_pattern;

        /**
         * @brief 
         * 
         * \todo Explain this.
         * 
         */
        double dash_phase;

        /**
         * @brief This determines how transparent the path drawn is. It must
         * be a real number between zero and one. One means it is fully opaque
         * and zero represents full transparency.
         * 
         */
        double opacity;
        
        /**
         * @brief This determines the thickness of the line drawn. For default
         * predefined values, you can use the class LineWidth.
         * 
         */
        double line_width;

        PGFConf(const DrawType draw_type = DrawType::DRAW,
            const Color& color = Color::BLACK,
            const double opacity = 1,
            const double line_width = LineWidth::SEMITHICK,
            const std::vector<double>& dash_pattern = std::vector<double>(),
            const double dash_phase = 0,
            const LineCap line_cap = LineCap::BUTT,
            const LineJoin line_join = LineJoin::BEVEL
            );
    };
}


#endif // TIKZCONF_HPP

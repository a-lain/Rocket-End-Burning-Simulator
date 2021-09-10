#ifndef TEXT_HPP
#define TEXT_HPP

#include "../AffineSpace2d/Point2d.hpp"
#include "../PGFBasics/PGFConf.hpp"

namespace CPGF
{
    /**
     * @brief 
     * 
     * Part of the following comments are copied from the TikZ/PGF manual.
     * 
     */
    enum class TextAlignment
    {
        /**
         * @brief The default TextAlignment
         * 
         */
        CENTER,

        /**
         * @brief The key causes the text box to be placed such that its left border is on the origin.
         * 
         */
        LEFT,

        /**
         * @brief The key causes the text box to be placed such that its right border is on the origin.
         * 
         */
        RIGHT,

        /**
         * @brief This key causes the text box to be placed such that its top is on the origin.
         * 
         */
        TOP,

        /**
         * @brief This key causes the text box to be placed such that its bottom is on the origin.
         * 
         */
        BOTTOM,

        /**
         * @brief This key causes the text box to be placed such that its baseline is on the origin.
         * 
         */
        BASE,

        /**
         * @brief This key causes the text box to be placed such that its top left corner is on the origin.
         * 
         */
        TOP_LEFT,

        /**
         * @brief This key causes the text box to be placed such that its top right corner is on the origin.
         * 
         */
        TOP_RIGHT,

        /**
         * @brief This key causes the text box to be placed such that its bottom left corner is on the origin.
         * 
         */
        BOTTOM_LEFT,

        /**
         * @brief This key causes the text box to be placed such that its bottom right corner is on the origin.
         * 
         */
        BOTTOM_RIGHT,

        /**
         * @brief This key causes the text box to be placed such that its base left corner is on the origin.
         * 
         */
        BASE_LEFT,

        /**
         * @brief This key causes the text box to be placed such that its base right corner is on the origin.
         * 
         */
        BASE_RIGHT,
    };

    /**
     * @brief
     * 
     * @warning Text cannot be scaled or rotated, only moved.
     * 
     */
    class Text
    {
        public:
        /**
         * @brief The text that is going to be printed.
         * 
         */
        std::string text;

        /**
         * @brief The text position.
         * 
         */
        AffineSpace::Point2d pos;

        /**
         * @brief The text rotation in degrees.
         * 
         */
        double rot;

        /**
         * @brief The text alignment form.
         * 
         */
        TextAlignment text_alignment;

        /**
         * @brief The text color.
         * 
         */
        Color color;

        Text(const AffineSpace::Point2d& pos = AffineSpace::Point2d(0,0),
            const std::string& text = "",
            const Color& color = Color::BLACK,
            const TextAlignment text_alignment = TextAlignment::CENTER,
            const double rot = 0);

        std::string render_to_string() const;
    };
}

#endif // TEXT_HPP
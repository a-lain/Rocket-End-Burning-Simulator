#ifndef TEXT_CPP
#define TEXT_CPP

#include "Text.hpp"

using namespace CPGF;
using namespace CPGF::AffineSpace;

Text::Text(const Point2d& pos, const std::string& text,
            const Color& color,
            const TextAlignment text_alignment, const double rot):
    text(text), pos(pos), rot(rot), text_alignment(text_alignment),
    color(color)
{

}

std::string Text::render_to_string() const
{
    // We start with previous commands for the path
    std::string text = "\\begin{pgfscope}\n";

    // Color.
    text += "\\definecolor{tempcolor}{rgb}{" + std::to_string(color.r) + ", " +
        std::to_string(color.g) + ", " + std::to_string(color.b) + "}\n" +
        "\\pgfsetcolor{tempcolor}\n";

    // Text position.
    text += "\\pgftext[at={\\pgfpoint{" + std::to_string(pos.x) + "cm}{"
        + std::to_string(pos.y) + "cm}}";

    // Text alignment.
    switch (text_alignment)
    {
        case TextAlignment::CENTER:
            text += ", ";
            break;
        case TextAlignment::LEFT:
            text += ", left, ";
            break;
        case TextAlignment::RIGHT:
            text += ", right, ";
            break;
        case TextAlignment::TOP:
            text += ", top, ";
            break;
        case TextAlignment::BOTTOM:
            text += ", bottom, ";
            break;
        case TextAlignment::BASE:
            text += ", base, ";
            break;
        case TextAlignment::TOP_LEFT:
            text += ", top, left, ";
            break;
        case TextAlignment::TOP_RIGHT:
            text += ", top, right, ";
            break;
        case TextAlignment::BOTTOM_LEFT:
            text += ", bottom, left, ";
            break;
        case TextAlignment::BOTTOM_RIGHT:
            text += ", bottom, right, ";
            break;
        case TextAlignment::BASE_LEFT:
            text += ", base, left, ";
            break;
        case TextAlignment::BASE_RIGHT:
            text += ", base, right, ";
            break;
    }

    // Rotation.
    text += "rotate=" + std::to_string(rot) + "]";

    // Actual text.
    text += "{" + this->text + "}\n";

    text += "\\end{pgfscope}\n";
    return text;
}

#endif // TEXT_CPP
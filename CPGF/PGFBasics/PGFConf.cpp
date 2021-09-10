#ifndef PGFCONF_CPP
#define PGFCONF_CPP

#include "PGFConf.hpp"

using namespace CPGF;

namespace CPGF
{
    const double LineWidth::ULTRA_THIN = 0.025;
    const double LineWidth::VERY_THIN = 0.05;
    const double LineWidth::THIN = 0.1;
    const double LineWidth::SEMITHICK = 0.2;
    const double LineWidth::THICK = 0.4;
    const double LineWidth::VERY_THICK = 0.8;
    const double LineWidth::ULTRA_THICK = 1.2;

    const Color Color::MAROON = Color::from_RGB(128, 0, 0);
    const Color Color::DARK_RED = Color::from_RGB(139, 0, 0);
    const Color Color::BROWN = Color::from_RGB(165, 42, 42);
    const Color Color::FIREBRICK = Color::from_RGB(178, 34, 34);
    const Color Color::CRIMSON = Color::from_RGB(220, 20, 60);
    const Color Color::RED = Color::from_RGB(255, 0 , 0);
    const Color Color::TOMATO = Color::from_RGB(255, 99, 71);
    const Color Color::CORAL = Color::from_RGB(255, 127, 80);
    const Color Color::INDIAN_RED = Color::from_RGB(205, 92, 92);
    const Color Color::LIGHT_CORAL = Color::from_RGB(240, 128, 128);
    const Color Color::DARK_SALMON = Color::from_RGB(233, 150, 122);
    const Color Color::SALMON = Color::from_RGB(250, 128, 114);
    const Color Color::LIGHT_SALMON = Color::from_RGB(255, 160, 122);
    const Color Color::ORANGE_RED = Color::from_RGB(255, 69, 0);
    const Color Color::DARK_ORANGE = Color::from_RGB(255, 140, 0);
    const Color Color::ORANGE = Color::from_RGB(255, 165, 0);
    const Color Color::GOLD = Color::from_RGB(255, 215, 0);
    const Color Color::DARK_GOLDEN_ROD = Color::from_RGB(184, 134, 11);
    const Color Color::GOLDEN_ROD = Color::from_RGB(218, 165, 32);
    const Color Color::PALE_GOLDEN_ROD = Color::from_RGB(238, 232, 170);
    const Color Color::DARK_KHAKI = Color::from_RGB(189, 183, 107);
    const Color Color::KHAKI = Color::from_RGB(240, 230, 140);
    const Color Color::OLIVE = Color::from_RGB(128, 128, 0);
    const Color Color::YELLOW = Color::from_RGB(255, 255, 0);
    const Color Color::YELLOW_GREEN = Color::from_RGB(154, 205, 50);
    const Color Color::DARK_OLIVE_GREEN = Color::from_RGB(85, 107, 47);
    const Color Color::OLIVE_DRAB = Color::from_RGB(107, 142, 35);
    const Color Color::LAWN_GREEN = Color::from_RGB(124, 252, 0);
    const Color Color::CHART_REUSE = Color::from_RGB(127, 255, 0);
    const Color Color::GREEN_YELLOW = Color::from_RGB(173, 255, 47);
    const Color Color::DARK_GREEN = Color::from_RGB(0, 100, 0);
    const Color Color::GREEN = Color::from_RGB(0, 128, 0);
    const Color Color::FOREST_GREEN = Color::from_RGB(34, 139, 34);
    const Color Color::LIME = Color::from_RGB(0, 255, 0);
    const Color Color::LIME_GREEN = Color::from_RGB(50, 205, 50);
    const Color Color::LIGHT_GREEN = Color::from_RGB(144, 238, 144);
    const Color Color::PALE_GREEN = Color::from_RGB(152, 251, 152);
    const Color Color::DARK_SEA_GREEN = Color::from_RGB(143, 188, 143);
    const Color Color::MEDIUM_SPRING_GREEN = Color::from_RGB(0, 250, 154);
    const Color Color::SPRING_GREEN = Color::from_RGB(0, 255, 127);
    const Color Color::SEA_GREEN = Color::from_RGB(46, 139, 87);
    const Color Color::MEDIUM_AQUA_MARINE = Color::from_RGB(102, 205, 170);
    const Color Color::MEDIUM_SEA_GREEN = Color::from_RGB(60, 179, 113);
    const Color Color::LIGHT_SEA_GREEN = Color::from_RGB(32, 178, 170);
    const Color Color::DARK_SLATE_GRAY = Color::from_RGB(47, 79, 79);
    const Color Color::TEAL = Color::from_RGB(0, 128, 128);
    const Color Color::DARK_CYAN = Color::from_RGB(0, 139, 139);
    const Color Color::CYAN = Color::from_RGB(0, 255, 255);
    const Color Color::LIGHT_CYAN = Color::from_RGB(244, 255, 255);
    const Color Color::DARK_TURQUOISE = Color::from_RGB(0, 206, 209);
    const Color Color::TURQUOISE = Color::from_RGB(64, 224, 208);
    const Color Color::MEDIUM_TURQUOISE = Color::from_RGB(72, 209, 204);
    const Color Color::PALE_TORQUOISE = Color::from_RGB(175, 238, 238);
    const Color Color::AQUA_MARINE = Color::from_RGB(127, 255, 212);
    const Color Color::POWDER_BLUE = Color::from_RGB(176, 224, 230);
    const Color Color::CADET_BLUE = Color::from_RGB(95, 158, 160);
    const Color Color::STEEL_BLUE = Color::from_RGB(70, 130, 180);
    const Color Color::CORN_FLOWER_BLUE = Color::from_RGB(100, 149, 237);
    const Color Color::DEEP_SKY_BLUE = Color::from_RGB(0, 191, 255);
    const Color Color::DODGER_BLUE = Color::from_RGB(30, 144, 255);
    const Color Color::LIGHT_BLUE = Color::from_RGB(173, 216, 230);
    const Color Color::SKY_BLUE = Color::from_RGB(135, 206, 235);
    const Color Color::LIGHT_SKY_BLUE = Color::from_RGB(135, 206, 250);
    const Color Color::MIDNIGHT_BLUE = Color::from_RGB(25, 25, 112);
    const Color Color::NAVY = Color::from_RGB(0, 0, 128);
    const Color Color::DARK_BLUE = Color::from_RGB(0, 0, 139);
    const Color Color::MEDIUM_BLUE = Color::from_RGB(0, 0, 205);
    const Color Color::BLUE = Color::from_RGB(0, 0, 255);
    const Color Color::ROYAL_BLUE = Color::from_RGB(65, 105, 225);
    const Color Color::BLUE_VIOLET = Color::from_RGB(138, 43, 226);
    const Color Color::INDIGO = Color::from_RGB(75, 0, 130);
    const Color Color::DARK_SLATE_BLUE = Color::from_RGB(72, 61, 139);
    const Color Color::SLATE_BLUE = Color::from_RGB(106, 90, 205);
    const Color Color::MEDIUM_SLATE_BLUE = Color::from_RGB(123, 104, 238);
    const Color Color::MEDIUM_PURPLE = Color::from_RGB(147, 112, 219);
    const Color Color::DARK_MAGENTA = Color::from_RGB(139, 0, 139);
    const Color Color::DARK_VIOLET = Color::from_RGB(148, 0, 211);
    const Color Color::DARK_ORCHID = Color::from_RGB(153, 50, 204);
    const Color Color::MEDIUM_ORCHID = Color::from_RGB(186, 85, 211);
    const Color Color::PURPLE = Color::from_RGB(128, 0, 128);
    const Color Color::THISTLE = Color::from_RGB(216, 191, 216);
    const Color Color::PLUM = Color::from_RGB(221, 160, 221);
    const Color Color::VIOLET = Color::from_RGB(238, 130, 238);
    const Color Color::MAGENTA = Color::from_RGB(255, 0, 255);
    const Color Color::ORCHID = Color::from_RGB(218, 112, 214);
    const Color Color::MEDIUM_VIOLET_RED = Color::from_RGB(199, 21, 133);
    const Color Color::PALE_VIOLET_RED = Color::from_RGB(219, 112, 147);
    const Color Color::DEEP_PINK = Color::from_RGB(255, 20, 147);
    const Color Color::HOT_PINK = Color::from_RGB(255, 105, 180);
    const Color Color::LIGHT_PINK = Color::from_RGB(255, 182, 193);
    const Color Color::PINK = Color::from_RGB(255, 192, 203);
    const Color Color::ANTIQUE_WHITE = Color::from_RGB(250, 235, 215);
    const Color Color::BEIGE = Color::from_RGB(245, 245, 220);
    const Color Color::BISQUE = Color::from_RGB(255, 228, 196);
    const Color Color::BLANCHED_ALMOND = Color::from_RGB(255, 235, 205);
    const Color Color::WHEAT = Color::from_RGB(245, 222, 179);
    const Color Color::CORN_SILK = Color::from_RGB(255, 248, 220);
    const Color Color::LEMON_CHIFFON = Color::from_RGB(255, 250, 205);
    const Color Color::LIGHT_GOLDEN_ROD_YELLOW = Color::from_RGB(250, 250, 210);
    const Color Color::LIGHT_YELLOW = Color::from_RGB(255, 255, 224);
    const Color Color::SADDLE_BROWN = Color::from_RGB(139, 69, 19);
    const Color Color::SIENNA = Color::from_RGB(160, 82, 45);
    const Color Color::CHOCOLATE = Color::from_RGB(210, 105, 30);
    const Color Color::PERU = Color::from_RGB(205, 133, 63);
    const Color Color::SANDY_BROWN = Color::from_RGB(244, 164, 96);
    const Color Color::BURLY_WOOD = Color::from_RGB(222, 184, 135);
    const Color Color::TAN = Color::from_RGB(210, 180, 140);
    const Color Color::ROSY_BROWN = Color::from_RGB(188, 143, 143);
    const Color Color::MOCCASIN = Color::from_RGB(255, 228, 181);
    const Color Color::NAVAJO_WHITE = Color::from_RGB(255, 222, 173);
    const Color Color::PEACH_STUFF = Color::from_RGB(255, 218, 185);
    const Color Color::MISTY_ROSE = Color::from_RGB(255, 228, 225);
    const Color Color::LAVENDER_BLUSH = Color::from_RGB(255, 240, 245);
    const Color Color::LINEN = Color::from_RGB(250, 240, 230);
    const Color Color::OLD_LACE = Color::from_RGB(253, 245, 230);
    const Color Color::PAPAYA_WHIP = Color::from_RGB(255, 239, 213);
    const Color Color::SEA_SHELL = Color::from_RGB(255, 245, 238);
    const Color Color::MINT_CREAM = Color::from_RGB(245, 255, 250);
    const Color Color::SLATE_GRAY = Color::from_RGB(112, 128, 144);
    const Color Color::LIGHT_SLATE_GRAY = Color::from_RGB(119, 136, 153);
    const Color Color::LIGHT_STEEL_BLUE = Color::from_RGB(176, 196, 222);
    const Color Color::LAVENDER = Color::from_RGB(230, 230, 250);
    const Color Color::FLORAL_WHITE = Color::from_RGB(255, 250, 240);
    const Color Color::ALICE_BLUE = Color::from_RGB(240, 248, 255);
    const Color Color::GHOST_WHITE = Color::from_RGB(248, 248, 255);
    const Color Color::HONEYDEW = Color::from_RGB(240, 255, 240);
    const Color Color::IVORY = Color::from_RGB(255, 255, 240);
    const Color Color::AZURE = Color::from_RGB(240, 255, 255);
    const Color Color::SNOW = Color::from_RGB(255, 250, 250);
    const Color Color::BLACK = Color::from_RGB(0, 0, 0);
    const Color Color::DIM_GRAY = Color::from_RGB(105, 105, 105);
    const Color Color::GRAY = Color::from_RGB(128, 128, 128);
    const Color Color::DARK_GRAY = Color::from_RGB(169, 169, 169);
    const Color Color::SILVER = Color::from_RGB(192, 192, 192);
    const Color Color::LIGHT_GRAY = Color::from_RGB(211, 211, 211);
    const Color Color::GAINSBORO = Color::from_RGB(220, 220, 220);
    const Color Color::WHITE_SMOKE = Color::from_RGB(245, 245, 245);
    const Color Color::WHITE = Color::from_RGB(255, 255, 255);

    const std::vector<double> DashPatterns::SOLID = std::vector<double>();
    const std::vector<double> DashPatterns::DASHED = std::vector<double>({0.3, 0.3});
    const std::vector<double> DashPatterns::DOTTED = std::vector<double>({0.1, 0.1});
}


Color::Color()
{
    this->r = 1;
    this->g = 0;
    this->b = 0;
}

Color::Color(const double r, const double g, const double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color Color::mix(const Color& A, const Color& B, const double alpha)
{
    if (alpha > 1)
    {
        throw "Not a valid mixture! α must be in the interval [0,1]";
    }
    return Color(alpha*A.r + (1 - alpha)*B.r, alpha*A.g + (1 - alpha)*B.g, alpha*A.b + (1 - alpha)*B.b);
}

Color Color::from_RGB(const unsigned char R, const unsigned char G, const unsigned char B)
{
    return Color((double)R / 255, (double)G / 255, (double)B / 255);
}

std::string Color::to_string()
{
    return "C(" + std::to_string(this->r) + ", " + std::to_string(this->g) + ", " + std::to_string(this->b) + ")";
}

PGFConf::PGFConf(const DrawType draw_type,
    const Color& color,
    const double opacity,
    const double line_width,
    const std::vector<double>& dash_pattern,
    const double dash_phase,
    const LineCap line_cap,
    const LineJoin line_join): draw_type(draw_type),
        line_cap(line_cap), line_join(line_join),
        color(color), dash_pattern(dash_pattern),
        dash_phase(dash_phase), opacity(opacity),
        line_width(line_width)
{

}


#endif // TIKZCONF_CPP
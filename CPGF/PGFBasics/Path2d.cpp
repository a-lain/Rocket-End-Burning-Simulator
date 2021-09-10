#ifndef PATH2D_CPP
#define PATH2D_CPP

#include "Path2d.hpp"

using namespace CPGF::AffineSpace;
using namespace CPGF::Basics;


Path2d::Path2d()
{
    
}

Path2d::Path2d(SimpleStroke2d& stroke, const PGFConf& conf):
    strokes(std::vector<SimpleStroke2d*>({&stroke})), conf(conf)
{
    
}

Path2d::Path2d(const std::vector<SimpleStroke2d*>& strokes, const PGFConf& conf):
    strokes(strokes), conf(conf)
{

}

Path2d::Path2d(const Path2d& path)
{
    this->conf = path.conf;
    for (unsigned int i = 0; i < path.strokes.size(); i++)
    {
        this->strokes.push_back(path.strokes[i]->clone());
    }
}

Path2d& Path2d::operator=(const Path2d& path)
{
    *this = Path2d(path);
    return *this;
}

Path2d::~Path2d()
{
    // Possible memory leak here!!!

    // for (unsigned int i = 0; i < strokes.size(); i++)
    // {
    //     // delete strokes[i];
    // }
}

Path2d& Path2d::add_stroke(SimpleStroke2d& stroke)
{
    strokes.push_back(&stroke);
    return *this;
}

Path2d& Path2d::operator+=(SimpleStroke2d& stroke)
{
    return this->add_stroke(stroke);
}

unsigned int Path2d::size() const
{
    return this->strokes.size();
}

Path2d& Path2d::translate(const Vector2d& v)
{
    for (unsigned int i = 0; i < strokes.size(); i++)
    {
        strokes[i]->translate(v);
    }
    return *this;
}

Path2d& Path2d::rotate_with_respect_to(const Point2d& Q, const double theta)
{
    for (unsigned int i = 0; i < strokes.size(); i++)
    {
        strokes[i]->rotate_with_respect_to(Q, theta);
    }
    return *this;
}

Path2d& Path2d::scale_with_respect_to(const Point2d& Q, const Vector2d& s)
{
    for (unsigned int i = 0; i < strokes.size(); i++)
    {
        strokes[i]->scale_with_respect_to(Q, s);
    }
    return *this;
}

Point2d& Path2d::start()
{
    return strokes[0]->start();
}

Point2d Path2d::start() const
{
    return strokes[0]->start();
}

Point2d& Path2d::end()
{
    return strokes[size() - 1]->end();
}

Point2d Path2d::end() const
{
    return strokes[size() - 1]->end();
}

double Path2d::length() const
{
    double res = 0;
    for (unsigned int i = 0; i < size(); i++)
    {
        res += strokes[i]->length();
    }
    return res;
}

double Path2d::area() const
{
    double res = 0;
    for (unsigned int i = 0; i < size(); i++)
    {
        res += strokes[i]->area();
    }
    return res;
}

std::string Path2d::render_to_string() const
{
    // We start with previous commands for the path
    std::string text = "\\begin{pgfscope}\n";

    switch (conf.line_cap)
    {
    case LineCap::BUTT:
        text += "\\pgfsetbuttcap\n";
        break;

    case LineCap::ROUND:
        text += "\\pgfsetroundcap\n";
        break;

    case LineCap::RECT:
        text += "\\pgfsetrectcap\n";
        break;
    }

    switch (conf.line_join)
    {
    case LineJoin::MITER:
        text += "\\pgfsetmiterjoin\n";
        break;

    case LineJoin::ROUND:
        text += "\\pgfsetroundjoin\n";
        break;

    case LineJoin::BEVEL:
        text += "\\pgfsetbeveljoin\n";
        break;
    }

    text += "\\definecolor{tempcolor}{rgb}{" + std::to_string(conf.color.r) + ", " +
        std::to_string(conf.color.g) + ", " + std::to_string(conf.color.b) + "}\n" +
        "\\pgfsetcolor{tempcolor}\n";

    if (conf.dash_pattern.size() != 0)
    {
        text += "\\pgfsetdash{";
        for (unsigned int i = 0; i < conf.dash_pattern.size(); i++)
        {
            text += "{" + std::to_string(conf.dash_pattern[i]) + "cm}";
        }
        text += "}{" + std::to_string(conf.dash_phase) + "cm}\n";

    }

    // We add opacity configuration.
    switch (conf.draw_type)
        {
            case DrawType::DRAW:
                text += "\\pgfsetstrokeopacity{" + std::to_string(conf.opacity) + "}\n";
                break;
            
            case DrawType::FILL:
                text += "\\pgfsetfillopacity{" + std::to_string(conf.opacity) + "}\n";
                break;
        }

    text += "\\pgfsetlinewidth{" + std::to_string(conf.line_width) + "cm}\n";


    // We know begin tracing the actual path
    if (strokes.size() != 0)
    {
        for (unsigned int i = 0; i < strokes.size(); i++)
        {
            text += strokes[i]->render_to_string();
        }

        // We either draw or fill the path
        switch (conf.draw_type)
        {
            case DrawType::DRAW:
                text += "\\pgfusepath{draw}\n";
                break;
            
            case DrawType::FILL:
                text += "\\pgfusepath{fill}\n";
                break;
        }
    }    

    text += "\\end{pgfscope}\n";
    return text;
}


#endif // PATH2D_CPP

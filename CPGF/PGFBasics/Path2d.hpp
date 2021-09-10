#ifndef PATH2D_HPP
#define PATH2D_HPP

#include "PGFConf.hpp"
#include "Strokes2d.hpp"
#include <vector>

namespace CPGF
{
    /**
     * @brief 
     * 
     * @todo include length and area calculation.
     * 
     */
    namespace Basics
    {

        class Path2d
        {
            public:
            std::vector<SimpleStroke2d*> strokes;
            PGFConf conf;

            Path2d();
            Path2d(SimpleStroke2d& stroke,
                const PGFConf& conf = PGFConf());
            Path2d(const std::vector<SimpleStroke2d*>& strokes,
                const PGFConf& conf = PGFConf());
            Path2d(const Path2d& path);
            Path2d& operator=(const Path2d& path);
            ~Path2d();

            Path2d& translate(const AffineSpace::Vector2d& v);
            Path2d& rotate_with_respect_to(const AffineSpace::Point2d& Q, const double theta);
            Path2d& scale_with_respect_to(const AffineSpace::Point2d& Q, const AffineSpace::Vector2d& v);

            AffineSpace::Point2d& start();
            AffineSpace::Point2d start() const;
            AffineSpace::Point2d& end();
            AffineSpace::Point2d end() const;

            double length() const;
            double area() const;

            Path2d& operator+=(SimpleStroke2d& stroke);
            Path2d& add_stroke(SimpleStroke2d& stroke);
            unsigned int size() const;
            std::string render_to_string() const;
        };
    }
}

#endif // PATH2D_HPP
#ifndef OBJECT2D_HPP
#define OBJECT2D_HPP

#include <string>
#include <vector>
#include "../PGFBasics/PGFConf.hpp"
#include "../PGFBasics/Path2d.hpp"

namespace CPGF
{
    namespace Objects2d
    {
        class Object2d
        {
            public:
            std::vector<Basics::Path2d> paths;

            Object2d();
            Object2d(const Basics::Path2d& path);
            Object2d(const std::vector<Basics::Path2d>& paths);

            unsigned int size() const;

            Object2d& translate(const AffineSpace::Vector2d& v);
            Object2d& rotate_with_respect_to(const AffineSpace::Point2d& Q, const double theta);
            Object2d& scale_with_respect_to(const AffineSpace::Point2d& Q, const AffineSpace::Vector2d& s);

            friend Object2d operator+(const Object2d& A, const Object2d& B);
            Object2d& operator+=(const Object2d& B);
            Object2d& operator+=(const Basics::Path2d& path);

            std::string render_to_string() const;
        };
    }
}

#endif // OBJECT2D_HPP
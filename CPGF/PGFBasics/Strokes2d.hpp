#ifndef STROKES2D_HPP
#define STROKES2D_HPP

#include "../AffineSpace2d/Point2d.hpp"
#include <vector>

namespace CPGF
{
    namespace Basics
    {
        class SimpleStroke2d
        {
            public:
            // Clone function and destructor.
            virtual SimpleStroke2d* clone() const = 0;
            virtual ~SimpleStroke2d() = default;

            // Start and end of the stroke.
            virtual AffineSpace::Point2d& start() = 0;
            virtual AffineSpace::Point2d start() const = 0;
            virtual AffineSpace::Point2d& end() = 0;
            virtual AffineSpace::Point2d end() const = 0;

            // Invert stroke.
            // virtual SimpleStroke2d& invert() = 0;

            // Basic affine transformations.
            virtual SimpleStroke2d& translate(const AffineSpace::Vector2d& v) = 0;
            virtual SimpleStroke2d& rotate_with_respect_to(const AffineSpace::Point2d& Q, const double theta) = 0;
            virtual SimpleStroke2d& scale_with_respect_to(const AffineSpace::Point2d& Q, const AffineSpace::Vector2d& s) = 0;

            // Length and area.
            virtual double length() const = 0;
            virtual double area() const = 0;

            // Intersections.
            virtual std::vector<AffineSpace::Point2d> operator/(const SimpleStroke2d& B) = 0;

            // Return point in specific position. Tangent and normal unitary vectors.
            /**
             * @brief 
             * 
             * @param alpha a number between zero and one.
             * @return AffineSpace::Point2d 
             */
            // virtual AffineSpace::Point2d point_at_position(const double alpha) const = 0;
            // virtual AffineSpace::Vector2d tangent_vector_at_position(const double alpha) const = 0;
            // virtual AffineSpace::Vector2d normal_vector_at_position(const double alpha) const = 0;

            // // Closest point.
            // virtual double closest_position(const AffineSpace::Point2d& P) const = 0;
            // virtual AffineSpace::Point2d closest_point(const AffineSpace::Point2d& P) const = 0; 

            // Rendering to string.
            virtual std::string render_to_string() const = 0;
        };

        class StraightStroke2d;
        class BezierStroke2d;

        class StraightStroke2d: public SimpleStroke2d
        {
            public:
            std::vector<AffineSpace::Point2d> points;

            StraightStroke2d();
            StraightStroke2d(const std::vector<AffineSpace::Point2d> points);
            StraightStroke2d* clone() const override;

            StraightStroke2d& operator+=(const AffineSpace::Point2d& P);
            StraightStroke2d& add_point(const AffineSpace::Point2d& P);
            unsigned int size() const;

            AffineSpace::Point2d& start() override;
            AffineSpace::Point2d start() const override;
            AffineSpace::Point2d& end() override;
            AffineSpace::Point2d end() const override;
            StraightStroke2d& translate(const AffineSpace::Vector2d& v) override;
            StraightStroke2d& rotate_with_respect_to(const AffineSpace::Point2d& Q, const double theta) override;
            StraightStroke2d& scale_with_respect_to(const AffineSpace::Point2d& Q, const AffineSpace::Vector2d& s) override;

            double length() const override;
            double area() const override;

            std::vector<AffineSpace::Point2d> operator/(const SimpleStroke2d& B) override;
            std::vector<AffineSpace::Point2d> operator/(const StraightStroke2d& B);
            // std::vector<AffineSpace::Point2d> operator/(const BezierStroke2d& B);

            std::string render_to_string() const override;
        };

        /**
         * @brief 
         * 
         * @todo Finish implementation.
         */
        class BezierStroke2d: public SimpleStroke2d
        {
            public:
            // Start and end points.
            /*! Starting point.*/
            AffineSpace::Point2d P1;
            /*! Ending point*/
            AffineSpace::Point2d P2;

            // Control Points.
            /*! First control point.*/
            AffineSpace::Point2d Q1;
            /*! Second control point.*/
            AffineSpace::Point2d Q2;

            BezierStroke2d();
            BezierStroke2d(const AffineSpace::Point2d& P1, const AffineSpace::Point2d& P2,
                const AffineSpace::Point2d& Q1, const AffineSpace::Point2d& Q2);
            
            BezierStroke2d* clone() const override;

            AffineSpace::Point2d& start() override;
            AffineSpace::Point2d start() const override;
            AffineSpace::Point2d& end() override;
            AffineSpace::Point2d end() const override;
            BezierStroke2d& translate(const AffineSpace::Vector2d& v) override;
            BezierStroke2d& rotate_with_respect_to(const AffineSpace::Point2d& Q, const double theta) override;
            BezierStroke2d& scale_with_respect_to(const AffineSpace::Point2d& Q, const AffineSpace::Vector2d& s) override;

            double length() const override;
            double area() const override;

            std::vector<AffineSpace::Point2d> operator/(const SimpleStroke2d& B) override;
            // std::vector<AffineSpace::Point2d> operator/(const StraightStroke2d& B);
            // std::vector<AffineSpace::Point2d> operator/(const BezierStroke2d& B);

            std::string render_to_string() const override;
        };
    }
}

#endif // STROKES2D_HPP
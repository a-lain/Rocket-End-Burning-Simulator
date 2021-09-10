/** \file This file contains all operations related to points.
 */

#ifndef POINT2D_HPP
#define POINT2D_HPP

#include "Vector2d.hpp"

namespace CPGF
{
    namespace AffineSpace
    {
        /**
         * @brief A point of a 2D affine space.
         * 
         * As usual, a vector may be added to a point and two points may be "substracted"
         * to obtain a vector.
         * 
         */
        class Point2d
        {
            public:
            /**
             * @brief The x coordinate of the point. 
             * 
             */
            double x;
            /**
             * @brief The y coordinate of the point.
             * 
             */
            double y;

            /**
             * @brief Function that adds a vector v to the point P.
             * 
             * Returns the position of the end of the vector v when its start is
             * placed at the point P.
             * 
             * @param P a 2D point.
             * @param v a 2D vector.
             * @return Point2D
             */
            friend Point2d operator+(const Point2d& P, const Vector2d& v);

            /**
             * @brief Function that returns the vector that joins two points.
             * 
             * Returns the vector whose start is at point P and whose end is at
             * point Q.
             * 
             * @param P a 2D point.
             * @param Q a 2D point.
             * @return Vector2d 
             */
            friend Vector2d operator-(const Point2d& P, const Point2d& Q);

            /**
             * @brief Returns whether P and Q are closer than 1e-10.
             * 
             * @param P 
             * @param Q 
             * @return true 
             * @return false 
             */
            friend bool operator==(const Point2d& P, const Point2d& Q);

            /**
             * @brief Traslates the point by the vector v.
             * 
             * Equivalent to *this = *this + v;
             * 
             * @param v a 2D vector.
             * @return Point2d& 
             */
            Point2d& operator+= (const Vector2d& v);

            /**
             * @brief Calculates the angle the point Q has with respect to P.
             * 
             * Returns the angle between the horizontal line that passes through P and the
             * line that joins P and Q.
             * 
             * @param Q a 2D point.
             * @return double. A real number between 0 and \f$2\pi\f$.
             */
            double angle_with_respect_to(const Point2d& Q);

            /**
             * @brief Rotates the point (*this) an angle theta around the point Q.
             * 
             * Using Q as the center of rotation, the line that joins P and Q is
             * rotates around Q an angle theta. As a consequence, the position of P
             * changes, although the distance PQ is preserved.
             * 
             * @param Q a 2D point.
             * @param theta a double. A real number between 0 and \f$2\pi\f$.
             * @return Point2d& 
             */
            Point2d& rotate_with_respect_to(const Point2d& Q, const double theta);

            /**
             * @brief Rotates the point (*this) to a fixed angle theta around the
             * point Q.
             * 
             * Using Q as the center of rotation, the line that joins P and Q is
             * rotated around Q until the angle between that line
             * and the horizontal line that passes through Q is exactly theta. As
             * a consequence, the position of P changes, although the distance PQ
             * is preserved.
             * 
             * @param Q a 2D point.
             * @param theta a double. A real number between 0 and \f$2\pi\f$.
             * @return Point2d& 
             */
            Point2d& rotate_to_with_respect_to(const Point2d& Q, const double theta);

            /**
             * @brief The position of *this with respect to the point Q is scaled
             * according to the components of the vector s.
             * 
             * First, we express the point P through its coordinates with respect to
             * the point Q. Then, those coordinates are multiplied componentwise by
             * the vector s. Finally, the new P is expressed through its coordinates
             * with respect to the origin.
             * 
             * @param Q 
             * @param s 
             * @return Point2d& 
             */
            Point2d& scale_with_respect_to(const Point2d& Q, const Vector2d& s);

            /**
             * @brief Retruns a string representation of the point.
             * 
             * @return std::string 
             */
            std::string to_string() const;

            /**
             * @brief Returns the origin of coordinates.
             * 
             */
            Point2d();

            /**
             * @brief Retruns the point (x,y).
             * 
             * @param x 
             * @param y 
             */
            Point2d(double x, double y);
        };
    }
    
}


#endif // POINT2D_HPP
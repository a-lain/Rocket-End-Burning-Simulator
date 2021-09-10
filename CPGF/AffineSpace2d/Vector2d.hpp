#ifndef VECTOR2D_HPP
#define VECTOR2D_HPP

#include <string>

/**
 * @brief Every component of this library is part of the CTikZ namespace.
 * 
 */
namespace CPGF
{
    /**
     * @brief This namespace contains all objects related to a mathematical
     * affine space, i.e., vectors and points. 
     */
    namespace AffineSpace
    {
        /**
         * @brief An object that represents a 2D vector.
         * 
         * All usual operations +,-,*,/ are defined component wise. The |
         * operator is used for the euclidean scalar product.
         */
        class Vector2d
        {
            public:
            /**
             * @brief The x component of the vector.
             * 
             */
            double x;

            /**
             * @brief The y component of the vector.
             * 
             */
            double y;

            /**
             * @brief Returns (v.x + w.x, v.y + w.y).
             * 
             * @param v a 2D vector.
             * @param w a 2D vector.
             * @return Vector2d 
             */
            friend Vector2d operator+(const Vector2d& v, const Vector2d& w);

            /**
             * @brief Retruns (v.x - w.x, v.y - w.y).
             * 
             * @param v a 2D vector.
             * @param w a 2D vector.
             * @return Vector2d 
             */
            friend Vector2d operator-(const Vector2d& v, const Vector2d& w);

            /**
             * @brief Returns (v.x * w.x, v.y * w.y).
             * 
             * @param v a 2D vector.
             * @param w a 2D vector.
             * @return Vector2d 
             */
            friend Vector2d operator*(const Vector2d& v, const Vector2d& w);

            /**
             * @brief Retruns (v.x / w.x, v.y / w.y).
             * 
             * @param v a 2D vector.
             * @param w a 2D vector.
             * @return Vector2d 
             */
            friend Vector2d operator/(const Vector2d& v, const Vector2d& w);

            /**
             * @brief Retruns the euclidean scalar product of the two vectors.
             * 
             * @param v a 2D vector.
             * @param w a 2D vector.
             * @return double 
             */
            friend double operator|(const Vector2d& v, const Vector2d& w);

            /**
             * @brief Returns itself.
             * 
             * @return Vector2d& 
             */
            Vector2d& operator+();

            /**
             * @brief Multiplies all components by -1.
             * 
             * @return Vector2d& 
             */
            Vector2d& operator-();

            /**
             * @brief Stores in *this the sum of *this and w.
             * 
             * Equivalent to *this = *this + w;
             * 
             * @param w a 2D vector.
             * @return Vector2d& 
             */
            Vector2d& operator+= (const Vector2d& w);

            /**
             * @brief Stores in *this the substraction of *this and w.
             * 
             * Equivalent to *this = *this - w;
             * 
             * @param w a 2D vector.
             * @return Vector2d& 
             */
            Vector2d& operator-= (const Vector2d& w);

            /**
             * @brief Stores in *this the multiplication of *this and w.
             * 
             * Equivalent to *this = *this * w;
             * 
             * @param w a 2D vector.
             * @return Vector2d& 
             */
            Vector2d& operator*= (const Vector2d& w);

            /**
             * @brief Stores in *this the division of *this and w.
             * 
             * Equivalent to *this = *this / w;
             * 
             * @param w 
             * @return Vector2d& 
             */
            Vector2d& operator/= (const Vector2d& w);

            /**
             * @brief Returns a perpendicular vector.
             * 
             * The returned vector is always positioned 90° degrees
             * anticlockwise.
             * 
             */
            Vector2d perp() const;

            /**
             * @brief Returns the euclidean norm of the vector.
             * 
             * @return double 
             */
            double norm() const;

            /**
             * @brief Retruns a string representation of the vector.
             * 
             * @return std::string 
             */
            std::string to_string() const;

            /**
             * @brief Returns the null vector (0,0).
             * 
             */
            Vector2d();

            /**
             * @brief Returns the vector (x,x).
             * 
             * @param x a real number.
             */
            Vector2d(double x);

            /**
             * @brief Returns the vector (x,y).
             * 
             * @param x a real number.
             * @param y a real number.
             */
            Vector2d(double x, double y);
        };
    }
}



#endif // VECTOR2D_HPP
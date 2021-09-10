#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>

namespace Math
{
    namespace Interpolation
    {
        class AverageLinearInterpolation
        {
            public:

            AverageLinearInterpolation(const std::vector<double>& averages,
                const std::vector<double>& partition);
            AverageLinearInterpolation(double* averages,
                double *partition, unsigned int N_cells);

            AverageLinearInterpolation(const AverageLinearInterpolation& f);
            AverageLinearInterpolation& operator=(const AverageLinearInterpolation& f);
            ~AverageLinearInterpolation();

            /**
             * @brief 
             * 
             * We do linear extrapolation at the ends.
             * 
             * @param x 
             * @return double 
             */
            double operator()(const double x) const;

            public:
            unsigned int N_cells;
            double* x_part;
            double* b;
            double* c;
        };

        class AverageQuadraticInterpolation
        {
            public:

            AverageQuadraticInterpolation(const std::vector<double>& averages,
                const std::vector<double>& partition);
            AverageQuadraticInterpolation(double* averages,
                double *partition, unsigned int N_cells);

            AverageQuadraticInterpolation(const AverageQuadraticInterpolation& f);
            AverageQuadraticInterpolation& operator=(const AverageQuadraticInterpolation& f);
            ~AverageQuadraticInterpolation();

            /**
             * @brief 
             * 
             * We do linear extrapolation at the ends.
             * 
             * @param x 
             * @return double 
             */
            double operator()(const double x) const;

            protected:
            unsigned int N_cells;
            double* x_part;
            double* a;
            double* b;
            double* c;
        };

        class LinearInterpolation
        {
            public:
            LinearInterpolation(const std::vector<double>& Y, const std::vector<double>& X);
            LinearInterpolation(const double* Y, const double* X, const unsigned int N);

            LinearInterpolation(const LinearInterpolation& I);
            LinearInterpolation& operator=(const LinearInterpolation& I);
            ~LinearInterpolation();

            double operator()(const double x) const;

            protected:
            unsigned int N;
            double* x_part;
            double* b;
            double* c;
        };
    }
}

#endif // INTERPOLATION_HPP
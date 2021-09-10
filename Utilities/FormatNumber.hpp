#ifndef FORMAT_NUMBER_HPP
#define FORMAT_NUMBER_HPP

#include <functional>
#include <string>

namespace Utilities
{
    int position_of_most_significant_digit(const double x);

    unsigned long long pow10(const int n);

    double round_to_precision(double x, const unsigned int precision);
    double ceil_to_precision(double x, const unsigned int precision);
    double floor_to_precision(double x, const unsigned int precision);

    std::string format_number(double x, const bool latex = true,
        const unsigned int precision = 3, const bool show_sign = false,
        const int lim_inf = -3, const int lim_sup = 3);
}

#endif // FORMAT_NUMBER_HPP
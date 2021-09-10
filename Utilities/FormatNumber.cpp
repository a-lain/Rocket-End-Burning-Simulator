#ifndef FORMAT_NUMBER_CPP
#define FORMAT_NUMBER_CPP

#include "FormatNumber.hpp"
#include <cmath>

using namespace Utilities;

int Utilities::position_of_most_significant_digit(const double x)
{
    if (x == 0)
    {
        return 0;
    }
    else
    {
        return (int)floor(log10(fabs(x)*(1 + 1e-14)));
    }
}

unsigned long long Utilities::pow10(const int n)
{
    static unsigned long long pow10_array[19] =
    {
        1ll, 10ll, 100ll, 1000ll, 10000ll, 
        100000ll, 1000000ll, 10000000ll, 100000000ll, 1000000000ll,
        10000000000ll, 100000000000ll, 1000000000000ll, 10000000000000ll, 100000000000000ll,
        1000000000000000ll, 10000000000000000ll, 100000000000000000ll, 1000000000000000000ll
    };

    return pow10_array[n]; 
}

double Utilities::round_to_precision(double x, const unsigned int precision)
{
    if (precision == 0)
    {
        return 0;
    }
    else
    {
        int exponent = (int) floor((x != 0) ? log10(fabs(x)*(1 + 1./pow10(precision+1))) : 0);
        return (exponent - (int)precision + 1 >= 0) ?
            round(x / pow10(exponent - precision + 1)) * pow10(exponent - precision + 1)
            : round(x * pow10(precision - exponent - 1)) / pow10(precision - exponent - 1);
    }
}

double Utilities::ceil_to_precision(double x, const unsigned int precision)
{
    if (precision == 0)
    {
        return 0;
    }
    else
    {
        int exponent = (int) floor((x != 0) ? log10(fabs(x)*(1 + 1./pow10(precision+1))) : 0);
        return (exponent - (int)precision + 1 >= 0) ?
            ceil(x / pow10(exponent - precision + 1)) * pow10(exponent - precision + 1)
            : ceil(x * pow10(precision - exponent - 1)) / pow10(precision - exponent - 1);
    }
}

double Utilities::floor_to_precision(double x, const unsigned int precision)
{
    if (precision == 0)
    {
        return 0;
    }
    else
    {
        int exponent = (int) floor((x != 0) ? log10(fabs(x)*(1 + 1./pow10(precision+1))) : 0);
        return (exponent - (int)precision + 1 >= 0) ?
            floor(x / pow10(exponent - precision + 1)) * pow10(exponent - precision + 1)
            : floor(x * pow10(precision - exponent - 1)) / pow10(precision - exponent - 1);
    }
}

std::string Utilities::format_number(const double x, const bool latex,
    const unsigned int precision, const bool show_sign, const int lim_inf, const int lim_sup)
{
    if (precision == 0)
    {
        return std::string("0");
    }

    else
    {
        // We calculate the decimal logarithm of 'x', unless it's zero. In that case,
        // we impose log10(0):=0. We add a small number to faciliate rounding up operations.
        double log10_x = (x != 0) ? log10(fabs(x)*(1 + 1./pow10(precision+1))) : 0;

        // If we apply the floor function to 'log10_x', we will obtain the exponent needed
        // for the number's scientific notation. Please realize that the exponent is
        // computed even if the conditions for scientific notation display are not met;
        // this calculation will be needed later. Explanation: let our number 'x' be
        // greater than 1, for example between 10^3 and 10^4. Then the required 
        // exponent for the number's scientific notation is 3. On the contrary, if our
        // number 'x' is less than 1, for instance between 10^(-3) and 10^(-4), the
        // searched exponent will be -4. Therefore, the exponent will always be the
        // floor function of the 'log10_x'.
        int exponent = (int) floor(log10_x);

        // We always round off the number to the wanted precision first.
        unsigned long long base = (exponent - (int)precision + 1 >= 0) ?
            (unsigned long long)round(fabs(x) / pow10(exponent - precision + 1))
            : (unsigned long long)round(fabs(x) * pow10(precision - exponent - 1));

        // We take into account the sign stuff.
        std::string sign = "";
        if (x > 0 && show_sign)
        {
            sign = "+";
        }
        else if (x < 0)
        {
            sign = "-";
        }

        // If the conditions for scientific notation display are satisfied.
        if (log10_x >= lim_sup || log10_x <= lim_inf)
        {
            // The special case of a power of ten.
            if (base == pow10(precision - 1))
            {
                return sign + ( (latex) ?
                "10^{" + std::to_string(exponent) + "}" : "10^(" + std::to_string(exponent) + ")");
            }
            else
            {
                std::string end = (latex) ?
                    "\\cdot 10^{" + std::to_string(exponent) + "}" : "*10^(" + std::to_string(exponent) + ")";
                return sign + std::to_string(base).insert(1, ".") + end;     
            }        
        }

        else
        {
            if (base == 0)
            {
                return "0." + std::string(precision - 1, '0');
            }
            else
            {
                if (exponent >= 0)
                {
                    if ((unsigned int)exponent + 1 >= precision)
                    {
                        return sign + std::to_string(base * pow10(exponent + 1 - precision));
                    }
                    else
                    {
                        return sign + std::to_string(base).insert(exponent + 1, ".");
                    }
                }
                else
                {
                    return sign + "0." + std::string(-exponent - 1, '0') + std::to_string(base);
                }
            }
        }
    }    
}

#endif // FORMAT_NUMBER_CPP
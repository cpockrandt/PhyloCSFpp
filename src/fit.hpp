#pragma once

#include <math.h>
#include <cassert>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

template <typename function_t>
auto fit_find_init(const uint32_t max_tries, const double init, const double lo, const double hi, function_t f)
{
    assert(lo < hi && lo > 0.0);

    srand(0);

    const double width = log(hi) - log(lo);
    const auto flo = f(lo);
    const auto fhi = f(hi);
    double x = init;
    auto fx = f(init);
    std::cout << "xxx: " << x << " ||| " << fx << '\n';

    uint32_t i = 0;
    while (i < max_tries && (std::get<1>(fx) <= std::get<1>(flo) || std::get<1>(fx) <= std::get<1>(fhi)))
    {
        const double r = width * ((float) rand()/RAND_MAX);
        x = exp(log(lo) + r);
        fx = f(x);
        std::cout << "xxx: " << x << " ||| " << fx << '\n';
        ++i;
    }
    std::cout << "------------------------------------------------\n";
    if (i == max_tries)
        return (flo > fhi) ? flo : fhi;
    return fx;
}

//template <typename function_t>
//auto fit_maximize(const double lo, const double hi, const std::tuple<double, double, double> & good_init, function_t f)
//{
//
//}

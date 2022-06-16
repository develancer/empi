/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cmath>
#include "GaussianFamily.h"

static const double NORM = std::sqrt(M_SQRT2);

//////////////////////////////////////////////////////////////////////////////

GaussianFamily::GaussianFamily(double min_max) : min_max(min_max) {
    assert(min_max > 0);
}

double GaussianFamily::max_arg() const {
    return min_max;
}

double GaussianFamily::min_arg() const {
    return -min_max;
}

const char *GaussianFamily::name() const {
    return "gauss";
}

double GaussianFamily::value(double t) const {
    return NORM * std::exp(-M_PI * t * t);
}

double GaussianFamily::scale_integral(double log_scale) const {
    return 1.0 / std::sqrt(std::cosh(log_scale));
}

double GaussianFamily::inv_scale_integral(double value) const {
    return std::acosh(1 / (value * value));
}

double GaussianFamily::freq_integral(double x) const {
    return std::exp(-M_PI_2 * x * x);
}

double GaussianFamily::inv_freq_integral(double value) const {
    return std::sqrt(-M_2_PI * std::log(value));
}

double GaussianFamily::skew_integral(double) const {
    return 0.0;
}

double GaussianFamily::time_integral(double x) const {
    return std::exp(-M_PI_2 * x * x);
}

double GaussianFamily::inv_time_integral(double value) const {
    return std::sqrt(-M_2_PI * std::log(value));
}

double GaussianFamily::optimality_factor_e2(double epsilon2) const {
    return 1 - 1.5 * epsilon2;
}

double GaussianFamily::optimality_factor_sf(double sf) const {
    return 1 - exp(-1.59 * sf - 2.11);
}

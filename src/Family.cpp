/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cmath>
#include "Family.h"

static const double NORM = std::sqrt(M_SQRT2);

//////////////////////////////////////////////////////////////////////////////

GaussianFamily::GaussianFamily(double min_max) : min_max(min_max) {
    assert(min_max > 0);
}

double GaussianFamily::max_arg() {
    return min_max;
}

double GaussianFamily::min_arg() {
    return -min_max;
}

double GaussianFamily::value(double t) {
    return NORM * std::exp(-M_PI * t * t);
}

double GaussianFamily::scale_integral(double log_scale) {
    return 1.0 / std::sqrt(std::cosh(log_scale));
}

double GaussianFamily::inv_scale_integral(double value) {
    return std::acosh(1 / (value * value));
}

double GaussianFamily::freq_integral(double x) {
    return std::exp(-M_PI_2 * x * x);
}

double GaussianFamily::inv_freq_integral(double value) {
    return std::sqrt(-M_2_PI * std::log(value));
}

double GaussianFamily::skew_integral(double) {
    return 0.0;
}

double GaussianFamily::time_integral(double x) {
    return std::exp(-M_PI_2 * x * x);
}

double GaussianFamily::inv_time_integral(double value) {
    return std::sqrt(-M_2_PI * std::log(value));
}

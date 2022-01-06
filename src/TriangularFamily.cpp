/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cmath>
#include "TriangularFamily.h"

static const double MIN_MAX = std::sqrt(2.5*M_1_PI);
static const double NORM = std::pow(0.9 * M_PI, 0.25);

//////////////////////////////////////////////////////////////////////////////

TriangularFamily::TriangularFamily() {
}

double TriangularFamily::max_arg() const {
    return MIN_MAX;
}

double TriangularFamily::min_arg() const {
    return -MIN_MAX;
}

const char *TriangularFamily::name() const {
    return "triangular";
}

double TriangularFamily::value(double t) const {
    const double abs_t = std::abs(t);
    return (abs_t < MIN_MAX) ? NORM * (1.0 - abs_t / MIN_MAX) : 0.0;
}

double TriangularFamily::scale_integral(double log_scale) const {
    const double exp_half_log_scale = std::exp(-0.5 * log_scale);
    return (3.0 - exp_half_log_scale * exp_half_log_scale) * exp_half_log_scale / 2;
}

double TriangularFamily::inv_scale_integral(double value) const {
    return solve_integral(&TriangularFamily::scale_integral, value);
}

double TriangularFamily::freq_integral(double x) const {
    if (std::abs(x) < 0.001) {
        // approximation to avoid numerical errors
        return 1.0 - M_PI_2 * x * x;
    }
    const double x_rel = std::sqrt(10 * M_PI) * x;
    return 6 / (x_rel * x_rel) * (1 - sin(x_rel) / x_rel);
}

double TriangularFamily::inv_freq_integral(double value) const {
    return solve_integral(&TriangularFamily::freq_integral, value);
}

double TriangularFamily::skew_integral(double) const {
    return 0.0;
}

double TriangularFamily::time_integral(double x) const {
    const double x_rel = std::abs(x) / MIN_MAX;
    if (x_rel <= 1.0) {
        return 1 + 0.75 * x_rel * x_rel * (x_rel - 2.0);
    }
    if (x_rel <= 2.0) {
        const double t = 1.0 - 0.5 * x_rel;
        return 2 * t * t * t;
    }
    return 0.0;
}

double TriangularFamily::inv_time_integral(double value) const {
    return solve_integral(&TriangularFamily::time_integral, value);
}

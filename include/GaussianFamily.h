/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_GAUSSIAN_FAMILY_H
#define EMPI_GAUSSIAN_FAMILY_H

#include "Family.h"

/**
 * Family implementation for Gaussian envelope (Gabor atoms).
 */
class GaussianFamily : public FamilyTemplate<GaussianFamily> {
    const double min_max;

public:
    /**
     * Create a new instance representing a Gaussian envelope.
     *
     * @param min_max half-width of envelope function i.e. for |t| larger than this value, f(t) will be assumed as zero
     */
    explicit GaussianFamily(double min_max = 3.0);

    double max_arg() const final;

    double min_arg() const final;

    double value(double t) const final;

    double scale_integral(double log_scale) const final;

    double inv_scale_integral(double value) const final;

    double freq_integral(double x) const final;

    double inv_freq_integral(double value) const final;

    double skew_integral(double x) const final;

    double time_integral(double x) const final;

    double inv_time_integral(double value) const final;
};

#endif //EMPI_GAUSSIAN_FAMILY_H

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_TRIANGULAR_FAMILY_H
#define EMPI_TRIANGULAR_FAMILY_H

#include "Family.h"

/**
 * Family implementation for triangular envelope.
 */
class TriangularFamily : public FamilyTemplate<TriangularFamily> {
public:
    /**
     * Create a new instance representing a triangular envelope.
     */
    explicit TriangularFamily() = default;

    double max_arg() const final;

    double min_arg() const final;

    const char *name() const final;

    double value(double t) const final;

    double scale_integral(double log_scale) const final;

    double inv_scale_integral(double value) const final;

    double freq_integral(double x) const final;

    double inv_freq_integral(double value) const final;

    double skew_integral(double x) const final;

    double time_integral(double x) const final;

    double inv_time_integral(double value) const final;
};

#endif //EMPI_TRIANGULAR_FAMILY_H

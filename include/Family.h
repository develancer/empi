/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_FAMILY_H
#define EMPI_FAMILY_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include "Types.h"

/**
 * Base class for all families of envelope functions.
 * Each family implementation should provide not only the function's values,
 * but also a set of integrals needed to construct an optimal dictionary.
 * The only direct descendant is FamilyTemplate, using curiously recurring template pattern.
 *
 * All integrals (A, B, C etc.) used in this class are defined in https://doi.org/10.1049/iet-spr.2019.0246
 */
class Family {
public:
    /**
     * Calculate number of samples needed to store all non-zero values of the realization of this envelope function
     * with a given scale, centered at a given position.
     *
     * @param center_position in samples
     * @param scale in samples
     * @param offset pointer to be optionally filled with the index of the first sample
     * (can be negative if center_position is zero or negative)
     * @return number of samples for the envelope realization
     */
    virtual index_t size_for_values(double center_position, double scale, index_t *offset) = 0;

    /**
     * Fill the array with all non-zero values of the realization of this envelope function
     * with a given scale, centered at a given position.
     *
     * @param center_position in samples
     * @param scale in samples
     * @param offset pointer to additionally return the index of the first sample
     * (can be negative if center_position is zero or negative)
     * @param values array to be filled with values;
     * in particular: values[i] will be filled with this->value((offset-center_position+i)/scale)
     * @param normalize true if values should be L²-normalized (so the sum of squares equals 1 exactly), false otherwise
     * @return value of normalization factor that was used (or would be, if normalize=false) for normalization
     */
    virtual double generate_values(double center_position, double scale, index_t *offset, double *values, bool normalize) = 0;

    /**
     * @return the maximal value of t for which this->value(t) returns a non-negligible (non-zero) value.
     */
    virtual double max_arg() = 0;

    /**
     * @return the minimal value of t for which this->value(t) returns a non-negligible (non-zero) value.
     */
    virtual double min_arg() = 0;

    /**
     * Return the value of the envelope function for any given argument.
     * This is the most important method of the Family classes.
     * All Family implementations should be L²-normalized, fulfilling ∫ f(t)² dt = 1
     * as well as ∫ t² f(t)² dt = 1/4π.
     *
     * @param t argument (time)
     * @return value of this envelope function
     */
    virtual double value(double t) = 0;

    /**
     * Return the value of integral A(Δλ) = ∫ f(t e^(Δλ/2)) f(t e^(−Δλ/2)) dt
     *
     * @param log_scale Δλ
     * @return value of this improper integral over ±∞
     */
    virtual double scale_integral(double log_scale) = 0;

    /**
     * Return the value of the function inverse to integral A(Δλ).
     *
     * @param value fulfilling A(Δλ) = value
     * @return Δλ
     */
    virtual double inv_scale_integral(double value) = 0;

    /**
     * Return the value of integral B(x) = ∫ f(t)² cos(2πxt) dt
     *
     * @param x frequency-related parameter
     * @return value of this improper integral over ±∞
     */
    virtual double freq_integral(double x) = 0;

    /**
     * Return the value of the function inverse to integral B(x).
     *
     * @param value fulfilling B(x) = value
     * @return x
     */
    virtual double inv_freq_integral(double value) = 0;

    /**
     * Return the value of integral S(x) = ∫ f(t)² sin(2πxt) dt
     *
     * @param x frequency-related parameter
     * @return value of this improper integral over ±∞
     */
    virtual double skew_integral(double x) = 0;

    /**
     * Return the value of integral C(x) = ∫ f(t+x/2) f(t−x/2) dt
     *
     * @param x time-related parameter
     * @return value of this improper integral over ±∞
     */
    virtual double time_integral(double x) = 0;

    /**
     * Return the value of the function inverse to integral C(x).
     *
     * @param value fulfilling C(x) = value
     * @return x
     */
    virtual double inv_time_integral(double value) = 0;

    virtual ~Family() = default;
};

/**
 * Direct ancestor of all implementations of envelope function families.
 * The curiously recurring template pattern (CRTP) is used here.
 * @tparam E concrete Family implementation class
 */
template<class E>
class FamilyTemplate : public Family {
public:
    /**
     * Calculate number of samples needed to store all non-zero values of the realization of this envelope function
     * with a given scale, centered at a given position.
     *
     * @param center_position in samples
     * @param scale in samples
     * @param offset pointer to be optionally filled with the index of the first sample
     * (can be negative if center_position is zero or negative)
     * @return number of samples for the envelope realization
     */
    index_t size_for_values(double center_position, double scale, index_t *offset) final {
        E *self = static_cast<E *>(this);
        assert(self);

        index_t first_sample_offset = Types::round<index_t>(center_position + scale * self->min_arg() + 0.5);
        index_t last_sample_offset = Types::round<index_t>(center_position + scale * self->max_arg() - 0.5);
        index_t sample_count = last_sample_offset + 1 - first_sample_offset;

        if (offset) {
            *offset = first_sample_offset;
        }
        return std::max<index_t>(0, sample_count);
    }

    /**
     * Fill the array with all non-zero values of the realization of this envelope function
     * with a given scale, centered at a given position.
     *
     * @param center_position in samples
     * @param scale in samples
     * @param offset pointer to additionally return the index of the first sample
     * (can be negative if center_position is zero or negative)
     * @param values array to be filled with values;
     * in particular: values[i] will be filled with this->value((offset-center_position+i)/scale)
     * @param normalize true if values should be L²-normalized (so the sum of squares equals 1 exactly), false otherwise
     * @return value of normalization factor that was used (or would be, if normalize=false) for normalization
     */
    double generate_values(double center_position, double scale, index_t *offset, double *values, bool normalize) final {
        E *self = static_cast<E *>(this);
        assert(self);

        index_t first_sample_offset;
        index_t sample_count = size_for_values(center_position, scale, &first_sample_offset);
        if (offset) {
            *offset = first_sample_offset;
        }

        double sum2 = 0.0;
        for (index_t i = 0; i < sample_count; ++i) {
            const double value = self->value(static_cast<double>(first_sample_offset + i - center_position) / scale);
            values[i] = value;
            sum2 += value * value;
        }

        const double norm = 1.0 / std::sqrt(sum2);
        if (normalize) {
            for (index_t i = 0; i < sample_count; ++i) {
                values[i] *= norm;
            }
            return norm;
        }

        return norm;
    }
};

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

    double max_arg() final;

    double min_arg() final;

    double value(double t) final;

    double scale_integral(double log_scale) final;

    double inv_scale_integral(double value) final;

    double freq_integral(double x) final;

    double inv_freq_integral(double value) final;

    double skew_integral(double x) final;

    double time_integral(double x) final;

    double inv_time_integral(double value) final;
};

#endif //EMPI_FAMILY_H

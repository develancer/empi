/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_ATOM_BASE_H
#define EMPI_BLOCK_ATOM_BASE_H

#include <array>
#include <cmath>
#include <utility>

/**
 * Simple structure with three fields, used in both BlockAtom and ExtendedBlockAtom.
 */
struct BlockAtomParams {
    double frequency; // between 0 and 0.5 (1 would be sampling frequency)
    double position; // in samples
    double scale; // in samples

    BlockAtomParams(double frequency, double position, double scale)
            : frequency(frequency), position(position), scale(scale) {}
};

/**
 * Auxiliary class converting atom parameters between actual values and scaled values used internally by the minimizer.
 */
class BlockAtomParamsConverter {
public:
    const double frequency_step;
    const double position_step;
    const double log_scale_step;
    const double max_scaled_frequency;

    BlockAtomParamsConverter() : frequency_step(1.0), position_step(1.0), log_scale_step(1.0), max_scaled_frequency(0.5) { }

    BlockAtomParamsConverter(double frequency_step, double position_step, double log_scale_step, double frequency_max);

    /**
     * Convert from actual values to scaled values used internally by the minimizer.
     * @param params in actual units (@see BlockAtomParams)
     * @return array of 3 parameters in unit related to dictionary steps
     */
    [[nodiscard]] std::array<double, 3> arrayFromParams(const BlockAtomParams &params) const;

    /**
     * Convert from scaled values used internally by the minimizer to actual values.
     * @param array of 3 parameters in unit related to dictionary steps
     * @return params in actual units (@see BlockAtomParams)
     */
    [[nodiscard]] std::pair<BlockAtomParams, double> paramsFromArray(const std::array<double, 3> &array) const;

private:
    double fix_frequency(double &scaled_frequency) const;

protected:
    virtual double fix_log_scale(double &scaled_log_scale) const;

    /**
     * Fix the value given as a reference so it will be between min and max, inclusive.
     *
     * @return if given value has already been between min and max, 0;
     * otherwise, its distance from this interval
     */
    static double fix_scaled_argument(double &value, double min, double max);
};

/**
 * Auxiliary class converting atom parameters between actual values and scaled values used internally by the minimizer.
 */
class BlockAtomParamsConverterBounded : public BlockAtomParamsConverter {
    double min_scaled_log_scale;
    double max_scaled_log_scale;

public:
    BlockAtomParamsConverterBounded(double frequency_step, double position_step, double log_scale_step, double frequency_max, double min_scale, double max_scale);

    double fix_log_scale(double &scaled_log_scale) const final;
};

#endif //EMPI_BLOCK_ATOM_BASE_H

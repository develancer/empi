/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_EXTRA_DATA_H
#define EMPI_EXTRA_DATA_H

/**
 * Plain data structure for channel-related additional atom data.
 */
struct ExtraData {
    /** atom's amplitude, must be multiplied by family->value(0.0) / std::sqrt(scale) */
    double amplitude;

    /** atom's energy, computed as sum of all samples squared */
    double energy;

    /** atom's phase between ±π */
    double phase;
};

#endif //EMPI_EXTRA_DATA_H

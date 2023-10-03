/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_EXTRACTOR_H
#define EMPI_EXTRACTOR_H

#include "ExtraData.h"
#include "Corrector.h"
#include "Types.h"

/**
 * Plain data structure storing a frequency index (bin) of the best-match atom as well as its energy.
 */
struct ExtractedMaximum {
    /** atom's energy, computed as sum of all samples squared */
    double energy;

    /** index of the bin in discrete Fourier transform where 0 corresponds to zero frequency and so forth */
    int bin_index;

    bool operator<(const ExtractedMaximum &other) const {
        return energy < other.energy;
    }
};

using Extractor = ExtractedMaximum (*)(int channel_count, int output_bins, complex *const *channels, const Corrector *correctors,
                                       double *bins_buffer, ExtraData *atom_data);

/**
 * Extractor implementation for single-channel decomposition.
 */
ExtractedMaximum extractorSingleChannel(int channel_count, int output_bins, complex *const *channels, const Corrector *correctors,
                                        double *bins_buffer, ExtraData *atom_data);

/**
 * Extractor implementation for multi-channel decomposition where atoms in all signals share position, frequency and scale,
 * but amplitudes as well as phases may differ across channels.
 */
ExtractedMaximum extractorVariablePhase(int channel_count, int output_bins, complex *const *channels, const Corrector *correctors,
                                        double *bins_buffer, ExtraData *atom_data);

/**
 * Extractor implementation for multi-channel decomposition where atoms in all signals share position, frequency, scale as well as phase,
 * but amplitudes may differ across channels.
 */
ExtractedMaximum extractorConstantPhase(int channel_count, int output_bins, complex *const *channels, const Corrector *correctors,
                                        double *bins_buffer, ExtraData *atom_data);

#endif //EMPI_EXTRACTOR_H

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cassert>
#include "Extractor.h"

ExtractedMaximum extractorSingleChannel(int channel_count, int output_bins, const complex *const *channels, const Corrector *correctors,
                                        double *bins_buffer, ExtraData *extra_data) {
    assert(channel_count == 1);
    const complex *channel = channels[0];
    double max_upper_bound = std::numeric_limits<double>::lowest();
    int i_max = 0;
    for (int i = 0; i < output_bins; ++i) {
        double energy_upper_bound = correctors[i].estimate_energy(channel[i]);
        bins_buffer[i] = energy_upper_bound;
        if (energy_upper_bound > max_upper_bound) {
            max_upper_bound = energy_upper_bound;
            i_max = i;
        }
    }
    double energy_guess = correctors[i_max].compute(channels[0][i_max]).energy();
    bins_buffer[i_max] = energy_guess;
    ExtractedMaximum max{energy_guess, i_max};
    for (int i = 0; i < output_bins; ++i)
        if (bins_buffer[i] >= max.energy) {
            CorrectorResult result = correctors[i].compute(channel[i]);
            double energy = result.energy();
            if (extra_data) {
                extra_data[0].amplitude = result.amplitude();
                extra_data[0].energy = energy;
                extra_data[0].phase = result.phase();
            }
            max = std::max(max, ExtractedMaximum{energy, i});
        }
    return max;
}

ExtractedMaximum extractorVariablePhase(int channel_count, int output_bins, const complex *const *channels, const Corrector *correctors,
                                        double *bins_buffer, ExtraData *extra_data) {
    double max_upper_bound = std::numeric_limits<double>::lowest();
    int i_max = 0;
    for (int i = 0; i < output_bins; ++i) {
        double energy_upper_bound = 0.0;
        for (int c = 0; c < channel_count; ++c) {
            energy_upper_bound += correctors[i].estimate_energy(channels[c][i]);
        }
        bins_buffer[i] = energy_upper_bound;
        if (energy_upper_bound > max_upper_bound) {
            max_upper_bound = energy_upper_bound;
            i_max = i;
        }
    }
    double energy_guess = 0;
    for (int c = 0; c < channel_count; ++c) {
        energy_guess += correctors[i_max].compute(channels[c][i_max]).energy();
    }
    bins_buffer[i_max] = energy_guess;
    ExtractedMaximum max{energy_guess, i_max};
    for (int i = 0; i < output_bins; ++i)
        if (bins_buffer[i] >= max.energy) {
            double energy = 0;
            for (int c = 0; c < channel_count; ++c) {
                CorrectorResult result = correctors[i].compute(channels[c][i]);
                double channel_energy = result.energy();
                energy += channel_energy;
                if (extra_data) {
                    extra_data[c].amplitude = result.amplitude();
                    extra_data[c].energy = channel_energy;
                    extra_data[c].phase = result.phase();
                }
            }
            max = std::max(max, ExtractedMaximum{energy, i});
        }
    return max;
}

ExtractedMaximum extractorConstantPhase(int, int, const complex *const *, const Corrector *, double *, ExtraData *) {
    ExtractedMaximum max{0, 0};
    // TODO
    return max;
}
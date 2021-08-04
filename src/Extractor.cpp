/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cassert>
#include "Extractor.h"

ExtractedMaximum extractorSingleChannel(int channel_count, int output_bins, complex *const *channels, const Corrector *correctors,
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

ExtractedMaximum extractorVariablePhase(int channel_count, int output_bins, complex *const *channels, const Corrector *correctors,
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

ExtractedMaximum extractorConstantPhase(int channel_count, int output_bins, complex *const *channels, const Corrector *correctors,
                                        double *bins_buffer, ExtraData *extra_data) {
    double max_energy = std::numeric_limits<double>::lowest();
    int i_max = 0;
    for (int i = 0; i < output_bins; ++i) {
        complex direction = 0;
        for (int c = 0; c < channel_count; ++c) {
            direction += channels[c][i] * channels[c][i];
        }
        double abs_direction = std::abs(direction);
        complex best_direction = (abs_direction > 0) ? std::sqrt(direction / abs_direction) : 1.0;
        auto corrector_result = correctors[i].compute(best_direction);
        double energy_correction = corrector_result.energy();

        double energy = 0.0;
        for (int c = 0; c < channel_count; ++c) {
            double norm_channel = std::norm(channels[c][i]);
            if (norm_channel > 0) {
                double cos_factor = (channels[c][i].real() * best_direction.real() + channels[c][i].imag() * best_direction.imag())
                                    / std::sqrt(norm_channel);
                energy += norm_channel * cos_factor * cos_factor;
            }
        }
        energy *= energy_correction;

        if (energy > max_energy) {
            max_energy = energy;
            i_max = i;
            if (extra_data) {
                const double amplitude_correction = corrector_result.amplitude();
                const double best_phase = corrector_result.phase();
                for (int c = 0; c < channel_count; ++c) {
                    double norm_channel = std::norm(channels[c][i]);
                    double abs_channel = std::sqrt(norm_channel);
                    extra_data[c].energy = norm_channel * energy_correction;
                    extra_data[c].phase = best_phase;
                    extra_data[c].amplitude = abs_channel * amplitude_correction;
                    if (abs_channel > 0) {
                        double cos_factor =
                                (channels[c][i].real() * best_direction.real() + channels[c][i].imag() * best_direction.imag()) / abs_channel;
                        extra_data[c].energy *= cos_factor * cos_factor;
                        extra_data[c].amplitude *= cos_factor;
                    }
                }
            }
        }
    }

    return {max_energy, i_max};
}

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <climits>
#include <stdexcept>
#include "BlockHelper.h"

//////////////////////////////////////////////////////////////////////////////

Block BlockHelper::create_block(PinnedArray2D<double> data, std::shared_ptr<Family> family, double scale,
                         std::shared_ptr<BlockAtomParamsConverter> converter, double booster,
                         int window_length, int output_bins, int input_shift, Extractor extractor, SpectrumCalculator& calculator, bool allow_overstep) {
    auto envelope = generate_envelope(family.get(), scale);
    auto correctors = generate_correctors(envelope, window_length, output_bins, calculator);
    return Block(std::move(data), std::move(family), scale, envelope, correctors, std::move(converter), booster, window_length, input_shift, extractor, allow_overstep);
}

PinnedArray1D<double> BlockHelper::generate_envelope(const Family* family, double scale) {
    index_t envelope_length = family->size_for_values(0.0, scale, nullptr);
    PinnedArray1D<double> envelope(envelope_length);
    family->generate_values(0.0, scale, nullptr, envelope.get(), true);
    return envelope;
}

PinnedArray1D<Corrector> BlockHelper::generate_correctors(const Array1D<double>& envelope, int window_length, int output_bins, SpectrumCalculator& calculator) {
    const index_t envelope_length = envelope.length();
    Array1D<double> envelope_squared(envelope_length);
    for (index_t i = 0; i < envelope_length; ++i) {
        envelope_squared[i] = envelope[i] * envelope[i];
    }

    PinnedArray1D<Corrector> correctors(output_bins);
    const complex *spectrum = calculator.computeSpectrum(envelope_squared, window_length);
    for (int k = 0; k < output_bins; ++k) {
        int k_twice = 2 * k;
        if (2 * k_twice > window_length) {
            correctors[k] = Corrector(std::conj(spectrum[window_length - k_twice]));
        } else {
            correctors[k] = Corrector(spectrum[k_twice]);
        }
    }
    return correctors;
}

std::map<double,int> BlockHelper::compute_scales_and_transform_sizes(const Family* family, double scale_min, double scale_max, double log_scale_step, double df_scale) {
    const double log_scale_min = std::log(scale_min);
    const double log_scale_max = std::log(scale_max);
    int il_max = (scale_min == scale_max) ? 0 :
            Types::ceil<int>((log_scale_max - log_scale_min) / log_scale_step); // TODO a może floor + 1 ?

    std::map<double,int> result;
    for (int il = 0; il <= il_max; ++il) {
        const double scale = std::exp(log_scale_min + (log_scale_max - log_scale_min) * il / il_max);
        int transform_size = round_transform_size(scale / df_scale);
        // TODO center position does not have to be 0.0 in case of sub-sample spacing
        index_t envelope_length = family->size_for_values(0.0, scale, nullptr);
        if (envelope_length > static_cast<index_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("max atom scale is too large");
        }
        result.emplace(std::make_pair(scale, std::max(transform_size, static_cast<int>(envelope_length))));
    }
    return result;
}

int BlockHelper::round_transform_size(double min_transform_size_as_float) {
    if (min_transform_size_as_float <= 0) {
        throw std::logic_error("invalid transform size passed to compute_window_length");
    }
    int transform_size_bits = std::max(0, Types::ceil<int>(log2(min_transform_size_as_float)));
    if (transform_size_bits >= sizeof(int) * CHAR_BIT - 1) {
        throw std::runtime_error("FFT transform size is too large");
    }
    return 1 << transform_size_bits;
}

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include "Block.h"

//////////////////////////////////////////////////////////////////////////////

Block::Block(PinnedArray2D<real> data_, std::shared_ptr<Family> family_, double scale, PinnedArray1D<double> envelope_,
             PinnedArray1D<Corrector> correctors_, std::shared_ptr<BlockAtomParamsConverter> converter_, double booster,
             int window_length, int input_shift, Extractor extractor, bool allow_overstep)
        : data(std::move(data_)),
          family(std::move(family_)),
          scale(scale),
          envelope(std::move(envelope_)),
          correctors(std::move(correctors_)),
          converter(std::move(converter_)),
          best_index(0) {
    assert(data.length() > 0);
    assert(data.height() > 0);
    assert(data.height() <= std::numeric_limits<int>::max());

    const int envelope_length = as_positive_int(envelope.length());
    const int output_bins = as_positive_int(correctors.length());

    assert(input_shift > 0);
    assert(window_length > 0);
    assert(window_length >= envelope_length); // the entire envelope must fit
    assert(output_bins <= window_length / 2 + 1); // cannot exceed FFT output

    index_t how_many_candidate = (data.length() - (allow_overstep ? 1 : envelope_length)) / static_cast<index_t>(input_shift) + 1;
    const int how_many = (how_many_candidate <= 0 ? 0 : as_positive_int(how_many_candidate));

    maxima = PinnedArray1D<ExtractedMaximum>(how_many);
    this->booster = Array1D<double>(output_bins);
    for (int k=0; k<output_bins; ++k) {
        const double scale_frequency = scale * std::min(k, window_length-k) / window_length;
        this->booster[k] = booster / family->optimality_factor_sf(scale_frequency);
    }

    total_request.data = data.get();
    total_request.channel_length = data.length();
    total_request.channel_count = data.height();
    total_request.input_offset = allow_overstep ? -envelope_length / 2 : 0;
    total_request.input_shift = input_shift;
    total_request.how_many = how_many;
    total_request.envelope = envelope.get();
    total_request.envelope_length = envelope_length;
    total_request.window_length = window_length;
    total_request.output_bins = output_bins;
    total_request.correctors = correctors.get();
    total_request.extractor = extractor;

    // to be filled later on
    total_request.maxima = nullptr;
    total_request.interface = nullptr;

    extended_atom_cache = std::make_shared<BlockAtomCache>();
}


SpectrogramRequest Block::buildRequest() {
    return buildRequest(0, data.length());
}

SpectrogramRequest Block::buildRequest(index_t first_sample_index, index_t end_sample_index) {
    assert(first_sample_index < end_sample_index);
    assert(end_sample_index <= total_request.channel_length);

    // removing stale elements from cache
    IndexRange updated_range{first_sample_index, end_sample_index};
    extended_atom_cache->remove_overlapping(updated_range);

    index_t samples_covered_by_first_envelope = total_request.input_offset + total_request.envelope_length;
    index_t first_index = (samples_covered_by_first_envelope > first_sample_index) ? 0 : 1 +
                                                                                         (first_sample_index - samples_covered_by_first_envelope) /
                                                                                         static_cast<index_t>(total_request.input_shift);
    index_t end_index = 1 +
                        (end_sample_index + static_cast<index_t>(total_request.envelope_length) / 2) /
                        static_cast<index_t>(total_request.input_shift);
    if (end_index > static_cast<index_t>(total_request.how_many)) {
        end_index = total_request.how_many;
    }

    SpectrogramRequest request(total_request);
    request.input_offset += first_index * static_cast<index_t>(request.input_shift);
    request.how_many = (end_index > first_index) ? static_cast<int>(end_index - first_index) : 0;
    request.maxima = maxima.get() + first_index;
    request.interface = const_cast<Block *>(this);
    return request;
}

void Block::notify() {
    // TODO not thread-safe in case of multiple requests per block
    // TODO recompute only part
    ExtractedMaximum *best = std::max_element(maxima.get(), maxima.get() + maxima.length());
    ptrdiff_t index = best - maxima.get();
    assert(index >= 0);
    assert(index < maxima.length());
    best_index = static_cast<int>(index);
}

size_t Block::get_atom_count() const {
    return static_cast<size_t>(total_request.how_many) * static_cast<size_t>(total_request.output_bins);
}

BlockAtom Block::get_best_match() const {
    return get_atom_from_index(best_index);
}

std::list<BlockAtom> Block::get_candidate_matches(double energy_to_exceed) const {
    std::list<BlockAtom> matches;
    for (int index = 0; index < total_request.how_many; ++index) {
        if (index != best_index && maxima[index].energy * booster[maxima[index].bin_index] > energy_to_exceed) {
            matches.push_back(get_atom_from_index(index));
        }
    }
    return matches;
}

BlockAtom Block::get_atom_from_index(int index) const {
    ExtractedMaximum extracted = maxima[index];
    index_t center_position =
            total_request.input_offset + static_cast<index_t>(index) * static_cast<index_t>(total_request.input_shift) + envelope.length() / 2;

    BlockAtom result = BlockAtom(
            data,
            extracted.energy,
            extracted.energy * booster[extracted.bin_index],
            family,
            static_cast<double>(extracted.bin_index) / static_cast<double>(total_request.window_length),
            static_cast<double>( center_position ),
            scale,
            total_request.extractor,
            converter
    );
    if (extended_atom_cache) {
        size_t cache_index = static_cast<size_t>(index) * static_cast<size_t>(total_request.output_bins) + static_cast<size_t>(extracted.bin_index);
        result.connect_cache(extended_atom_cache, cache_index);
    }
    return result;
}

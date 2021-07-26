/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "BlockDictionary.h"

//////////////////////////////////////////////////////////////////////////////

void BlockDictionary::add_block(double scale, int window_length, int input_shift, int output_bins,
                                Extractor extractor, SpectrumCalculator &calculator) {
    assert(window_length > 0);
    assert(input_shift > 0);
    assert(output_bins > 0);
    assert(output_bins <= window_length/2 + 1);

    index_t first_sample_index;
    index_t envelope_length = family->size_for_values(0.0, scale, &first_sample_index);
    assert(envelope_length <= static_cast<index_t>(window_length));
    PinnedArray1D<double> envelope(envelope_length);
    Array1D<double> envelope_squared(envelope_length);
    family->generate_values(0.0, scale, nullptr, envelope.get(), true);
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

    /*
    printf(
            "<block uses=\"GAUSS-WINDOW\"><param name=\"type\" value=\"gabor\"/><param name=\"windowLen\" value=\"%d\"/><param name=\"windowShift\" value=\"%d\"/><param name=\"fftSize\" value=\"%d\"/></block>\n",
            (int) envelope_length, input_shift, window_length
    );
    */
    Block block(data, family, std::move(envelope), std::move(correctors), window_length, input_shift, extractor, allow_overstep);
    if (block.get_how_many() > 0) {
        blocks.emplace(scale, std::move(block));
    }
}

size_t BlockDictionary::get_atom_count() {
    size_t result = 0;
    for (auto &kv : blocks) {
        result += kv.second.get_atom_count();
    }
    return result;
}

BasicAtomPointer BlockDictionary::get_best_match() {
    std::optional<BlockAtom> result;
    for (auto &kv : blocks) {
        BlockAtom atom = kv.second.get_best_match();
        atom.scale = kv.first;
        if (!result || *result < atom) {
            result = atom;
        }
    }
    if (!result) {
        return nullptr;
    }
    return std::make_shared<BlockAtom>(*std::move(result));
}

int BlockDictionary::get_max_window_length() const {
    int max_window_length = 0;
    for (auto &kv: blocks) {
        max_window_length = std::max(max_window_length, kv.second.get_window_length());
    }
    return max_window_length;
}

void BlockDictionary::fetch_proto_requests(std::list<ProtoRequest> &requests) {
    for (auto &kv: blocks) {
        requests.push_back(kv.second.buildRequest(0, data.length()));
    }
}

void BlockDictionary::fetch_requests(IndexRange signal_range, std::list<SpectrogramRequest> &requests) {
    for (auto &kv: blocks) {
        SpectrogramRequest request = kv.second.buildRequest(signal_range.first_index, signal_range.end_index);
        if (request.how_many > 0) {
            requests.push_back(std::move(request));
        }
    }
}

IndexRange BlockDictionary::subtract_from_signal(const ExtendedAtom &atom) {
    auto block_atom = dynamic_cast<const BlockExtendedAtom &>(atom);
    index_t first_sample_offset;
    const index_t sample_count = family->size_for_values(block_atom.position, block_atom.scale, &first_sample_offset);
    index_t end_sample_offset = first_sample_offset + sample_count;
    if (end_sample_offset <= 0) {
        return {0, 0};
    }

    const index_t first_valid_sample_offset = std::max<index_t>(0, first_sample_offset);
    const index_t end_valid_sample_offset = std::min<index_t>(data.length(), end_sample_offset);

    Array1D<double> samples(sample_count);
    family->generate_values(block_atom.position, block_atom.scale, nullptr, samples.get(), false);
    double common_amplitude_factor = 1.0 / family->value(0.0);

    const int channel_count = data.height();
    const double omega = 2 * M_PI * block_atom.frequency;
    for (int c = 0; c < channel_count; ++c) {
        const double channel_amplitude_factor = block_atom.extra[c].amplitude * common_amplitude_factor;
        for (index_t i = first_valid_sample_offset; i < end_valid_sample_offset; ++i) {
            data[c][i] -= channel_amplitude_factor * samples[i - first_sample_offset] *
                          std::cos(omega * (static_cast<double>(i) - block_atom.position) + block_atom.extra[c].phase);
        }
    }

    return {first_valid_sample_offset, end_valid_sample_offset};
}

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "BlockAtomObjective.h"

//////////////////////////////////////////////////////////////////////////////

BlockAtomObjective::BlockAtomObjective(std::shared_ptr<Family> family, Array2D<double> data_, Extractor extractor,
                                       std::shared_ptr<BlockAtomParamsConverter> converter_) :
        family(std::move(family)),
        channel_count(data_.height()),
        data(std::move(data_)),
        products(channel_count, 1),
        extractor(extractor),
        converter(std::move(converter_)) {}

double BlockAtomObjective::calculate_energy(const std::array<double,3> &array, double *out_norm, ExtraData *out_extra_data) {
    auto params_pair = converter->paramsFromArray(array);
    BlockAtomParams& params = params_pair.first;
    double penalty = params_pair.second;

    index_t envelope_length = family->size_for_values(params.position, params.scale, nullptr);
    index_t envelope_offset;
    envelope.resize(envelope_length);
    oscillating.resize(envelope_length);
    double norm = family->generate_values(params.position, params.scale, &envelope_offset, envelope.data(), true);
    if (out_norm) {
        *out_norm = norm;
    }

    const double omega = 2 * M_PI * params.frequency;
    for (index_t i=0; i < envelope_length; ++i) {
        const double t = static_cast<double>(envelope_offset + i) - params.position;
        oscillating[i] = std::polar(envelope[i], -omega * t);
    }

    complex FT = 0.0;
    for (index_t i = 0; i < envelope_length; ++i) {
        FT += oscillating[i] * oscillating[i];
    }
    index_t first_sample_offset = std::max<index_t>(0, envelope_offset);
    index_t last_sample_offset = std::max<index_t>(0, envelope_offset + envelope_length - 1);

    first_sample_offset = std::min(data.length() - 1, first_sample_offset);
    last_sample_offset = std::min(data.length() - 1, last_sample_offset);

    products.fill(0.0);
    for (index_t i = first_sample_offset; i <= last_sample_offset; ++i) {
        complex z = oscillating[i - envelope_offset];
        for (int c = 0; c < channel_count; ++c) {
            products[c][0] += data[c][i] * z;
        }
    }

    Corrector corrector(FT);
    double tmp_for_extractor;
    double result = std::exp(-penalty) * extractor(channel_count, 1, products.get(), &corrector, &tmp_for_extractor, out_extra_data).energy;
    return result;
}

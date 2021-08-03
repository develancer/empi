/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cmath>
#include <memory>
#include <vector>
#include "Array.h"
#include "Atom.h"
#include "BlockAtom.h"
#include "Corrector.h"

//////////////////////////////////////////////////////////////////////////////

BlockAtom::BlockAtom(Array2D<double> data, double energy, std::shared_ptr<Family> family_, double frequency, double position, Extractor extractor)
        : BasicAtom(std::move(data), energy), BlockAtomTrait(frequency, position, NAN), extractor(extractor), family(std::move(family_)) {}

ExtendedAtomPointer BlockAtom::extend() const {
    index_t envelope_length;
    index_t envelope_offset;
    std::vector<double> envelope;

    // TODO local parameter optimization

    envelope_length = family->size_for_values(position, scale, nullptr);
    envelope.resize(envelope_length);
    double norm = family->generate_values(position, scale, &envelope_offset, envelope.data(), true);

    complex FT = 0.0;
    const double arg_for_integrals = 4 * M_PI * frequency;
    for (index_t i = 0; i < envelope_length; ++i) {
        const double t = static_cast<double>(envelope_offset + i) - position;
        const double value2 = envelope[i] * envelope[i];
        FT += value2 * std::polar(1.0, -arg_for_integrals * t);
    }
    index_t first_sample_offset = std::max<index_t>(0, envelope_offset);
    index_t last_sample_offset = std::max<index_t>(0, envelope_offset + envelope_length - 1);

    first_sample_offset = std::min(data().length() - 1, first_sample_offset);
    last_sample_offset = std::min(data().length() - 1, last_sample_offset);

    const int channel_count = data().height();
    Array2D<complex> products(channel_count, 1);
    Array1D<ExtraData> extra_data(channel_count);
    products.fill(0.0);

    for (index_t i = first_sample_offset; i <= last_sample_offset; ++i) {
        double t = static_cast<double>(i) - position;
        for (int c = 0; c < channel_count; ++c) {
            products[c][0] += data()[c][i] * envelope[i - envelope_offset] * std::polar(1.0, -2 * M_PI * frequency * t);
        }
    }

    Corrector corrector(FT);
    const double amplitude_factor = norm * family->value(0.0);
    double tmp_for_extractor;
    extractor(channel_count, 1, products.get(), &corrector, &tmp_for_extractor, extra_data.get());
    for (int i = 0; i < channel_count; ++i) {
        extra_data[i].amplitude *= amplitude_factor;
    }

    return std::static_pointer_cast<ExtendedAtom>(
            std::make_shared<BlockExtendedAtom>(
                    data(), get_energy(),
                    frequency, position, scale,
                    std::move(extra_data)
            )
    );
}

//////////////////////////////////////////////////////////////////////////////

BlockExtendedAtom::BlockExtendedAtom(Array2D<double> data, double energy, double frequency, double position, double scale,
                                     Array1D<ExtraData> extra)
        : ExtendedAtom(std::move(data), energy), BlockAtomTrait(frequency, position, scale), extra(std::move(extra)) {}

void BlockExtendedAtom::export_atom(std::list<ExportedAtom> *atoms) {
    for (index_t c = 0; c < extra.length(); ++c) {
        const ExtraData &e = extra[c];
        ExportedAtom atom(e.energy);
        atom.amplitude = e.amplitude;
        atom.envelope = "gauss"; // TODO
        atom.frequency = frequency;
        atom.scale = scale;
        atom.phase = e.phase;
        atom.position = position;
        atoms[c].push_back(std::move(atom));
    }
}

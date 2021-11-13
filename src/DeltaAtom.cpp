/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "DeltaAtom.h"

//////////////////////////////////////////////////////////////////////////////

DeltaAtom::DeltaAtom(Array2D<double> data, double energy, index_t position)
        : BasicAtom(std::move(data), energy), position(position) {}

ExtendedAtomPointer DeltaAtom::extend(bool allow_optimization) {
    const int channel_count = data.height();
    std::vector<double> amplitudes(channel_count);
    for (int c = 0; c < channel_count; ++c) {
        amplitudes[c] = data[c][position];
    }
    return std::make_shared<DeltaExtendedAtom>(data, energy, position, std::move(amplitudes));
}

[[nodiscard]] double DeltaAtom::get_energy_upper_bound() const {
    return energy;
}

//////////////////////////////////////////////////////////////////////////////

DeltaExtendedAtom::DeltaExtendedAtom(Array2D<double> data, double energy, index_t position, std::vector<double> &&amplitudes)
        : ExtendedAtom(std::move(data), energy), position(position), amplitudes(amplitudes) {}

void DeltaExtendedAtom::export_atom(std::list<ExportedAtom> *atoms) {
    const int channel_count = data.height();
    for (int c = 0; c < channel_count; ++c) {
        ExportedAtom atom(amplitudes[c] * amplitudes[c]);
        atom.amplitude = amplitudes[c];
        atom.envelope = "delta";
        atom.position = static_cast<double>(position);
        atoms[c].push_back(std::move(atom));
    }
}

IndexRange DeltaExtendedAtom::subtract_from_signal() const {
    const int channel_count = data.height();
    for (int c = 0; c < channel_count; ++c) {
        data[c][position] -= amplitudes[c];
    }
    return {position, position + 1};
}

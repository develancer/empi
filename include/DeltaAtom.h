/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_DELTA_ATOM_H
#define EMPI_DELTA_ATOM_H

#include <vector>
#include "Atom.h"

/**
 * Basic atom implementation for delta (single sample) dictionaries.
 */
class DeltaAtom : public BasicAtom {
    index_t position;

public:
    /**
     * Create a new basic atom for delta (single sample) dictionary.
     *
     * @param position position of the sample
     * @param amplitudes amplitudes for the sample for all channels
     */
    DeltaAtom(Array2D<double> data, double energy, index_t position);

    /**
     * Create a new BlockExtendedAtom instance as a fully-defined representation of this atom.
     *
     * @return smart pointer to a newly created BlockExtendedAtom
     */
    [[nodiscard]] ExtendedAtomPointer extend(bool allow_optimization) final;

    /**
     * @return estimate for maximum possible energy of the atom that can be obtained by locally optimizing the coefficients
     */
    [[nodiscard]] double get_energy_upper_bound() const final;
};

/**
 * Extended atom implementation for delta (single sample) dictionaries.
 * It defines position and amplitudes for all channels.
 * Some of the amplitudes may be negative.
 */
class DeltaExtendedAtom : public ExtendedAtom {
    index_t position;
    std::vector<double> amplitudes;

public:
    /**
     * Create a new basic atom for delta (single sample) dictionary.
     *
     * @param position position of the sample
     * @param amplitudes amplitudes for the sample for all channels
     */
    DeltaExtendedAtom(Array2D<double> data, double energy, index_t position, std::vector<double>&& amplitudes);

    /**
     * Export all parameters of this atom to a uniform representation.
     *
     * @param atoms array of lists of atoms (one list per channel) to which atoms' parameters should be appended
     */
    void export_atom(std::list<ExportedAtom> atoms[]) final;

    /**
     * Subtract this atom from the multi-channel signal being analyzed.
     *
     * @return range of the signal samples that are being changed
     * (can be later passed to fetch_requests in all dictionaries)
     */
    IndexRange subtract_from_signal() const final;
};

#endif //EMPI_DELTA_ATOM_H

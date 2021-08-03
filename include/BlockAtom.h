/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_ATOM_H
#define EMPI_BLOCK_ATOM_H

#include <cmath>
#include <list>
#include <memory>
#include "Array.h"
#include "Atom.h"
#include "Extractor.h"
#include "Family.h"

/**
 * Simple structure with three fields, used in both BlockAtom and ExtendedBlockAtom.
 */
struct BlockAtomTrait {
    double frequency; // between 0 and 0.5 (1 would be sampling frequency)
    double position; // in samples
    double scale; // in samples

    BlockAtomTrait(double frequency, double position, double scale)
            : frequency(frequency), position(position), scale(scale) {}
};

/**
 * Basic atom implementation for block dictionaries.
 * It defines position, frequency and scale, but neither amplitude nor phase.
 */
class BlockAtom : public BasicAtom, public BlockAtomTrait {
    Extractor extractor;
    std::shared_ptr<Family> family;

public:
    /**
     * Create a new basic atom for block dictionary.
     *
     * @param data reference to multi-channel data of the analysed signal
     * @param energy energy of the atom, computed as sum of all samples squared
     * @param family Family object describing properties of the envelope function
     * @param frequency frequency in range between 0 and 0.5 (1 would be sampling frequency)
     * @param position center position in samples
     * @param extractor  function that will be used to extract information from multi-channel results,
     * according to a particular mode of multi-channel operation (constant vs variable phase)
     */
    BlockAtom(Array2D<double> data, double energy, std::shared_ptr<Family> family, double frequency, double position, Extractor extractor);

    /**
     * Create a new BlockExtendedAtom instance as a fully-defined representation of this atom.
     *
     * @return smart pointer to a newly created BlockExtendedAtom
     */
    ExtendedAtomPointer extend() const final;
};

/**
 * Extended atom implementation for block dictionaries.
 * It defines all parameters, including amplitudes, phases and energies for all channels.
 */
class BlockExtendedAtom : public ExtendedAtom, public BlockAtomTrait {
public:
    /**
     * Channel-specific data (amplitudes, phases and energies) for all channels of the analysed signal.
     */
    Array1D<ExtraData> extra;

    /**
     * Create a new extended atom for block dictionary.
     *
     * @param data reference to multi-channel data of the analysed signal
     * @param energy energy of the atom, computed as sum of all samples squared
     * @param frequency frequency in range between 0 and 0.5 (1 would be sampling frequency)
     * @param position center position in samples
     * @param scale atom's scale (related to its duration) in samples
     * @param extra channel-specific data (amplitudes, phases and energies) for all channels of the analysed signal
     */
    BlockExtendedAtom(Array2D<double> data, double energy, double frequency, double position, double scale,
                      Array1D<ExtraData> extra);

    /**
     * Export all parameters of this atom to a uniform representation.
     *
     * @param atoms array of lists of atoms (one list per channel) to which atoms' parameters should be appended
     */
    void export_atom(std::list<ExportedAtom> atoms[]) final;
};

#endif //EMPI_BLOCK_ATOM_H

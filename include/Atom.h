/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_ATOM_H
#define EMPI_ATOM_H

#include <list>
#include <memory>
#include "Array.h"
#include "ExportedAtom.h"
#include "IndexRange.h"

/**
 * Base class for all atom representations.
 * Each atom stores a set of parameters in a dictionary parameter space,
 * has a well-defined energy and can be used to generate a waveform.
 */
class Atom {
protected:
    explicit Atom(Array2D<double> data, double energy) : data(std::move(data)), energy(energy) {}

public:
    /**
     * reference to multi-channel data of the analysed signal
     */
    const Array2D<double> data;

    /**
     * energy of the atom, computed as sum of all samples squared
     */
    const double energy;

    /**
     * Compares atoms according to their energies.
     * @return true if atom has smaller energy than the other atom, false otherwise
     */
    bool operator<(const Atom &other) const {
        return energy < other.energy;
    }

    virtual ~Atom() = default;
};

/**
 * Base class for all atom representations for which all parameters are already known.
 * Instances of this class can be obtained from BasicAtoms' method extend().
 */
class ExtendedAtom : public Atom {
protected:
    explicit ExtendedAtom(Array2D<double> data, double energy) : Atom(std::move(data), energy) {}

public:
    /**
     * Export all parameters of this atom to a uniform representation.
     *
     * @param atoms array of lists of atoms (one list per channel) to which atoms' parameters should be appended
     */
    virtual void export_atom(std::list<ExportedAtom> atoms[]) = 0;

    /**
     * Subtract this atom from the multi-channel signal being analyzed.
     *
     * @return range of the signal samples that are being changed
     * (can be later passed to fetch_requests in all dictionaries)
     */
    virtual IndexRange subtract_from_signal() const = 0;
};

using ExtendedAtomPointer = std::shared_ptr<ExtendedAtom>;

/**
 * Base class for all atom representations for which some of the parameters may not be known yet.
 * Instance of ExtendedAtom can be obtained by calling the method extend().
 */
class BasicAtom : public Atom {
public:
    explicit BasicAtom(Array2D<double> data, double energy) : Atom(std::move(data), energy) {}

    /**
     * Create a new ExtendedAtom instance as a fully-defined representation of this atom.
     *
     * @return smart pointer to a newly created ExtendedAtom
     */
    [[nodiscard]] virtual ExtendedAtomPointer extend(bool allow_optimization) = 0;

    /**
     * @return estimate for maximum possible energy of the atom that can be obtained by locally optimizing the coefficients
     */
    [[nodiscard]] virtual double get_energy_upper_bound() const = 0;
};

using BasicAtomPointer = std::shared_ptr<BasicAtom>;

#endif //EMPI_ATOM_H

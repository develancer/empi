/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_ATOM_H
#define EMPI_ATOM_H

#include <list>
#include <memory>
#include "Array.h"
#include "ExportedAtom.h"

/**
 * Base class for all atom representations.
 * Each atom stores a set of parameters in a dictionary parameter space,
 * has a well-defined energy and can be used to generate a waveform.
 */
class Atom {
    Array2D<double> data_;
    double energy;

protected:
    explicit Atom(Array2D<double> data, double energy) : data_(std::move(data)), energy(energy) {}

public:
    /**
     * Compares atoms according to their energies.
     * @return true if atom has smaller energy than the other atom, false otherwise
     */
    bool operator<(const Atom &other) const {
        return energy < other.energy;
    }

    /**
     * @return reference to multi-channel data of the analysed signal
     */
    [[nodiscard]] const Array2D<double> &data() const {
        return data_;
    }

    /**
     * @return energy of the atom, computed as sum of all samples squared
     */
    [[nodiscard]] double get_energy() const {
        return energy;
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
    virtual ExtendedAtomPointer extend() = 0;
};

using BasicAtomPointer = std::shared_ptr<BasicAtom>;

#endif //EMPI_ATOM_H

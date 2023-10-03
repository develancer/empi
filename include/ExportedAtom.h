/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_EXPORTED_ATOM_H
#define EMPI_EXPORTED_ATOM_H

#include <string>

/**
 * Uniform plain data representation for all different types of dictionary atoms.
 */
class ExportedAtom {
public:
    /** atom's envelope type (e.g. "gauss") */
    std::string envelope;

    /** atom's amplitude */
    double amplitude;

    /** atom's energy, computed as sum of all samples squared */
    double energy;

    /** atom's frequency between 0 and 0.5 (1 would be sampling frequency) */
    double frequency;

    /** atom's phase between ±π */
    double phase;

    /** atom's scale in samples */
    double scale;

    /** atom's central position in samples */
    double position;

    /**
     * Create a new ExportedAtom with a given energy.
     * All other fields will be initialized to NAN.
     *
     * @param energy atom's energy, computed as sum of all samples squared
     */
    explicit ExportedAtom(double energy);
};

#endif //EMPI_EXPORTED_ATOM_H

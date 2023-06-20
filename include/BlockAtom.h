/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_ATOM_H
#define EMPI_BLOCK_ATOM_H

#include <atomic>
#include <cmath>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include "Array.h"
#include "Atom.h"
#include "BlockAtomBase.h"
#include "BlockAtomCache.h"
#include "Extractor.h"
#include "Family.h"

/**
 * Extended atom implementation for block dictionaries.
 * It defines all parameters, including amplitudes, phases and energies for all channels.
 */
class BlockExtendedAtom : public ExtendedAtom {
    std::shared_ptr<Family> family;

public:
    const BlockAtomParams params;

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
    BlockExtendedAtom(Array2D<double> data, double energy, std::shared_ptr<Family> family, double frequency, double position, double scale,
                      Array1D<ExtraData> extra);

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

/**
 * Basic atom implementation for block dictionaries.
 * It defines position, frequency and scale, but neither amplitude nor phase.
 */
class BlockAtom : public BasicAtom {
    double energy_upper_bound;
    Extractor extractor;
    std::shared_ptr<Family> family;
    std::shared_ptr<BlockAtomParamsConverter> converter;

    std::optional<BlockAtomCacheSlot> cache_slot;

    static std::atomic<size_t> failed_optimization_count;
    static std::atomic<size_t> total_optimization_count;

public:
    const BlockAtomParams params;

    static double get_failed_optimization_percent()
    {
        return 100.0 * failed_optimization_count / total_optimization_count;
    }

    /**
     * Create a new basic atom for block dictionary.
     *
     * @param data reference to multi-channel data of the analysed signal
     * @param energy energy of the atom, computed as sum of all samples squared
     * @param family Family object describing properties of the envelope function
     * @param frequency frequency in range between 0 and 0.5 (1 would be sampling frequency)
     * @param position center position in samples
     * @param scale in samples
     * @param extractor  function that will be used to extract information from multi-channel results,
     * according to a particular mode of multi-channel operation (constant vs variable phase)
     */
    BlockAtom(Array2D<double> data, double energy, double energy_upper_bound, std::shared_ptr<Family> family,
              double frequency, double position, double scale,
              Extractor extractor, std::shared_ptr<BlockAtomParamsConverter> converter);

    void connect_cache(std::shared_ptr<BlockAtomCache> cache, size_t key);

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

#endif //EMPI_BLOCK_ATOM_H

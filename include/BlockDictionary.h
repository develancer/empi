/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_DICTIONARY_H
#define EMPI_BLOCK_DICTIONARY_H

#include <list>
#include <memory>
#include <set>

#include "Atom.h"
#include "Block.h"
#include "BlockDictionaryStructure.h"
#include "Family.h"
#include "Dictionary.h"
#include "IndexRange.h"
#include "PinnedArray.h"
#include "SpectrogramRequest.h"
#include "SpectrumCalculator.h"

/**
 * Dictionary for envelope-based oscillating atoms, consisting of one or more Block instances.
 */
class BlockDictionary : public Dictionary {

    std::list<Block> blocks;

public:
    /**
     * Create a new dictionary for a given envelope function.
     *
     * @param data reference to multi-channel data of the analysed signal
     * @param family Family object describing properties of the envelope function
     * @param scale atom's scale (for scaling an envelope function's argument) in samples
     * @param extractor function that will be used to extract information from multi-channel results,
     * according to a particular mode of multi-channel operation (constant vs variable phase)
     * @param calculator needed for calculating corrections to account for normalization issues
     */
    BlockDictionary(const BlockDictionaryStructure& structure, const PinnedArray2D<double>& data,
                    Extractor extractor, SpectrumCalculator &calculator, bool allow_overstep);

    explicit BlockDictionary(Block&& block);

    /**
     * Calculate how many atoms are represented by this dictionary.
     *
     * @return number of atoms in this dictionary
     */
    size_t get_atom_count() final;

    /**
     * Compute and return the atom that is currently a best match for the analyzed signal.
     *
     * @return best matching atom
     */
    BasicAtomPointer get_best_match() final;

    std::list<BasicAtomPointer> get_candidate_matches(double energy_to_exceed) final;

    /**
     * Create all request templates that can be requested by this dictionary.
     * The generated proto-requests will be added to the given list.
     *
     * @param requests list for the requests to be appended to
     */
    void fetch_proto_requests(std::list<ProtoRequest> &requests) final;

    /**
     * Create all actual recalculation requests for given updated signal range.
     * The generated requests will be added to the given list.
     *
     * @param signal_range range of the signal samples that have been changed
     * (e.g. as returned from the last call to subtract_from_signal)
     * @param requests list for the requests to be appended to
     */
    void fetch_requests(IndexRange signal_range, std::list<SpectrogramRequest> &requests) final;
};

#endif //EMPI_BLOCK_DICTIONARY_H

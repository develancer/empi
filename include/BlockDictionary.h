/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_DICTIONARY_H
#define EMPI_BLOCK_DICTIONARY_H

#include <list>
#include <map>
#include <memory>

#include "Atom.h"
#include "Block.h"
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

    PinnedArray2D<double> data;
    std::shared_ptr<Family> family;
    std::map<double, Block> blocks;
    const bool allow_overstep;

public:
    /**
     * Create a new dictionary for a given envelope function.
     *
     * @param data reference to multi-channel data of the analysed signal
     * @param family Family object describing properties of the envelope function
     * @param allow_overstep whether we can assume that samples before and after the actual signal range are equal to zero
     */
    BlockDictionary(PinnedArray2D<double> data, std::shared_ptr<Family> family, bool allow_overstep = true)
            : data(std::move(data)), family(std::move(family)), allow_overstep(allow_overstep) {}

    /**
     * Add a block for a given scale to the dictionary.
     *
     * @param scale atom's scale (for scaling an envelope function's argument) in samples
     * @param window_length number of samples of the FFT to which data will be zero-padded, should be a power of 2
     * @param input_shift shift (in samples) between consecutive time-shifted envelope realizations
     * @param output_bins number of meaningful frequency bins from the FFT, should not exceed window_length/2+1
     * @param extractor function that will be used to extract information from multi-channel results,
     * according to a particular mode of multi-channel operation (constant vs variable phase)
     * @param calculator needed for calculating corrections to account for normalization issues
     */
    void add_block(double scale, int window_length, int input_shift, int output_bins, Extractor extractor, SpectrumCalculator &calculator);

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

    /**
     * @return largest number of samples (including zero-padding) of FFTs from all blocks in the dictionary
     */
    [[nodiscard]] int get_max_window_length() const;

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

    /**
     * Subtract given atom from the multi-channel signal being analyzed.
     *
     * @param atom atom to be subtracted from signal
     * @return range of the signal samples that are being changed
     * (can be later passed to fetch_requests)
     */
    IndexRange subtract_from_signal(const ExtendedAtom &atom) final;
};

#endif //EMPI_BLOCK_DICTIONARY_H

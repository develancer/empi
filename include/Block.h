/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_H
#define EMPI_BLOCK_H

#include <list>
#include <map>
#include "BlockAtom.h"
#include "BlockAtomBase.h"
#include "BlockAtomCache.h"
#include "BlockInterface.h"
#include "Extractor.h"
#include "PinnedArray.h"
#include "SpectrogramRequest.h"
#include "Types.h"

/**
 * Basic building block for all dictionaries.
 * Each block includes a well-defined envelope waveform of a given width,
 * and is responsible for computing scalar products between signal and all atoms
 * that can be obtained by time-shifting and modulating that waveform.
 */
class Block : public BlockInterface {

    PinnedArray2D<real> data;
    std::shared_ptr<Family> family;
    double scale;
    PinnedArray1D<double> envelope;
    PinnedArray1D<Corrector> correctors;
    std::shared_ptr<BlockAtomParamsConverter> converter;

    PinnedArray1D<ExtractedMaximum> maxima;
    Array1D<double> booster;
    SpectrogramRequest total_request;

    int best_index;

    std::shared_ptr<BlockAtomCache> extended_atom_cache;

public:
    /**
     * Create a new block for the dictionary.
     *
     * @param data reference to multi-channel data of the analysed signal
     * @param family Family object describing properties of the envelope function
     * @param envelope L²-normalized samples of the particular realization of the envelope function
     * @param correctors array of corrections to be applied on the computed spectrum to account for normalization issues
     * @param steps TODO
     * @param window_length number of samples of the FFT to which data will be zero-padded, should be a power of 2
     * @param input_shift shift (in samples) between consecutive time-shifted envelope realizations
     * @param extractor function that will be used to extract information from multi-channel results,
     * according to a particular mode of multi-channel operation (constant vs variable phase)
     * @param allow_overstep whether we can assume that samples before and after the actual signal range are equal to zero
     */
    Block(PinnedArray2D<double> data, std::shared_ptr<Family> family, double scale, PinnedArray1D<double> envelope,
          PinnedArray1D<Corrector> correctors, std::shared_ptr<BlockAtomParamsConverter> converter,
          int window_length, int input_shift, Extractor extractor, bool allow_overstep = true);

    /**
     * Create a SpectrogramRequest that can be used to recompute the entire block.
     */
    SpectrogramRequest buildRequest();

    /**
     * Create a SpectrogramRequest that can be used to recompute part of this block.
     *
     * @param first_sample_index zero-based index of the first sample in the interval that needs to be updated
     * @param end_sample_index zero-based index of the first sample _following_ the interval that needs to be updated (last sample index + 1)
     * @return
     */
    SpectrogramRequest buildRequest(index_t first_sample_index, index_t end_sample_index);

    /**
     * This method is called internally by Computer instances
     * to notify that the request obtained by the last call to buildRequest is completed.
     * It is then the block's responsibility to update its internal cache.
     */
    void notify() final;

    /**
     * Calculate how many atoms are represented by this block.
     *
     * @return number of atoms in this block
     */
    [[nodiscard]] size_t get_atom_count() const;

    /**
     * Compute and return the atom that is currently a best match for the analyzed signal.
     *
     * @return best matching atom
     */
    [[nodiscard]] BlockAtom get_best_match() const;

    [[nodiscard]] std::list<BlockAtom> get_candidate_matches(double energy_to_exceed) const;

private:
    [[nodiscard]] BlockAtom get_atom_from_index(int index) const;
};

#endif //EMPI_BLOCK_H

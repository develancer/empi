/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_DELTA_DICTIONARY_H
#define EMPI_DELTA_DICTIONARY_H

#include <list>
#include "Dictionary.h"
#include "PinnedArray.h"

/**
 * Base class for dictionaries of all types.
 * One or more instances of any class derived from Dictionary
 * can be associated with a single Computer instance.
 */
class DeltaDictionary : public Dictionary {
    PinnedArray2D<double> data;

public:
    explicit DeltaDictionary(PinnedArray2D<double> data);

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

    std::list<BasicAtomPointer> get_candidate_matches(double energy_to_exceed) final {
        return {};
    }

    /**
     * Create all request templates that can be requested by this dictionary.
     * The generated proto-requests will be added to the given list.
     *
     * @param requests list for the requests to be appended to
     */
    void fetch_proto_requests(std::list<ProtoRequest> &requests) final {};

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

#endif //EMPI_DELTA_DICTIONARY_H

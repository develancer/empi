/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "DeltaAtom.h"
#include "DeltaDictionary.h"

//////////////////////////////////////////////////////////////////////////////

DeltaDictionary::DeltaDictionary(PinnedArray2D<double> data) : data(std::move(data))
{ }

size_t DeltaDictionary::get_atom_count() {
    return data.length();
}

BasicAtomPointer DeltaDictionary::get_best_match() {
    const int channel_count = data.height();
    const index_t sample_count = data.length();

    std::pair<index_t, double> best_match;
    for (index_t i=0; i<sample_count; ++i) {
        double energy = 0.0;
        for (int c=0; c<channel_count; ++c) {
            energy += data[c][i] * data[c][i];
        }
        if (energy > best_match.second) {
            best_match = {i, energy};
        }
    }
    return std::make_shared<DeltaAtom>(data, best_match.second, best_match.first);
}

void DeltaDictionary::fetch_requests(IndexRange signal_range, std::list<SpectrogramRequest> &requests) {
    // TODO recalculate only part
}

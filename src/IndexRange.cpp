/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include "IndexRange.h"

//////////////////////////////////////////////////////////////////////////////

IndexRange IndexRange::overlap(const IndexRange &other) const {
    index_t new_first_index = std::max(first_index, other.first_index);
    index_t new_end_index = std::min(end_index, other.end_index);
    if (new_first_index < new_end_index) {
        return {new_first_index, new_end_index};
    } else {
        return {};
    }
}

bool IndexRange::includes(index_t index) const {
    return first_index <= index && index < end_index;
}

bool IndexRange::operator!() const {
    return first_index >= end_index;
}

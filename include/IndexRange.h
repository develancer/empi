/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_INDEX_RANGE_H
#define EMPI_INDEX_RANGE_H

#include "Types.h"

/**
 * Simple data object representing a range of samples in a signal.
 * Represented range is from first_index (inclusive) to end_index (exclusive),
 * which means the last sample in the range has an index of end_index-1.
 */
struct IndexRange {
    /**
     * Index of the first sample in the range.
     */
    index_t first_index;

    /**
     * Index of the sample following the last sample in the range.
     */
    index_t end_index;

    /**
     * Create a new range starting at the beginning of the signal, with length (in samples) equal to end_index.
     *
     * @param end_index length of the range = index of the sample following the last sample in the range
     */
    IndexRange(index_t end_index = 0) : first_index(0), end_index(end_index) {} // NOLINT

    /**
     * Create a new range.
     *
     * @param first_index index of the first sample in the range
     * @param end_index index of the sample following the last sample in the range
     */
    IndexRange(index_t first_index, index_t end_index) : first_index(first_index), end_index(end_index) {}

    /**
     * Calculate the intersection between two ranges.
     *
     * @param other the other range
     * @return intersection of the two ranges, or IndexRange(0) if the ranges don't overlap
     */
    [[nodiscard]] IndexRange overlap(const IndexRange &other) const;

    [[nodiscard]] bool includes(index_t index) const;

    /**
     * @return true if range is empty, false otherwise
     */
    bool operator!() const;
};

#endif //EMPI_INDEX_RANGE_H

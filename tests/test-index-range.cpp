/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "IndexRange.h"
#include "Testing.h"

void assertNoOverlap(index_t start, index_t end, index_t other_length) {
    auto overlap = IndexRange(start, end).overlap(other_length);
    ASSERT(!overlap);
}

void assertOverlap(index_t start, index_t end, index_t other_length, index_t first_index, index_t end_index) {
    auto overlap = IndexRange(start, end).overlap(other_length);
    ASSERT_EQUALS(first_index, overlap.first_index);
    ASSERT_EQUALS(end_index, overlap.end_index);
}

int main() {
    assertOverlap(10, 30, 50, 10, 30);
    assertOverlap(10, 50, 50, 10, 50);
    assertOverlap(10, 60, 50, 10, 50);
    assertOverlap(-10, 10, 30, 0, 10);
    assertOverlap(-10, 40, 30, 0, 30);
    assertNoOverlap(-10, 0, 30);
    assertNoOverlap(30, 40, 30);
    puts("OK");
}

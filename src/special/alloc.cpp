/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdlib>
#include "PinnedArray.h"

// This file is meant as a replacement for CUDA-based pinned memory allocation
// in case it isn't needed (for some unit tests) or CUDA is not available.
// Each executable should link with either alloc.cpp or alloc.cu, but not both at once.

void *cuda_host_malloc(size_t length) {
    void *result = malloc(length);
    if (!result) {
        throw std::bad_alloc();
    }
    return result;
}

void cuda_host_free(void *pointer) {
    free(pointer);
}

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "PinnedArray.h"

// This file provides allocation of pinned (page-locked) host memory
// so it can be safely used in some CUDA routines.
// Each executable should link with either alloc.cpp or alloc.cu, but not both at once.

static bool CUDA_ALLOC_DISABLED = false;

void cuda_host_disable() {
    CUDA_ALLOC_DISABLED = true;
}

void *cuda_host_malloc(size_t length) {
    void *result;
    if (CUDA_ALLOC_DISABLED) {
        result = malloc(length);
        if (!result) {
            throw std::bad_alloc();
        }
    } else {
        if (cudaMallocHost(&result, length)) {
            throw std::bad_alloc();
        }
    }
    return result;
}

void cuda_host_free(void *pointer) {
    if (CUDA_ALLOC_DISABLED) {
        free(pointer);
    } else {
        cudaFreeHost(pointer);
    }
}

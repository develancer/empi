/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_PINNED_ARRAY_H
#define EMPI_PINNED_ARRAY_H

#include "Array.h"

/**
 * Allocate pinned memory area of a given size.
 * This defaults to ordinary malloc in case of CUDA-less compilation.
 *
 * @param length requested size of the memory area
 * @return pointer to memory area
 * @throw std::bad_alloc if allocation fails
 */
void *cuda_host_malloc(size_t length);

/**
 * Free given area of pinned memory, previously allocated with cuda_host_malloc.
 * This defaults to ordinary free in case of CUDA-less compilation.
 *
 * @param pointer pointer to memory area
 */
void cuda_host_free(void *pointer);

template<typename T>
T *cuda_host_alloc(size_t length) {
    return reinterpret_cast<T *>(cuda_host_malloc(sizeof(T) * length));
}

/**
 * Subclass of Array1D template for allocating one-dimensional arrays
 * using page-locked (pinned) memory which can be safely used with a wider range of CUDA calls.
 *
 * @tparam T type of elements to be stored in the array
 */
template<typename T>
class PinnedArray1D : public Array1D<T> {
public:
    /**
     * Create an empty 1-D array. No memory will be allocated.
     */
    PinnedArray1D() : Array1D<T>() {}

    /**
     * Create a 1-D array of the requested size, using pinned memory allocators.
     * Array's values won't be initialized.
     *
     * @param length number of elements to be stored in array
     */
    explicit PinnedArray1D(index_t length)
            : Array1D<T>(length, cuda_host_alloc<T>, cuda_host_free) {}
};

/**
 * Subclass of Array1D template for allocating two-dimensional arrays
 * using page-locked (pinned) memory which can be safely used with a wider range of CUDA calls.
 *
 * @tparam T type of elements to be stored in the array
 */
template<typename T>
class PinnedArray2D : public Array2D<T> {
public:
    /**
     * Create an empty 2-D array. No memory will be allocated.
     */
    PinnedArray2D() : Array2D<T>() {}

    /**
     * Create a 2-D array of the requested dimensions, using pinned memory allocators.
     * Array's values won't be initialized.
     *
     * @param height first dimension of the array (number of sub-arrays)
     * @param length number of elements to be stored in each sub-array
     */
    explicit PinnedArray2D(int height, index_t length)
            : Array2D<T>(height, length, cuda_host_alloc<T>, cuda_host_free) {}
};

#endif //EMPI_PINNED_ARRAY_H

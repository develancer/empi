/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_ARRAY_H
#define EMPI_ARRAY_H

#include <cstddef>
#include <memory>
#include "Types.h"

template<typename T>
T *default_creator(size_t length) {
    return new T[length];
}

template<typename T>
void default_deleter(T *pointer) {
    delete[] pointer;
}

template<typename T, typename C>
T *safe_create(C creator, index_t length) {
    if (length < 0) {
        throw std::logic_error("array of negative length requested");
    }
    T *created = creator(length);
    if (!created) {
        throw std::logic_error("failed to allocate array");
    }
    return created;
}

/**
 * Wrapper for one-dimensional array of constant size selected at runtime.
 * Internally based on std::shared_ptr.
 *
 * @tparam T type of elements to be stored in the array
 */
template<typename T>
class Array1D : public std::shared_ptr<T[]> {
    index_t length_;

public:
    /**
     * Create an empty 1-D array. No memory will be allocated.
     */
    Array1D() : std::shared_ptr<T[]>(), length_(0) {}

    /**
     * Create a 1-D array of the requested size, using standard allocators.
     * Array's values won't be initialized.
     *
     * @param length number of elements to be stored in array
     */
    explicit Array1D(index_t length)
            : std::shared_ptr<T[]>(safe_create<T>(default_creator<T>, length), default_deleter<T>), length_(length) {}

    /**
     * Create a 1-D array of the requested size, using custom allocators.
     * Array's values won't be initialized, unless custom allocator does so.
     *
     * @param length number of elements to be stored in array
     * @param creator custom allocator taking number of elements as a parameter, and returning T*
     * @param deleter custom deleter taking T* as a parameter
     */
    template<typename C, typename D>
    explicit Array1D(index_t length, C creator, D deleter)
            : std::shared_ptr<T[]>(safe_create<T, C>(creator, length), deleter), length_(length) {}

    /**
     * @return length (number of elements) of the array
     */
    [[nodiscard]] index_t length() const {
        return length_;
    }

    /**
     * Set all items in the array to a given value.
     */
    void fill(T value) {
        T *const pointer = this->get();
        std::fill(pointer, pointer + length_, value);
    }
};

template<typename T>
class Array2D : public std::shared_ptr<T *const[]> {
    int height_;
    index_t length_;

    template<typename C, typename D>
    static std::shared_ptr<T *[]> allocate(int height, index_t length, C creator, D deleter) {
        // no const at this step
        auto result = std::shared_ptr<T *[]>(
                new T *[height](),
                [height, deleter](T **pointer) {
                    for (int h = 0; h < height; ++h) {
                        if (pointer[h]) {
                            deleter(pointer[h]);
                        }
                    }
                    delete[] pointer;
                }
        );
        for (int h = 0; h < height; ++h) {
            result[h] = safe_create<T, C>(creator, length);
        }
        return result;
    }

public:
    /**
     * Create an empty 2-D array. No memory will be allocated.
     */
    Array2D() : std::shared_ptr<T *const[]>(), height_(0), length_(0) {}

    /**
     * Create a 2-D array of the requested dimensions, using standard allocators.
     * Array's values won't be initialized.
     *
     * @param height first dimension of the array (number of sub-arrays)
     * @param length number of elements to be stored in each sub-array
     */
    explicit Array2D(int height, index_t length)
            : std::shared_ptr<T *const[]>(allocate(height, length, default_creator<T>, default_deleter<T>)),
              height_(height), length_(length) {}

    /**
     * Create a 2-D array of the requested dimensions, using custom allocators.
     * Array's values won't be initialized, unless custom allocator does so.
     *
     * @param height first dimension of the array (number of sub-arrays)
     * @param length number of elements to be stored in each sub-array
     */
    template<typename C, typename D>
    explicit Array2D(int height, index_t length, C creator, D deleter)
            : std::shared_ptr<T *const[]>(allocate(height, length, creator, deleter)),
              height_(height), length_(length) {}

    /**
     * @return number of sub-arrays in the array
     */
    [[nodiscard]] inline int height() const {
        return height_;
    }

    /**
     * @return length (number of elements) of each sub-array
     */
    [[nodiscard]] inline index_t length() const {
        return length_;
    }

    /**
     * Set all items in all sub-arrays to a given value.
     */
    void fill(T value) {
        T *const *const pointer = this->get();
        for (int h = 0; h < height_; ++h) {
            std::fill(pointer[h], pointer[h] + length_, value);
        }
    }
};

#endif //EMPI_ARRAY_H

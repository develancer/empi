/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_TYPES_H
#define EMPI_TYPES_H

#include <cassert>
#include <cfenv>
#include <complex>
#include <cstddef>
#include <limits>

using index_t = std::ptrdiff_t;
using real = double; // TODO #ifdef SINGLE_PRECISION
using complex = std::complex<double>;

template<typename X, typename Y>
inline size_t mulx(X x, Y y) {
    // TODO range checking
    return static_cast<size_t>(x) * static_cast<size_t>(y);
}

template<typename T>
int as_positive_int(T x) {
    assert(x > 0);
    assert(static_cast<size_t>(x) <= static_cast<size_t>(std::numeric_limits<int>::max()));
    return static_cast<int>(x);
}

class Types {
    template<typename T>
    static T rint(double x);

public:
    template<typename T>
    static T round(double x) {
        if (!(x > static_cast<double>(std::numeric_limits<T>::lowest()) && x < static_cast<double>(std::numeric_limits<T>::max()))) {
            throw std::runtime_error("rounding overflow detected");
        }
        return rint<T>(x);
    }

private:
    template<typename T>
    static T round_with_mode(double x, int rounding_mode) {
        const int default_mode = fegetround();
        fesetround(rounding_mode);
        T result = round<T>(x);
        fesetround(default_mode);
        return result;
    }

public:
    template<typename T>
    static T ceil(double x) {
        return round_with_mode<T>(x, FE_UPWARD);
    }

    template<typename T>
    static T floor(double x) {
        return round_with_mode<T>(x, FE_DOWNWARD);
    }
};

template<>
int Types::rint(double x);

template<>
long Types::rint(double x);

template<>
long long Types::rint(double x);

#endif //EMPI_TYPES_H

/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cmath>
#include "Types.h"

template<>
int Types::rint(double x) {
    return static_cast<int>(std::lrint(x));
}

template<>
long Types::rint(double x) {
    return std::lrint(x);
}

template<>
long long Types::rint(double x) {
    return std::llrint(x);
}

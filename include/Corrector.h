/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_CORRECTOR_H
#define EMPI_CORRECTOR_H

#include "ExtraData.h"
#include "Types.h"

#ifdef _WIN32
__declspec(align(16))
#endif
class Corrector {
    complex ft;
    complex re_factor;
    complex im_factor;

public:
    explicit Corrector();

    explicit Corrector(complex envelope_ft);

    [[nodiscard]] inline double estimate_energy(complex value) const {
        return std::norm(value.real() * re_factor + value.imag() * im_factor);
    }

    double compute(complex value, ExtraData* extra = nullptr) const;

}
#ifndef  _WIN32
__attribute__ ((aligned (16)))
#endif
;

#endif //EMPI_CORRECTOR_H

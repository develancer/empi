/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "Corrector.h"

//////////////////////////////////////////////////////////////////////////////

Corrector::Corrector() : ft(NAN), re_factor(NAN), im_factor(NAN) {}

Corrector::Corrector(complex envelope_ft) : ft(envelope_ft) {
    const double norm_ft = std::norm(ft);
    if (norm_ft + 1.0e-8 >= 1) {
        re_factor = 2.0 / (1 + std::sqrt(norm_ft));
        im_factor = re_factor * complex(0.0, 1.0);
    } else {
        const double common_factor = 2 / (1 - norm_ft);
        re_factor = (1.0 - ft) * common_factor;
        im_factor = (1.0 + ft) * common_factor * complex(0.0, 1.0);
    }
}

double Corrector::compute(complex value, ExtraData* extra) const {
    const complex corrected = value.real() * re_factor + value.imag() * im_factor;
    const complex corrected2 = corrected * corrected;
    const double norm_corrected = std::norm(corrected);
    const double energy = 0.5 * (norm_corrected + corrected2.real() * ft.real() + corrected2.imag() * ft.imag());
    if (extra) {
        extra->amplitude = std::sqrt(norm_corrected);
        extra->energy = energy;
        extra->phase = std::arg(corrected);
    }
    return energy;
}

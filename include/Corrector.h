#ifndef EMPI_CORRECTOR_H
#define EMPI_CORRECTOR_H

#include "Types.h"

class CorrectorResult {
    complex ft;
    complex corrected;
    double norm_corrected;

public:
    CorrectorResult(complex ft, complex corrected) : ft(ft), corrected(corrected), norm_corrected(std::norm(corrected)) {}

    // result has to be multiplied by family->value(0.0) / std::sqrt(scale);
    double amplitude() const {
        return std::sqrt(norm_corrected);
    }

    double energy() const {
        complex corrected2 = corrected * corrected;
        return 0.5 * (norm_corrected + corrected2.real() * ft.real() + corrected2.imag() * ft.imag());
    }

    complex modulus() const {
        return corrected * std::sqrt(energy() / norm_corrected);
    }

    double phase() const {
        return std::arg(corrected);
    }
};

class Corrector {
    complex ft;
    double common_factor;
    double estimate_factor;

public:
    explicit Corrector() : ft(), common_factor(NAN), estimate_factor(NAN) {}

    explicit Corrector(complex envelope_ft) : ft(envelope_ft), common_factor(2 / (1 - std::norm(ft))) {
        const double norm_ft = std::norm(ft);
        if (norm_ft + 1.0e-6 >= 1) {
            estimate_factor = 1.0;
        } else {
            estimate_factor = 2 / (1 - std::sqrt(norm_ft));
        }
    }

    double estimate_energy(complex value) const {
        return estimate_factor * std::norm(value);
    }

    CorrectorResult compute(complex value) const {
        const bool is_value_real = std::abs(value.imag()) <= 1.0e-10 * std::abs(value.real());
        const complex corrected = is_value_real ? value * 2.0 / (1.0 + ft) :
                (value - std::conj(value) * ft) * common_factor;
        return CorrectorResult(ft, corrected);
    }
} __attribute__ ((aligned (16)));

#endif //EMPI_CORRECTOR_H

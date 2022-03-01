/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_SPECTROGRAM_CALCULATOR_DUMMY_H
#define EMPI_SPECTROGRAM_CALCULATOR_DUMMY_H

#include "SpectrogramCalculator.h"

/**
 * Dummy short-time Fourier transform calculator, used only for testing purposes.
 */
class SpectrogramCalculatorDummy : public SpectrogramCalculator {
public:
    /**
     * Compute a spectrogram of a given signal, according to the given specification.
     *
     * @param request specification for the spectrogram to be computed
     */
    void compute(const SpectrogramRequest &request) final;
};

#endif //EMPI_SPECTROGRAM_CALCULATOR_DUMMY_H

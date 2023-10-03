/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_SPECTRUM_CALCULATOR_H
#define EMPI_SPECTRUM_CALCULATOR_H

#include "Array.h"
#include "Types.h"

/**
 * Interface for objects calculating spectra (discrete Fourier transforms) of real signals.
 */
class SpectrumCalculator {
public:
    /**
     * Compute spectrum of a given signal with optional zero-padding.
     * Return the pointer to the internal buffer which may be overwritten on next calls to computeSpectrum.
     *
     * @param input input signal
     * @param window_length length of the window for FFT, should be at least the size of input
     * @return pointer to the internal buffer with size of at last window_length/2+1
     */
    virtual const complex *computeSpectrum(Array1D<double> input, int window_length) = 0;

    virtual ~SpectrumCalculator() = default;
};

#endif //EMPI_SPECTRUM_CALCULATOR_H

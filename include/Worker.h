/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_WORKER_H
#define EMPI_WORKER_H

#include "SpectrogramRequest.h"

/**
 * Interface for objects calculating spectrograms (short-time Fourier transforms) of real signals.
 */
class Worker {
public:
    /**
     * Compute a spectrogram of a given signal, according to the given specification.
     *
     * @param request specification for the spectrogram to be computed
     */
    virtual void compute(const SpectrogramRequest &request) = 0;

    virtual ~Worker() = default;
};

#endif //EMPI_WORKER_H

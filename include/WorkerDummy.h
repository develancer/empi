/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_WORKER_DUMMY_H
#define EMPI_WORKER_DUMMY_H

#include "Worker.h"

/**
 * Dummy short-time Fourier transform calculator, used only for testing purposes.
 */
class WorkerDummy : public Worker {
public:
    /**
     * Compute a spectrogram of a given signal, according to the given specification.
     *
     * @param request specification for the spectrogram to be computed
     */
    void compute(const SpectrogramRequest &request) final;
};

#endif //EMPI_WORKER_DUMMY_H

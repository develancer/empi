/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_WORKER_CUDA_H
#define EMPI_WORKER_CUDA_H

#include <list>
#include <memory>
#include "Worker.h"

class CudaTask;

/**
 * Short-time Fourier transform calculator based on CUDA and CuFFT.
 */
class WorkerCUDA : public Worker {
    std::shared_ptr<CudaTask> task;

public:
    /**
     * Create a new CUDA worker for multi-channel signals with given number of channels,
     * able to serve pre-defined set of request templates.
     *
     * @param channel_count number of channels in the signal
     * @param proto_requests list of request templates that will be passed to compute() later on
     */
    WorkerCUDA(int channel_count, const std::list<ProtoRequest>& proto_requests);

    /**
     * Compute a spectrogram of a given signal, according to the given specification.
     *
     * @param request specification for the spectrogram to be computed
     */
    void compute(const SpectrogramRequest &request) final;
};

#endif //EMPI_WORKER_CUDA_H

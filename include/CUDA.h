#ifndef EMPI_CUDA_H
#define EMPI_CUDA_H

#ifndef __CUDACC__
#error this file can be included only when compiling with nvcc
#endif

#include <cstdio>
#include <stdexcept>
#include <cufftXt.h>
#include "Types.h"

using cucomplex = cuDoubleComplex; // TODO #ifdef SINGLE_PRECISION

struct CudaCallbackInfo {
    real *envelope;
    void *correctors;
    uint32_t window_length;
    uint32_t spectrum_length;
    uint32_t envelope_length;
    uint32_t input_shift;
    uint32_t output_bins;
    unsigned window_length_bits;
    size_t window_length_mask;
};

class CudaCallback {
    void *hostCopyOfInputCallback = nullptr;
    void *hostCopyOfOutputCallback = nullptr;

public:
    void initialize();

    void associate(cufftHandle plan, CudaCallbackInfo *dev_info);
};

class CudaException : public std::runtime_error {
    static char buffer[256];

    template<typename... Args>
    const char *prepare(const char *format, Args... args) {
        snprintf(buffer, sizeof buffer, format, args...);
        return buffer;
    }

public:
    explicit CudaException(cudaError_t error)
            : std::runtime_error(prepare("GPU error (CUDA) %s (%s)", cudaGetErrorName(error), cudaGetErrorString(error))) {}

    explicit CudaException(cufftResult_t result)
            : std::runtime_error(prepare("GPU error (cuFFT) #%d", static_cast<int>(result))) {}
};

class CudaMemoryException : public std::exception {
    const char *what() const noexcept override {
        return "GPU ran out of memory";
    }
};

void cuda_check(cudaError_t error);

void cufft_check(cufftResult_t result);

#endif //EMPI_CUDA_H

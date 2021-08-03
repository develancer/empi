/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>
#include "Array.h"
#include "Corrector.h"
#include "CUDA.h"
#include "IndexRange.h"
#include "PinnedArray.h"
#include "WorkerCUDA.h"

char CudaException::buffer[256];

void cuda_check(cudaError_t error) {
    if (error) {
        if (error == cudaErrorMemoryAllocation) {
            throw CudaMemoryException();
        }
        throw CudaException(error);
    }
}

void cufft_check(cufftResult_t result) {
    if (result) {
        if (result == CUFFT_ALLOC_FAILED) {
            throw CudaMemoryException();
        }
        throw CudaException(result);
    }
}

//////////////////////////////////////////////////////////////////////////////

class CudaCorrectorResult {
    cucomplex ft;
    cucomplex corrected;
    double norm_corrected;

public:
    __device__ CudaCorrectorResult(cucomplex ft, cucomplex corrected) : ft(ft), corrected(corrected),
                                                                        norm_corrected(corrected.x * corrected.x + corrected.y * corrected.y) {}

    // result has to be multiplied by family->value(0.0) / std::sqrt(scale);
    __device__ double amplitude() const {
        return sqrt(norm_corrected);
    }

    __device__ double energy() const {
        cucomplex corrected2 = cuCmul(corrected, corrected);
        return 0.5 * (norm_corrected + corrected2.x * ft.x + corrected2.y * ft.y);
    }

    __device__ cucomplex modulus() const {
        double factor = sqrt(energy() / norm_corrected);
        return make_cuDoubleComplex(corrected.x * factor, corrected.y * factor);
    }
};

//////////////////////////////////////////////////////////////////////////////

class CudaCorrector {
    cucomplex ft;
    double common_factor;
    double estimate_factor;

    CudaCorrector() = delete;

public:
    __device__ CudaCorrectorResult compute(cucomplex value) const {
        const bool is_value_real = std::abs(cuCimag(value)) <= 1.0e-10 * std::abs(cuCreal(value));
        if (is_value_real) {
            return CudaCorrectorResult(ft, cuCdiv(value, make_cuDoubleComplex(0.5 + cuCreal(ft) / 2, cuCimag(ft) / 2)));
        }
        const cucomplex corrected = cuCsub(value, cuCmul(cuConj(value), ft));
        return CudaCorrectorResult(ft, make_cuDoubleComplex(corrected.x * common_factor, corrected.y * common_factor));
    }
} __attribute__ ((aligned (16)));

//////////////////////////////////////////////////////////////////////////////

static __device__ __inline__ void oneReduce(volatile ExtractedMaximum &dest, const volatile ExtractedMaximum &other) {
    if (other.energy > dest.energy) {
        dest.energy = other.energy;
        dest.bin_index = other.bin_index;
    }
}

template<unsigned THREADS_IN_BLOCK>
static __device__ __inline__ void warpReduce(volatile ExtractedMaximum *sdata, unsigned tid) {
    if (THREADS_IN_BLOCK >= 64) oneReduce(sdata[tid], sdata[tid + 32]);
    if (THREADS_IN_BLOCK >= 32) oneReduce(sdata[tid], sdata[tid + 16]);
    if (THREADS_IN_BLOCK >= 16) oneReduce(sdata[tid], sdata[tid + 8]);
    if (THREADS_IN_BLOCK >= 8) oneReduce(sdata[tid], sdata[tid + 4]);
    if (THREADS_IN_BLOCK >= 4) oneReduce(sdata[tid], sdata[tid + 2]);
    if (THREADS_IN_BLOCK >= 2) oneReduce(sdata[tid], sdata[tid + 1]);
}

template<unsigned THREADS_IN_BLOCK>
static __global__ void
extractorSingleChannelCUDA(unsigned channel_count, unsigned output_bins, cucomplex **spectra,
                           ExtractedMaximum *maxima, CudaCorrector *correctors) {
    __shared__ ExtractedMaximum sdata[THREADS_IN_BLOCK];
    const unsigned int tid = threadIdx.x;
    sdata[tid].energy = 0.0;
    sdata[tid].bin_index = 0;
    assert(channel_count == 1);

    unsigned offset = blockIdx.x * output_bins;
    for (unsigned i = tid; i < output_bins; i += THREADS_IN_BLOCK) {
        // tutaj każdy wątek zbiera po wszystkich elementach którymi ma się zająć
        real energy = correctors[i].compute(spectra[0][offset + i]).energy();
        if (energy > sdata[tid].energy) {
            sdata[tid].energy = energy;
            sdata[tid].bin_index = i;
        }
    }
    __syncthreads();

    // tutaj mamy już tylko tyle elementów ile wątków
    if (THREADS_IN_BLOCK >= 512) {
        if (tid < 256) oneReduce(sdata[tid], sdata[tid + 256]);
        __syncthreads();
    }
    if (THREADS_IN_BLOCK >= 256) {
        if (tid < 128) oneReduce(sdata[tid], sdata[tid + 128]);
        __syncthreads();
    }
    if (THREADS_IN_BLOCK >= 128) {
        if (tid < 64) oneReduce(sdata[tid], sdata[tid + 64]);
        __syncthreads();
    }
    if (tid < 32) warpReduce<THREADS_IN_BLOCK>(sdata, tid);
    if (tid == 0) {
        maxima[blockIdx.x].bin_index = sdata[0].bin_index;
        maxima[blockIdx.x].energy = sdata[0].energy;
    }
}

template<unsigned THREADS_IN_BLOCK>
static __global__ void
extractorConstantPhaseCUDA(unsigned channel_count, unsigned output_bins, cucomplex **spectra,
                           ExtractedMaximum *maxima, CudaCorrector *correctors) {
    __shared__ ExtractedMaximum sdata[THREADS_IN_BLOCK];
    const unsigned int tid = threadIdx.x;
    sdata[tid].energy = 0.0;
    sdata[tid].bin_index = 0;

    unsigned offset = blockIdx.x * output_bins;
    for (unsigned i = tid; i < output_bins; i += THREADS_IN_BLOCK) {
        // tutaj każdy wątek zbiera po wszystkich elementach którymi ma się zająć
        cucomplex direction = make_cuDoubleComplex(0, 0);
        for (unsigned c = 0; c < channel_count; ++c) {
            direction = cuCfma(spectra[c][offset + i], spectra[c][offset + i], direction);
        }
        double angle = 0.5 * atan2(cuCimag(direction), cuCreal(direction));
        double sin_angle, cos_angle;
        sincos(angle, &sin_angle, &cos_angle);
        auto corrector_result = correctors[i].compute(make_cuDoubleComplex(cos_angle, sin_angle));

        double energy = 0;
        for (unsigned c = 0; c < channel_count; ++c) {
            const cucomplex v = spectra[c][offset + i];
            double cos_factor = cos(atan2(cuCimag(v), cuCreal(v)) - angle);
            energy += (cuCreal(v) * cuCreal(v) + cuCimag(v) * cuCimag(v)) * cos_factor * cos_factor;
        }
        energy *= corrector_result.energy();
        if (energy > sdata[tid].energy) {
            sdata[tid].energy = energy;
            sdata[tid].bin_index = i;
        }
    }
    __syncthreads();

    // tutaj mamy już tylko tyle elementów ile wątków
    if (THREADS_IN_BLOCK >= 512) {
        if (tid < 256) oneReduce(sdata[tid], sdata[tid + 256]);
        __syncthreads();
    }
    if (THREADS_IN_BLOCK >= 256) {
        if (tid < 128) oneReduce(sdata[tid], sdata[tid + 128]);
        __syncthreads();
    }
    if (THREADS_IN_BLOCK >= 128) {
        if (tid < 64) oneReduce(sdata[tid], sdata[tid + 64]);
        __syncthreads();
    }
    if (tid < 32) warpReduce<THREADS_IN_BLOCK>(sdata, tid);
    if (tid == 0) {
        maxima[blockIdx.x].bin_index = sdata[0].bin_index;
        maxima[blockIdx.x].energy = sdata[0].energy;
    }
}

template<unsigned THREADS_IN_BLOCK>
static __global__ void
extractorVariablePhaseCUDA(unsigned channel_count, unsigned output_bins, cucomplex **spectra,
                           ExtractedMaximum *maxima, CudaCorrector *correctors) {
    __shared__ ExtractedMaximum sdata[THREADS_IN_BLOCK];
    const unsigned int tid = threadIdx.x;
    sdata[tid].energy = 0.0;
    sdata[tid].bin_index = 0;

    unsigned offset = blockIdx.x * output_bins;
    for (unsigned i = tid; i < output_bins; i += THREADS_IN_BLOCK) {
        // tutaj każdy wątek zbiera po wszystkich elementach którymi ma się zająć
        real energy = 0;
        for (unsigned c = 0; c < channel_count; ++c) {
            const cucomplex v = spectra[c][offset + i];
            energy += correctors[i].compute(v).energy();
        }
        if (energy > sdata[tid].energy) {
            sdata[tid].energy = energy;
            sdata[tid].bin_index = i;
        }
    }
    __syncthreads();

    // tutaj mamy już tylko tyle elementów ile wątków
    if (THREADS_IN_BLOCK >= 512) {
        if (tid < 256) oneReduce(sdata[tid], sdata[tid + 256]);
        __syncthreads();
    }
    if (THREADS_IN_BLOCK >= 256) {
        if (tid < 128) oneReduce(sdata[tid], sdata[tid + 128]);
        __syncthreads();
    }
    if (THREADS_IN_BLOCK >= 128) {
        if (tid < 64) oneReduce(sdata[tid], sdata[tid + 64]);
        __syncthreads();
    }
    if (tid < 32) warpReduce<THREADS_IN_BLOCK>(sdata, tid);
    if (tid == 0) {
        maxima[blockIdx.x].bin_index = sdata[0].bin_index;
        maxima[blockIdx.x].energy = sdata[0].energy;
    }
}

void callExtractorKernel(const SpectrogramRequest &request, cudaStream_t stream, cucomplex **spectra,
                         ExtractedMaximum *maxima, Corrector *correctors) {
    const unsigned threads_in_block = 64;
    if (request.extractor == extractorVariablePhase) {
        extractorVariablePhaseCUDA<threads_in_block><<<request.how_many, threads_in_block, 0, stream>>>
        (request.channel_count, request.output_bins, spectra, maxima, reinterpret_cast<CudaCorrector *>(correctors));
    } else if (request.extractor == extractorConstantPhase) {
        extractorConstantPhaseCUDA<threads_in_block><<<request.how_many, threads_in_block, 0, stream>>>
        (request.channel_count, request.output_bins, spectra, maxima, reinterpret_cast<CudaCorrector *>(correctors));
    } else if (request.extractor == extractorSingleChannel) {
        extractorSingleChannelCUDA<threads_in_block><<<request.how_many, threads_in_block, 0, stream>>>
                (request.channel_count, request.output_bins, spectra, maxima, reinterpret_cast<CudaCorrector *>(correctors));
    } else {
        fprintf(stderr, "ERROR: unsupported extractor\n");
        exit(1);
    }
}

void *cuda_dev_malloc(size_t length) {
    void *result;
    cuda_check(cudaMalloc(&result, length));
    return result;
}

void cuda_dev_free(void *pointer) {
    cuda_check(cudaFree(pointer));
}

template<typename T>
T *cuda_dev_alloc(size_t length) {
    return reinterpret_cast<T *>(cuda_dev_malloc(sizeof(T) * length));
}

template<typename T>
Array1D<T> cuda_dev_array_1d(size_t n) {
    return Array1D<T>(n, cuda_dev_alloc<T>, cuda_dev_free);
}

template<typename T>
Array2D<T> cuda_dev_array_2d(size_t n, size_t m) {
    return Array2D<T>(n, m, cuda_dev_alloc<T>, cuda_dev_free);
}

std::shared_ptr<CUstream_st> cuda_create_stream() {
    cudaStream_t stream;
    cuda_check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    return std::shared_ptr<CUstream_st>(stream, cudaStreamDestroy);
}

//////////////////////////////////////////////////////////////////////////////

class CudaTask {
    std::map<std::pair<int, int>, cufftHandle> plans;

    std::shared_ptr<CudaCallbackInfo> info;
    std::shared_ptr<CudaCallbackInfo> dev_info;
    Array1D<real> dev_envelope;
    Array1D<Corrector> dev_correctors;

    std::shared_ptr<CUstream_st> stream;

    Array1D<char> dev_workspace;
    Array1D<real> dev_input;
    Array2D<cucomplex> outputs;
    Array1D<cucomplex *> dev_outputs;
    Array1D<ExtractedMaximum> dev_maxima;

    void process(const SpectrogramRequest &request, cufftHandle handle) {
        index_t total_input_length = request.envelope_length + (request.how_many - 1) * request.input_shift;
        IndexRange overlap = IndexRange(-request.input_offset, request.channel_length - request.input_offset).overlap(total_input_length);
        if (!overlap) {
            // all samples consist of zero-padding
            std::fill(request.maxima, request.maxima + request.how_many, ExtractedMaximum{0, 0});
            return;
        }

        for (int c = 0; c < request.channel_count; ++c) {
            if (overlap.first_index) {
                cuda_check(cudaMemsetAsync(dev_input.get(), 0, overlap.first_index * sizeof(real), stream.get()));
            }
            cuda_check(cudaMemcpyAsync(dev_input.get() + overlap.first_index,
                                       request.data[c] + request.input_offset + overlap.first_index,
                                       (overlap.end_index - overlap.first_index) * sizeof(real), cudaMemcpyHostToDevice,
                                       stream.get()));
            if (overlap.end_index < total_input_length) {
                cuda_check(cudaMemsetAsync(dev_input.get() + overlap.end_index, 0,
                                           (total_input_length - overlap.end_index) * sizeof(real), stream.get()));
            }
            cufft_check(cufftExecD2Z(handle, dev_input.get(), outputs[c]));
        }
        callExtractorKernel(request, stream.get(), dev_outputs.get(), dev_maxima.get(), dev_correctors.get());
        cuda_check(cudaMemcpyAsync(request.maxima, dev_maxima.get(), request.how_many * sizeof(ExtractedMaximum),
                                   cudaMemcpyDeviceToHost, stream.get()));
    }

public:
    CudaTask(int channel_count, const std::list<ProtoRequest> &proto_requests) {
        int max_envelope_length = 0;
        int max_output_bins = 0;
        int max_window_length = 0;
        for (const auto &proto_request : proto_requests) {
            max_envelope_length = std::max(max_envelope_length, proto_request.envelope_length);
            max_output_bins = std::max(max_output_bins, proto_request.output_bins);
            max_window_length = std::max(max_window_length, proto_request.window_length);
        }

        stream = cuda_create_stream();

        info.reset(cuda_host_alloc<CudaCallbackInfo>(1), cuda_host_free);
        dev_info.reset(cuda_dev_alloc<CudaCallbackInfo>(1), cuda_dev_free);
        dev_envelope = cuda_dev_array_1d<real>(max_envelope_length);
        dev_correctors = cuda_dev_array_1d<Corrector>(max_output_bins);
        dev_outputs = cuda_dev_array_1d<cucomplex *>(channel_count);

        size_t free_memory, total_memory;
        cuda_check(cudaMemGetInfo(&free_memory, &total_memory));

        size_t max_dev_input_length = 0;
        size_t max_dev_output_length = 0;
        size_t max_work_size = 0;
        int max_how_many = 0;
        std::map<int, std::set<int>> how_many_for_various_windows;
        for (const auto &proto_request : proto_requests) {
            int window_length = proto_request.window_length;
            int spectrum_length = window_length / 2 + 1;

            for (int how_many = 1; how_many <= proto_request.how_many; how_many *= 2) {
                size_t dev_input_length = static_cast<size_t>(proto_request.envelope_length)
                                          + static_cast<size_t>(how_many - 1) * static_cast<size_t>(proto_request.input_shift);
                size_t dev_output_length = static_cast<size_t>(how_many) * static_cast<size_t>(proto_request.output_bins);
                size_t work_size;
                cufft_check(cufftEstimateMany(
                        1, &window_length,
                        &window_length, 1, window_length,
                        &spectrum_length, 1, spectrum_length,
                        CUFFT_D2Z, how_many, &work_size
                ));

                size_t memory_needed = work_size
                                       + sizeof(double) * dev_input_length
                                       + sizeof(cucomplex) * static_cast<size_t>(channel_count) * dev_output_length
                                       + sizeof(ExtractedMaximum) * static_cast<size_t>(how_many);
                if (memory_needed > free_memory) {
                    if (how_many == 1) {
                        throw std::runtime_error("not enough CUDA memory");
                    }
                    break;
                }

                max_dev_input_length = std::max(max_dev_input_length, dev_input_length);
                max_dev_output_length = std::max(max_dev_output_length, dev_output_length);
                max_work_size = std::max(max_work_size, work_size);
                max_how_many = std::max(max_how_many, how_many);

                how_many_for_various_windows[window_length].insert(how_many);
            }
        }

        dev_input = cuda_dev_array_1d<double>(max_dev_input_length);

        outputs = cuda_dev_array_2d<cucomplex>(channel_count, max_dev_output_length);
        cuda_check(cudaMemcpy(dev_outputs.get(), outputs.get(), channel_count * sizeof(cucomplex *), cudaMemcpyHostToDevice));

        dev_maxima = cuda_dev_array_1d<ExtractedMaximum>(max_how_many);

        max_work_size = 0; // this will now be calculated more precisely
        for (const auto &pair : how_many_for_various_windows) {
            int window_length = pair.first;
            int spectrum_length = window_length / 2 + 1;
            for (int how_many : pair.second) {
                cufftHandle handle;
                size_t work_size = 0;
                cufft_check(cufftCreate(&handle));
                cufft_check(cufftSetAutoAllocation(handle, 0));
                cufft_check(cufftMakePlanMany(
                        handle, 1, &window_length,
                        &window_length, 1, window_length,
                        &spectrum_length, 1, spectrum_length,
                        CUFFT_D2Z, how_many, &work_size
                ));
                max_work_size = std::max(max_work_size, work_size);
                associateCallbackWithPlan(handle, dev_info.get());
                cufft_check(cufftSetStream(handle, stream.get()));
                plans.emplace(std::make_pair(window_length, how_many), std::move(handle));
            }
        }

        dev_workspace = cuda_dev_array_1d<char>(max_work_size);
        for (auto &pair : plans) {
            cufft_check(cufftSetWorkArea(pair.second, dev_workspace.get()));
        }
    }

    ~CudaTask() {
        for (auto &pair : plans) {
            cufftDestroy(pair.second);
        }
    }

    void compute(const SpectrogramRequest &request) {
        const index_t total_input_length = request.envelope_length + (request.how_many - 1) * request.input_shift;
        IndexRange overlap = IndexRange(-request.input_offset, request.channel_length - request.input_offset).overlap(total_input_length);
        if (!overlap) {
            // all samples consist of zero-padding
            std::fill(request.maxima, request.maxima + request.how_many, ExtractedMaximum{0, 0});
            return;
        }

        info->window_length = request.window_length;
        info->spectrum_length = request.window_length / 2 + 1;
        info->envelope_length = request.envelope_length;
        info->input_shift = request.input_shift;
        info->output_bins = request.output_bins;
        info->envelope = dev_envelope.get();
        info->correctors = dev_correctors.get();
        // only supported for powers of 2
        for (info->window_length_bits = 0, info->window_length_mask = 1;
             request.window_length != info->window_length_mask; ++info->window_length_bits) {
            info->window_length_mask <<= 1;
        }
        info->window_length_mask--;

        cuda_check(cudaMemcpyAsync(dev_envelope.get(), request.envelope, request.envelope_length * sizeof(real),
                                   cudaMemcpyHostToDevice, stream.get()));
        cuda_check(cudaMemcpyAsync(dev_correctors.get(), request.correctors,
                                   request.output_bins * sizeof(Corrector), cudaMemcpyHostToDevice, stream.get()));
        cuda_check(cudaMemcpyAsync(dev_info.get(), info.get(), sizeof(CudaCallbackInfo), cudaMemcpyHostToDevice, stream.get()));

        SpectrogramRequest part(request);
        int how_many = 1, how_many_left = request.how_many;
        while (how_many <= how_many_left) {
            how_many *= 2;
        }
        assert(how_many > how_many_left);
        while (how_many_left > 0 && (how_many /= 2) > 0) {
            auto plan_iterator = plans.find({request.window_length, how_many});
            if (plan_iterator != plans.end()) {
                while (how_many_left >= how_many) {
                    cufftHandle handle = plan_iterator->second;
                    part.how_many = how_many;

                    process(part, handle);

                    part.input_offset += static_cast<index_t>(request.input_shift) * static_cast<index_t>(how_many);
                    part.maxima += how_many;
                    how_many_left -= how_many;
                }
            }
        }
        cuda_check(cudaStreamSynchronize(stream.get()));
    }
};

//////////////////////////////////////////////////////////////////////////////

WorkerCUDA::WorkerCUDA(int channel_count, const std::list<ProtoRequest> &proto_requests) {
    task = std::make_shared<CudaTask>(channel_count, proto_requests);
}

void WorkerCUDA::compute(const SpectrogramRequest &request) {
    request.assertCorrectness();
    task->compute(request);
}

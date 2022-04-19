/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "CUDA.h"

static __device__  real myInputCallback(void *dataIn, size_t offset, void *callerInfo, void *) {
    auto *info = reinterpret_cast<CudaCallbackInfo *>(callerInfo);
    size_t index_in_batch = offset & info->window_length_mask;
    real ret = (index_in_batch < info->envelope_length)
               ? ((real *) dataIn)[(offset >> info->window_length_bits) * info->input_shift + index_in_batch] * info->envelope[index_in_batch]
               : 0;
    return ret;
}

static __device__  void
myOutputCallback(void *dataOut, size_t offset, cucomplex element, void *callerInfo, void *) {
    auto *info = reinterpret_cast<CudaCallbackInfo *>(callerInfo);
    size_t index_in_batch = offset % info->spectrum_length;
    if (index_in_batch < info->output_bins) {
        ((cucomplex *) dataOut)[offset / info->spectrum_length * info->output_bins + index_in_batch] = element;
    }
}

static __device__ cufftCallbackLoadD myInputPtr = myInputCallback;

static __device__ cufftCallbackStoreZ myOutputPtr = myOutputCallback;

void CudaCallback::initialize() {
    cuda_check(cudaMemcpyFromSymbol(&hostCopyOfInputCallback, myInputPtr, sizeof(void*)));
    cuda_check(cudaMemcpyFromSymbol(&hostCopyOfOutputCallback, myOutputPtr, sizeof(void*)));
}

void CudaCallback::associate(cufftHandle plan, CudaCallbackInfo *dev_info) {
    cufft_check(cufftXtSetCallback(plan, (void **) &hostCopyOfInputCallback, CUFFT_CB_LD_REAL_DOUBLE, (void **) &dev_info));
    cufft_check(cufftXtSetCallback(plan, (void **) &hostCopyOfOutputCallback, CUFFT_CB_ST_COMPLEX_DOUBLE, (void **) &dev_info));
}

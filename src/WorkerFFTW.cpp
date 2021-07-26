/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <list>
#include <fftw3.h>
#include "Corrector.h"
#include "IndexRange.h"
#include "WorkerFFTW.h"

template<typename TX, typename TY>
static void mul(size_t count, TX *__restrict x, const TY *__restrict y) noexcept {
    while (count--) {
        (*x++) *= (*y++);
    }
}

//////////////////////////////////////////////////////////////////////////////

WorkerFFTW::WorkerFFTW(int channel_count, const std::set<int> &window_lengths)
        : max_window_length(*window_lengths.rbegin()),
          input_buffers(channel_count, max_window_length, fftw_alloc_real, fftw_free),
          output_buffers(channel_count, max_window_length / 2 + 1, fftw_alloc_complex, fftw_free) {
    for (int window_length : window_lengths) {
        fftw_plan plan = fftw_plan_dft_r2c_1d(
                window_length,
                input_buffers[0],
                output_buffers[0],
                FFTW_ESTIMATE | FFTW_DESTROY_INPUT
        );
        plans[window_length].reset(plan, fftw_destroy_plan);
    }
}

WorkerFFTW::WorkerFFTW(const WorkerFFTW &source)
        : max_window_length(source.max_window_length),
          input_buffers(source.input_buffers.height(), source.input_buffers.length(), fftw_alloc_real, fftw_free),
          output_buffers(source.output_buffers.height(), source.output_buffers.length(), fftw_alloc_complex, fftw_free),
          plans(source.plans) {}

void WorkerFFTW::compute(const SpectrogramRequest &request) {
    request.assertCorrectness();
    fftw_plan plan = getPlan(request.window_length);

    for (int c = 0; c < request.channel_count; ++c) {
        real *const input_buffer = input_buffers[c];
        std::fill(input_buffer + request.envelope_length, input_buffer + request.window_length, 0);
    }

    for (int h = 0; h < request.how_many; ++h) {
        index_t offset = request.input_offset + mulx(h, request.input_shift);
        IndexRange overlap = IndexRange(-offset, request.channel_length - offset).overlap(request.envelope_length);
        if (!overlap) {
            // all samples consist of zero-padding
            request.maxima[h] = ExtractedMaximum{0, 0};
            continue;
        }

        for (int c = 0; c < request.channel_count; ++c) {
            real *const input_buffer = input_buffers[c];
            const real *input = request.data[c] + offset;

            if (overlap.first_index > 0) {
                std::fill(input_buffer, input_buffer + overlap.first_index, 0);
            }
            std::copy(input + overlap.first_index, input + overlap.end_index, input_buffer + overlap.first_index);
            if (overlap.end_index < static_cast<index_t>(request.window_length)) {
                std::fill(input_buffer + overlap.end_index, input_buffer + request.window_length, 0);
            }

            mul(request.envelope_length, input_buffer, request.envelope);
            fftw_execute_dft_r2c(plan, input_buffer, output_buffers[c]);
        }
        request.maxima[h] = request.extractor(request.channel_count, request.output_bins,
                                              reinterpret_cast<complex *const *>(output_buffers.get()), request.correctors, input_buffers[0],
                                              nullptr);
    }
}

const complex *WorkerFFTW::computeSpectrum(Array1D<double> input, int window_length) {
    fftw_plan plan = getPlan(window_length);

    double *input_buffer = input_buffers[0];
    fftw_complex *output_buffer = output_buffers[0];

    std::copy(input.get(), input.get() + input.length(), input_buffer);
    if (input.length() < window_length) {
        std::fill(input_buffer + input.length(), input_buffer + window_length, 0);
    }

    fftw_execute_dft_r2c(plan, input_buffer, output_buffer);
    return reinterpret_cast<complex *>(output_buffer);
}

fftw_plan WorkerFFTW::getPlan(int window_length) const {
    auto it = plans.find(window_length);
    if (it == plans.end()) {
        throw std::logic_error("invalid window_length passed to FFTW calculator");
    }
    return it->second.get();
}

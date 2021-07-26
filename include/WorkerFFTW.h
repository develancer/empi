/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_WORKER_FFTW_H
#define EMPI_WORKER_FFTW_H

#include <map>
#include <memory>
#include <set>
#include <fftw3.h>
#include "Array.h"
#include "SpectrumCalculator.h"
#include "Worker.h"

/**
 * Short-time Fourier transform calculator based on FFTW.
 */
class WorkerFFTW : public Worker, public SpectrumCalculator {
    const int max_window_length;
    Array2D<real> input_buffers;
    Array2D<fftw_complex> output_buffers;
    std::map<int, std::shared_ptr<fftw_plan_s>> plans;

public:
    /**
     * Create a new FFTW worker for multi-channel signals with given number of channels,
     * able to serve pre-defined set of possible window lengths.
     *
     * @param channel_count number of channels in the signal
     * @param window_lengths list of window lengths that will be used with compute() calls later on
     */
    WorkerFFTW(int channel_count, const std::set<int> &window_lengths);

    /**
     * Create a new FFTW worker based on plan wisdom from an existing worker.
     * Both workers will be independent, won't share any internal buffers and will be able to run in parallel.
     *
     * @param source existing worker to take plan wisdom from
     */
    WorkerFFTW(const WorkerFFTW &source);

    /**
     * Compute a spectrogram of a given signal, according to the given specification.
     *
     * @param request specification for the spectrogram to be computed
     */
    void compute(const SpectrogramRequest &request) final;

    /**
     * Compute spectrum of a given signal with optional zero-padding.
     * Return the pointer to the internal buffer which may be overwritten on next calls to computeSpectrum.
     *
     * @param input input signal
     * @param window_length length of the window for FFT, should be at least the size of input
     * @return pointer to the internal buffer with size of at last window_length/2+1
     */
    const complex *computeSpectrum(Array1D<double> input, int window_length) final;

    void operator=(const WorkerFFTW &) = delete;

private:
    [[nodiscard]] fftw_plan getPlan(int window_length) const;
};

#endif //EMPI_WORKER_FFTW_H

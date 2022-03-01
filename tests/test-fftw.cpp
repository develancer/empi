/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <random>
#include "Array.h"
#include "Testing.h"
#include "SpectrogramCalculatorDummy.h"
#include "SpectrogramCalculatorFFTW.h"

double rand_float() {
    return (double) rand() / (double) RAND_MAX - 0.5;
}

void benchmark(SpectrogramCalculator &calculator, const SpectrogramRequest &request, const ExtractedMaximum *reference_maxima) {
    calculator.compute(request);

    if (reference_maxima) {
        double max_diff = 0;
        for (int i = 0; i < request.how_many; ++i) {
            ASSERT_EQUALS(reference_maxima[i].bin_index, request.maxima[i].bin_index);
            max_diff = std::max(max_diff, std::abs(request.maxima[i].energy - reference_maxima[i].energy));
        }
        ASSERT(max_diff > 0);
        ASSERT_NEAR_ZERO(max_diff);
    }
}

void run_test(Extractor extractor) {
    const int envelope_length = 800;
    Array1D<real> envelope(envelope_length);
    for (int i = 0; i < envelope_length; ++i) {
        envelope[i] = rand_float();
    }

    const int output_bins = 300;
    Array1D<Corrector> correctors(output_bins);
    for (int i = 0; i < output_bins; ++i) {
        correctors[i] = Corrector(0.5 * complex(rand_float(), rand_float()));
    }

    const index_t channel_length = 12000;
    const int channel_count = 10;
    Array2D<real> data(channel_count, channel_length);
    for (int c = 0; c < channel_count; ++c) {
        real *channel = data[c];
        for (index_t i = 0; i < channel_length; ++i) {
            channel[i] = rand_float();
        }
    }

    SpectrogramRequest request;
    request.data = data.get();
    request.channel_length = 12000;
    request.channel_count = 10;
    request.input_offset = -1000;
    request.input_shift = 200;
    request.how_many = 64;
    request.envelope = envelope.get();
    request.envelope_length = envelope_length;
    request.window_length = 1024;
    request.output_bins = output_bins;
    request.correctors = correctors.get();
    request.extractor = extractor;

    Array1D<ExtractedMaximum> maxima_dummy(request.how_many);
    Array1D<ExtractedMaximum> maxima_fftw(request.how_many);

    SpectrogramCalculatorFFTW fftw(channel_count, {1024});
    SpectrogramCalculatorDummy dummy;

    request.maxima = maxima_dummy.get();
    benchmark(dummy, request, nullptr);

    request.maxima = maxima_fftw.get();
    benchmark(fftw, request, maxima_dummy.get());
}

int main() {
    run_test(extractorVariablePhase);
    run_test(extractorConstantPhase);
    puts("OK");
}

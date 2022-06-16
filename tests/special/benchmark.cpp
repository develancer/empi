#include "Array.h"
#include "PinnedArray.h"
#include "SpectrogramCalculatorCUDA.h"
#include "SpectrogramCalculatorFFTW.h"
#include "Timer.h"

const int REPEATS = 10;

double rand_float() {
    return (real) rand() / (real) RAND_MAX - 0.5;
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr, "USAGE: %s total_input_length window_length\n", argv[0]);
        return 1;
    }
    srand(time(nullptr));
    const long total_input_length = atol(argv[1]);
    const int window_length = atoi(argv[2]);
    const int channel_count = 10;

    const int input_shift = window_length / 8;
    const int how_many = (total_input_length - window_length) / input_shift + 1;
    const int spectrum_length = window_length/2 + 1;

    PinnedArray1D<double> envelope(window_length);
    PinnedArray1D<Corrector> correctors(spectrum_length);
    PinnedArray1D<ExtractedMaximum> maximaCUDA(how_many);
    Array1D<ExtractedMaximum> maximaFFTW(how_many);

    SpectrogramRequest request;
    request.channel_length = total_input_length;
    request.channel_count = channel_count;
    request.input_offset = 0;
    request.input_shift = input_shift;
    request.how_many = how_many;
    request.envelope = envelope.get();
    request.envelope_length = window_length;
    request.window_length = window_length;
    request.output_bins = spectrum_length;
    request.correctors = correctors.get();
    request.extractor = extractorVariablePhase;

    PinnedArray2D<real> inputs(channel_count, total_input_length);
    request.data = inputs.get();

    WorkerFFTW fftw(channel_count, {window_length });
    WorkerCUDA cuda(channel_count, {request });

    double fftwTime = 0.0, cudaTime = 0.0;
    for (int repeat=0; repeat<REPEATS; ++repeat) {
        Timer fftwTimer, cudaTimer;
        for (int i = 0; i<window_length; ++i) {
            envelope[i] = rand_float();
        }
        for (int i = 0; i<spectrum_length; ++i) {
            correctors[i] = Corrector(rand_float());
        }
        for (int c=0; c<channel_count; ++c) {
            for (long i = 0; i<total_input_length; ++i) {
                inputs[c][i] = rand_float();
            }
        }

        request.maxima = maximaFFTW.get();
        fftwTimer.start();
        fftw.start(request);
        fftw.waitForCompletion();
        fftwTime += fftwTimer.time();
        fftwTimer.stop();

        request.maxima = maximaCUDA.get();
        cudaTimer.start();
        cuda.start(request);
        cuda.waitForCompletion();
        cudaTime += cudaTimer.time();
        cudaTimer.stop();

        double max_diff = 0;
        for (int i=0; i<how_many; ++i) {
            if (maximaFFTW[i].bin_index != maximaCUDA[i].bin_index) {
                printf("ERROR: MISMATCH at %d! %d != %d\n", i, maximaFFTW[i].bin_index, maximaCUDA[i].bin_index);
                return 1;
            }
            max_diff = std::max(max_diff, std::abs(maximaFFTW[i].energy - maximaCUDA[i].energy));
        }
        if (max_diff > 1.0e-10) {
            printf("ERROR: max diff = %le\n", max_diff);
            return 1;
        }
    }

    printf("%.6lf %.6lf\n", fftwTime/REPEATS, cudaTime/REPEATS);

    return 0;
}

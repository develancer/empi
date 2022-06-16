/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "Array.h"
#include "BlockAtom.h"
#include "BlockDictionary.h"
#include "BlockHelper.h"
#include "Extractor.h"
#include "GaussianFamily.h"
#include "SpectrogramCalculatorFFTW.h"
#include "SpectrogramRequest.h"
#include "Testing.h"

const int N = 300;

void test(SpectrogramCalculatorFFTW &calculator, double frequency, int position, double scale, double phase0, double phase1, double amplitude0, double amplitude1) {
    std::shared_ptr<GaussianFamily> family = std::make_shared<GaussianFamily>();
    index_t envelope_offset;
    index_t envelope_length = family->size_for_values(position, scale, nullptr);
    PinnedArray1D<double> envelope(envelope_length);
    family->generate_values(position, scale, &envelope_offset, envelope.get(), true);

    PinnedArray2D<double> data(2, N);
    const double value_max = envelope[envelope.length() / 2];

    double energies[2];
    energies[0] = energies[1] = 0.0;
    for (int i = 0; i < N; ++i) {
        // a single dual-channel Gabor atom with specified phases and amplitudes
        const int io = i - envelope_offset;
        const double v = (io >= 0 && io < envelope_length) ? envelope[io] / value_max : 0.0;
        const double phi = 2 * M_PI * frequency * (i - position);
        data[0][i] = amplitude0 * v * std::cos(phi + phase0);
        data[1][i] = amplitude1 * v * std::cos(phi + phase1);
        energies[0] += data[0][i] * data[0][i];
        energies[1] += data[1][i] * data[1][i];
    }
    auto converter = std::make_shared<BlockAtomParamsConverter>();
    auto correctors = BlockHelper::generate_correctors(envelope, 256, 129, calculator);
    BlockDictionary dictionary(Block(data, family, scale, envelope, correctors, converter, NAN, 256, 1, extractorVariablePhase));

    std::list<SpectrogramRequest> requests;
    dictionary.fetch_requests({0, N}, requests);
    for (const auto &request : requests) {
        calculator.compute(request);
        request.interface->notify();
    }

    BasicAtomPointer atom = dictionary.get_best_match();
    ASSERT(atom);

    const BlockAtom &block_atom = static_cast<BlockAtom &>(*atom);

    ASSERT_NEAR_ZERO(block_atom.params.scale - scale);
    ASSERT_NEAR_ZERO(block_atom.params.frequency - frequency);
    ASSERT_NEAR_ZERO(block_atom.params.position - position);
    ASSERT_NEAR_ZERO(block_atom.energy - (energies[0] + energies[1]));
}

int main() {
    SpectrogramCalculatorFFTW fftw(2, {256});

    test(fftw,
         9.0 / 256, 144, 10,
         0.71529, 0.92517,
         1.5, 2.5);
    test(fftw,
         119.0 / 256, 148, 10,
         0.71529, 0.92517,
         1.5, 2.5);
    test(fftw,
         0.25, 100, 1.0,
         0.0, 0.0,
         1.5, 2.5);

    puts("OK");
}

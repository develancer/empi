#include <cassert>
#include <cstdio>

#include "Array.h"
#include "BlockAtom.h"
#include "BlockDictionary.h"
#include "Extractor.h"
#include "Family.h"
#include "WorkerFFTW.h"
#include "SpectrogramRequest.h"

const int N = 300;

void test(WorkerFFTW& calculator, double frequency, int position, double scale, double phase0, double phase1, double amplitude0, double amplitude1)
{
    std::shared_ptr<GaussianFamily> family = std::make_shared<GaussianFamily>();
    index_t envelope_offset;
    index_t envelope_length = family->size_for_values(position, scale, nullptr);
    Array1D<double> envelope(envelope_length);
    family->generate_values(position, scale, &envelope_offset, envelope.get(), true);

    PinnedArray2D<double> data(2, N);
    const double value_max = envelope[envelope.length() / 2];

    double energies[2];
    energies[0] = energies[1] = 0.0;
    for (int i=0; i<N; ++i) {
        // a single dual-channel Gabor atom with specified phases and amplitudes
        const int io = i - envelope_offset;
        const double v = (io >= 0 && io < envelope_length) ? envelope[io] / value_max : 0.0;
        const double phi = 2 * M_PI * frequency * (i - position);
        data[0][i] = amplitude0 * v * std::cos(phi + phase0);
        data[1][i] = amplitude1 * v * std::cos(phi + phase1);
        energies[0] += data[0][i] * data[0][i];
        energies[1] += data[1][i] * data[1][i];
    }
    BlockDictionary dictionary(data, family);
    dictionary.add_block(scale, 256, 1, 129, extractorVariablePhase, calculator);

    std::list<SpectrogramRequest> requests;
    dictionary.fetch_requests({0, N}, requests);
    for (const auto& request : requests) {
        calculator.compute(request);
        request.interface->notify();
    }

    BasicAtomPointer atom = dictionary.get_best_match();
    assert(atom);

    const BlockAtom& block_atom = static_cast<BlockAtom&>(*atom);

    assert(std::abs(block_atom.scale - scale) < 1.0e-10);
    assert(std::abs(block_atom.frequency - frequency) < 1.0e-10);
    assert(std::abs(block_atom.position - position) < 1.0e-10);
    assert(std::abs(block_atom.get_energy() - (energies[0]+energies[1])) < 1.0e-10);
}

int main() {
    WorkerFFTW fftw(2, {256 });

    test(fftw,
         9.0/256, 144, 10,
         0.71529, 0.92517,
         1.5, 2.5);
    test(fftw,
         119.0/256, 148, 10,
         0.71529, 0.92517,
         1.5, 2.5);
    test(fftw,
         0.25, 100, 1.0,
         0.0, 0.0,
         1.5, 2.5);

    puts("OK");
}
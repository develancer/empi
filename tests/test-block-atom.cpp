#include <cassert>
#include <cstdio>

#include "Array.h"
#include "BlockAtom.h"
#include "BlockDictionary.h"
#include "Family.h"

void test(double frequency, double position, double scale, double phase0, double phase1, double amplitude0, double amplitude1)
{
    std::shared_ptr<GaussianFamily> family = std::make_shared<GaussianFamily>();
    index_t envelope_offset;
    index_t envelope_length = family->size_for_values(position, scale, nullptr);
    Array1D<double> envelope(envelope_length);
    family->generate_values(position, scale, &envelope_offset, envelope.get(), true);

    PinnedArray2D<double> data(2, 100);
    const double value_max = envelope[envelope.length() / 2];

    double energies[2];
    energies[0] = energies[1] = 0.0;
    for (int i=0; i<100; ++i) {
        // a single dual-channel Gabor atom with specified phases and amplitudes
        const int io = i - envelope_offset;
        const double v = (io >= 0 && io < envelope_length) ? envelope[io] / value_max : 0.0;
        const double phi = 2 * M_PI * frequency * (i - position);
        data[0][i] = amplitude0 * v * std::cos(phi + phase0);
        data[1][i] = amplitude1 * v * std::cos(phi + phase1);
        energies[0] += data[0][i] * data[0][i];
        energies[1] += data[1][i] * data[1][i];
    }

    BlockAtom block_atom(data, NAN, family, frequency, position, extractorVariablePhase);
    block_atom.scale = scale;

    std::shared_ptr<BlockExtendedAtom> extended = std::dynamic_pointer_cast<BlockExtendedAtom>(block_atom.extend());
    assert(extended);
    assert(extended->extra.length() == 2);

    assert(std::fmod(std::abs(extended->extra[0].phase - phase0), 2*M_PI) < 1.0e-10);
    assert(std::fmod(std::abs(extended->extra[1].phase - phase1), 2*M_PI) < 1.0e-10);

    assert(std::abs(extended->extra[0].amplitude - amplitude0) < 1.0e-10);
    assert(std::abs(extended->extra[1].amplitude - amplitude1) < 1.0e-10);

    assert(std::abs(extended->extra[0].energy - energies[0]) < 1.0e-10);
    assert(std::abs(extended->extra[1].energy - energies[1]) < 1.0e-10);

    BlockDictionary dictionary(data, family);
    dictionary.subtract_from_signal(*extended);

    for (int c=0; c<2; ++c) {
        for (int i=0; i<100; ++i) {
            assert(std::abs(data[c][i]) < 1.0e-10);
        }
    }
}

int main() {
    test(0.051, 40, 10,
         0.71529, 0.92517,
         1.5, 2.5);
    test(0.451, 50, 10,
         0.71529, 0.92517,
         1.5, 2.5);
    test(0.25, 60, 0.7,
         0.1, 0.2,
         1.3, 2.8);
    test(0.0001, 50, 5,
         0.0, M_PI,
         2.0, 3.0);
    test(0.0, 50, 5,
         0.0, M_PI,
         2.0, 3.0);
    test(0.5, 50, 5,
         0.0, M_PI,
         2.0, 3.0);

    puts("OK");
}

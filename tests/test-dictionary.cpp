#include <cassert>
#include <cstdio>
#include <list>

#include "Array.h"
#include "BlockAtom.h"
#include "BlockDictionary.h"
#include "Family.h"
#include "SpectrumCalculator.h"
#include "SpectrogramRequest.h"

class SpectrumCalculatorForTest : public SpectrumCalculator {
    Array1D<complex> result;

public:
    SpectrumCalculatorForTest() : result(9) {}

    const complex *computeSpectrum(Array1D<double> input, int window_length) final {
        assert(input.length() == 5);
        assert(window_length == 16);
        return result.get();
    }
};

int main() {
    SpectrumCalculatorForTest calculator;
    PinnedArray2D<double> data(1, 11);
    auto family = std::make_shared<GaussianFamily>();

    BlockDictionary dictionary(data, family);
    dictionary.add_block(1.0, 16, 1, 4, nullptr, calculator);

    std::list<SpectrogramRequest> requests;
    dictionary.fetch_requests({0, 11}, requests);

    assert(requests.size() == 1);

    const SpectrogramRequest &request = requests.front();
    assert(request.channel_length == 11);
    assert(request.channel_count == 1);
    assert(request.input_offset == -2);
    assert(request.input_shift == 1);
    assert(request.how_many == 11);
    assert(request.envelope_length == 5);
    assert(request.window_length == 16);
    assert(request.output_bins == 4);

    for (int i = 0; i < 11; ++i) {
        request.maxima[i].energy = std::min(i + 1, 11 - i);
        request.maxima[i].bin_index = i % 4;
    }
    request.interface->notify();

    BasicAtomPointer atom_pointer = dictionary.get_best_match();
    BlockAtom atom = *dynamic_cast<BlockAtom *>(atom_pointer.get());
    assert(atom.position == 5);

    puts("OK");
}

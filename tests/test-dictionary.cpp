/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <list>
#include "Array.h"
#include "BlockAtom.h"
#include "BlockDictionary.h"
#include "BlockHelper.h"
#include "GaussianFamily.h"
#include "SpectrumCalculator.h"
#include "SpectrogramRequest.h"
#include "Testing.h"

class SpectrumCalculatorForTest : public SpectrumCalculator {
    Array1D<complex> result;

public:
    SpectrumCalculatorForTest() : result(9) {}

    const complex *computeSpectrum(Array1D<double> input, int window_length) final {
        ASSERT_EQUALS(7, input.length());
        ASSERT_EQUALS(16, window_length);
        return result.get();
    }
};

int main() {
    SpectrumCalculatorForTest calculator;
    PinnedArray2D<double> data(1, 11);
    auto family = std::make_shared<GaussianFamily>();

    auto converter = std::make_shared<BlockAtomParamsConverter>();
    BlockDictionary dictionary(BlockHelper::create_block(data, family, 2.0, converter, NAN, 16, 4, 1, 0.0, extractorVariablePhase, calculator));

    std::list<SpectrogramRequest> requests;
    dictionary.fetch_requests({0, 11}, requests);

    ASSERT_EQUALS(1, requests.size());

    const SpectrogramRequest &request = requests.front();
    ASSERT_EQUALS(11, request.channel_length);
    ASSERT_EQUALS(1, request.channel_count);
    ASSERT_EQUALS(-3, request.input_offset);
    ASSERT_EQUALS(1, request.input_shift);
    ASSERT_EQUALS(11, request.how_many);
    ASSERT_EQUALS(7, request.envelope_length);
    ASSERT_EQUALS(16, request.window_length);
    ASSERT_EQUALS(4, request.output_bins);

    for (int i = 0; i < 11; ++i) {
        request.maxima[i].energy = std::min(i + 1, 11 - i);
        request.maxima[i].bin_index = i % 4;
    }
    request.interface->notify();

    BasicAtomPointer atom_pointer = dictionary.get_best_match();
    BlockAtom atom = *dynamic_cast<BlockAtom *>(atom_pointer.get());
    ASSERT_EQUALS(5, atom.params.position);

    puts("OK");
}

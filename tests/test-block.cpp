/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "Block.h"
#include "Corrector.h"
#include "Extractor.h"
#include "GaussianFamily.h"
#include "SpectrogramRequest.h"
#include "Testing.h"
#include "Types.h"

ExtractedMaximum extractorForTest(int, int, complex *const *, const Corrector *, double *, ExtraData *) {
    return ExtractedMaximum();
}

int main() {
    PinnedArray2D<real> data(10, 300);
    PinnedArray1D<real> envelope(129);
    PinnedArray1D<Corrector> correctors(100);

    auto converter = std::make_shared<BlockAtomParamsConverter>();
    Block block(data, std::make_shared<GaussianFamily>(), 10.0, envelope, correctors, converter, 256, 30, extractorForTest);

    SpectrogramRequest request = block.buildRequest(0, 300);
    ASSERT_EQUALS(data.get(), request.data);
    ASSERT_EQUALS(300, request.channel_length);
    ASSERT_EQUALS(10, request.channel_count);
    ASSERT_EQUALS(-64, request.input_offset); // half of envelope
    ASSERT_EQUALS(30, request.input_shift);
    ASSERT_EQUALS(10, request.how_many); // centered at: 0, 30, 60, 90, 120, 150, 180, 210, 240 and 270
    ASSERT_EQUALS(envelope.get(), request.envelope);
    ASSERT_EQUALS(129, request.envelope_length);
    ASSERT_EQUALS(256, request.window_length);
    ASSERT_EQUALS(100, request.output_bins);
    ASSERT_EQUALS(correctors.get(), request.correctors);
    ASSERT_EQUALS(extractorForTest, request.extractor);

    request = block.buildRequest(64, 300);
    ASSERT_EQUALS(-64, request.input_offset); // half of envelope
    ASSERT_EQUALS(10, request.how_many); // still the same

    request = block.buildRequest(65, 300);
    ASSERT_EQUALS(-34, request.input_offset); // just got outside the first envelope [-64..64]
    ASSERT_EQUALS(9, request.how_many); // one less

    request = block.buildRequest(65, 206);
    ASSERT_EQUALS(-34, request.input_offset);
    ASSERT_EQUALS(9, request.how_many); // still the same

    request = block.buildRequest(65, 205);
    ASSERT_EQUALS(-34, request.input_offset);
    ASSERT_EQUALS(8, request.how_many); // got outside the last envelope [206..334]

    puts("OK");
}

#include <cassert>
#include <cstdio>

#include "Array.h"
#include "Block.h"
#include "Corrector.h"
#include "Extractor.h"
#include "SpectrogramRequest.h"
#include "Types.h"

ExtractedMaximum extractorForTest(int, int, complex *const *, const Corrector *, double *, ExtraData *)
{
    return ExtractedMaximum();
}

int main() {
    PinnedArray2D<real> data(10, 300);
    PinnedArray1D<real> envelope(129);
    PinnedArray1D<Corrector> correctors(100);

    Block block(data, nullptr, envelope, correctors, 256, 30, extractorForTest);

    SpectrogramRequest request = block.buildRequest(0, 300);
    assert(request.data == data.get());
    assert(request.channel_length == 300);
    assert(request.channel_count == 10);
    assert(request.input_offset == -64); // half of envelope
    assert(request.input_shift == 30);
    assert(request.how_many == 10); // centered at: 0, 30, 60, 90, 120, 150, 180, 210, 240 and 270
    assert(request.envelope == envelope.get());
    assert(request.envelope_length == 129);
    assert(request.window_length == 256);
    assert(request.output_bins == 100);
    assert(request.correctors == correctors.get());
    assert(request.extractor == extractorForTest);

    request = block.buildRequest(64, 300);
    assert(request.input_offset == -64); // half of envelope
    assert(request.how_many == 10); // still the same

    request = block.buildRequest(65, 300);
    assert(request.input_offset == -34); // just got outside the first envelope [-64..64]
    assert(request.how_many == 9); // one less

    request = block.buildRequest(65, 206);
    assert(request.input_offset == -34);
    assert(request.how_many == 9); // still the same

    request = block.buildRequest(65, 205);
    assert(request.input_offset == -34);
    assert(request.how_many == 8); // got outside the last envelope [206..334]

    puts("OK");
}

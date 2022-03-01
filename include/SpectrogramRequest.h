/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_SPECTROGRAM_REQUEST_H
#define EMPI_SPECTROGRAM_REQUEST_H

#include <cassert>
#include <complex>
#include "BlockInterface.h"
#include "Corrector.h"
#include "Extractor.h"
#include "ProtoRequest.h"
#include "Types.h"

/**
 * Plain data object for a single request for spectrogram calculation.
 * Such request can be served by SpectrogramCalculator implementations by
 * 1) computing a windowed spectrogram of each channel in given multi-channel input data, and
 * 2) using a specified Extractor implementation to find the maximum value of each spectrum.
 */
struct SpectrogramRequest : public ProtoRequest {

    /**
     * pointer to the multi-channel data, of at least channel_count×channel_length samples
     * (in case of CUDA calculator, each channel MUST point to page-locked host memory)
     */
    const real *const *data;

    /**
     * length (in samples) of each channel (at least input_offset+envelope_length+(how_many-1)*input_shift)
     */
    index_t channel_length;

    /**
     * offset of the first sample of the input for the first transform
     */
    index_t input_offset;

    /**
     * envelope (window function) for calculating transforms (at least envelope_length items)
     * (in case of CUDA calculator, each channel MUST point to page-locked host memory)
     */
    const real *envelope;

    /**
     * objects capable of transforming FFT output into normalized energy values (at least output_bins items)
     * (in case of CUDA calculator, each channel MUST point to page-locked host memory)
     */
    const Corrector *correctors;

    /**
     * extractor implementation to reduce each spectrum to a single value
     */
    Extractor extractor;

    /**
     * buffer for extracted maximum values, one for each transform (space for at least how_many items)
     * (in case of CUDA calculator, this MUST point to page-locked host memory)
     */
    ExtractedMaximum *maxima;

    /**
     * listener to be notified after the request is complete
     */
    BlockInterface *interface;

    void assertCorrectness() const {
        assert(data);
        assert(channel_length > 0);
        assert(channel_count > 0);
        assert(input_shift > 0);
        assert(how_many > 0);

        assert(envelope);
        assert(envelope_length > 0);
        assert(window_length > 0);
        assert(envelope_length <= window_length);

        assert(output_bins > 0);
        assert(output_bins <= window_length / 2 + 1);
        assert(correctors);
        assert(extractor);
        assert(maxima);
    }
};

#endif //EMPI_SPECTROGRAM_REQUEST_H

/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_PROTO_REQUEST_H
#define EMPI_PROTO_REQUEST_H

/**
 * Plain data object consisting of the subset of fields from SpectrogramRequest
 * which can be known and obtained in advance from Dictionary instances.
 * This way, SpectrogramCalculator instances may pre-allocate all needed internal buffers.
 */
struct ProtoRequest {

    /** number of channels in input data */
    int channel_count;

    /** shift (in samples) between inputs for consecutive transforms */
    int input_shift;

    /** number of transforms to perform */
    int how_many;

    /** length (in samples) of the envelope */
    int envelope_length;

    /** number of samples for calculating each transform, can be larger than envelope_length (zero-padding) */
    int window_length;

    /** number of meaningful output bins from each transform (must not exceed window_length/2+1) */
    int output_bins;

};

#endif //EMPI_PROTO_REQUEST_H

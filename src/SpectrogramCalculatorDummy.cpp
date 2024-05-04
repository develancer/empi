/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <vector>
#include "Array.h"
#include "SpectrogramCalculatorDummy.h"

//////////////////////////////////////////////////////////////////////////////

void SpectrogramCalculatorDummy::compute(const SpectrogramRequest &request) {
    std::vector<double> tmp_for_extractor(request.output_bins);
    request.assertCorrectness();

    Array2D<complex> tmp_spectra(request.channel_count, request.output_bins);
    index_t input_offset = request.input_offset;
    for (int m = 0; m < request.how_many; ++m) {
        for (int c = 0; c < request.channel_count; ++c) {
            for (int k = 0; k < request.output_bins; ++k) {
                complex sum = 0.0;
                for (int n = 0; n < request.envelope_length; ++n) {
                    index_t index = input_offset + n;
                    if (index >= 0 && index < request.channel_length) {
                        sum += request.data[c][index] * request.envelope[n] *
                               std::polar(1.0, -2 * M_PI * k * n / request.window_length);
                    }
                }
                tmp_spectra[c][k] = sum;
            }
        }
        input_offset += request.input_shift;
        request.maxima[m] = request.extractor(request.channel_count, request.output_bins, tmp_spectra.get(),
                                              request.correctors, tmp_for_extractor.data(), nullptr);
    }
}

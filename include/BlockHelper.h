/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_HELPER_H
#define EMPI_BLOCK_HELPER_H

#include <memory>
#include "Array.h"
#include "Block.h"
#include "BlockAtomBase.h"
#include "Corrector.h"
#include "Family.h"
#include "PinnedArray.h"
#include "SpectrumCalculator.h"

class BlockHelper {
public:
    static Block create_block(PinnedArray2D<double> data, std::shared_ptr<Family> family, double scale,
                             std::shared_ptr<BlockAtomParamsConverter> converter, double booster,
                             int window_length, int output_bins, int input_shift, Extractor extractor, SpectrumCalculator& calculator, bool allow_overstep = true);

    static PinnedArray1D<double> generate_envelope(const Family* family, double scale);

    static PinnedArray1D<Corrector> generate_correctors(const Array1D<double>& envelope, int window_length, int output_bins, SpectrumCalculator& calculator);

    static std::map<double,int> compute_scales_and_transform_sizes(const Family* family, double scale_min, double scale_max, double log_scale_step, double df_scale);

    static int round_transform_size(double min_transform_size_as_float);
};

#endif //EMPI_BLOCK_HELPER_H

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_HELPER_H
#define EMPI_BLOCK_HELPER_H

#include <memory>
#include <utility>
#include "Array.h"
#include "Block.h"
#include "BlockAtomBase.h"
#include "BlockStructure.h"
#include "Corrector.h"
#include "Family.h"
#include "PinnedArray.h"
#include "SpectrumCalculator.h"

class BlockHelper {
    static void add_block_structures_to_list(std::list<BlockStructure>& result, const Family* family, double scale, double df_scale, double dt_scale);

public:
    static Block create_block(PinnedArray2D<double> data, std::shared_ptr<Family> family, double scale,
                             std::shared_ptr<BlockAtomParamsConverter> converter, double booster,
                             int window_length, int output_bins, int input_shift, double subsample_offset, Extractor extractor, SpectrumCalculator& calculator, bool allow_overstep = true);

    static std::pair<PinnedArray1D<double>, double> generate_envelope(const Family* family, double scale, double subsample_offset);

    static PinnedArray1D<Corrector> generate_correctors(const Array1D<double>& envelope, int window_length, int output_bins, SpectrumCalculator& calculator);

    static std::list<BlockStructure> compute_block_structures(const Family* family, double scale_min, double scale_max, double log_scale_step, double df_scale, double dt_scale);

    static int round_transform_size(double min_transform_size_as_float);
};

#endif //EMPI_BLOCK_HELPER_H

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "BlockDictionaryStructure.h"
#include "BlockHelper.h"

//////////////////////////////////////////////////////////////////////////////

BlockDictionaryStructure::BlockDictionaryStructure(std::shared_ptr<Family> family_, double energy_error, double scale_min, double scale_max, double frequency_max)
    : family(std::move(family_)), energy_error(energy_error), scale_min(scale_min), scale_max(scale_max), frequency_max(frequency_max),
    log_scale_step(family->inv_scale_integral(1 - energy_error)),
    dt_scale(family->inv_time_integral(1 - energy_error)),
    df_scale(family->inv_freq_integral(1 - energy_error)),
    scales_and_transform_sizes(BlockHelper::compute_scales_and_transform_sizes(family.get(), scale_min, scale_max, log_scale_step, df_scale))
{ }

[[nodiscard]] std::set<int> BlockDictionaryStructure::get_transform_sizes() const
{
    std::set<int> transform_sizes;
    for (const auto& pair : scales_and_transform_sizes) {
        transform_sizes.insert(transform_sizes.end(), pair.second);
    }
    return transform_sizes;
}
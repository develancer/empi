/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_DICTIONARY_STRUCTURE_H
#define EMPI_BLOCK_DICTIONARY_STRUCTURE_H

#include <memory>
#include <map>
#include <set>
#include "Family.h"

struct BlockDictionaryStructure {

    const std::shared_ptr<Family> family;
    const double energy_error, scale_min, scale_max, frequency_max, log_scale_step, dt_scale, df_scale;
    const std::map<double,int> scales_and_transform_sizes;

    /**
     * @param family Family object describing properties of the envelope function
     */
    BlockDictionaryStructure(std::shared_ptr<Family> family, double energy_error, double scale_min, double scale_max, double frequency_max);

    [[nodiscard]] std::set<int> get_transform_sizes() const;
};

#endif //EMPI_BLOCK_DICTIONARY_STRUCTURE_H

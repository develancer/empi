/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_ATOM_PRODUCT_CALCULATOR_H
#define EMPI_BLOCK_ATOM_PRODUCT_CALCULATOR_H

#include <memory>
#include "BlockAtomBase.h"
#include "Family.h"

class BlockAtomProductCalculator {
    std::shared_ptr<Family> family;

    [[nodiscard]] index_t calculate_mean_index(const BlockAtomParams &p, const BlockAtomParams &q, const IndexRange &pq_range) const;

public:
    BlockAtomProductCalculator(std::shared_ptr<Family> family);

    double calculate_product(const BlockAtomParams &p, const BlockAtomParams &q, double phase = 0.0) const;
};

#endif //EMPI_BLOCK_ATOM_PRODUCT_CALCULATOR_H

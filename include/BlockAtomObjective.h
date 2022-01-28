/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_ATOM_OBJECTIVE_H
#define EMPI_BLOCK_ATOM_OBJECTIVE_H

#include <array>
#include <memory>
#include <vector>
#include "Array.h"
#include "BlockAtomBase.h"
#include "Extractor.h"
#include "Family.h"

/**
 * Objective function used by multidimensional optimization to find optimal atom parameters.
 */
class BlockAtomObjective {

    std::shared_ptr<Family> family;
    const int channel_count;
    Array2D<double> data;
    Array2D<complex> products;
    Extractor extractor;
    std::vector<double> envelope;
    std::vector<complex> oscillating;
    std::shared_ptr<BlockAtomParamsConverter> converter;

public:
    BlockAtomObjective(std::shared_ptr<Family> family, Array2D<double> data, Extractor extractor, std::shared_ptr<BlockAtomParamsConverter> converter);

    double calculate_energy(const std::array<double, 3> &array, double *out_norm, ExtraData *out_extra_data);

    double calculate_energy(const BlockAtomParams &params, double *out_norm, ExtraData *out_extra_data);

    double operator()(const std::array<double, 3> &array) {
        return -calculate_energy(array, nullptr, nullptr);
    }
};

#endif //EMPI_BLOCK_ATOM_OBJECTIVE_H

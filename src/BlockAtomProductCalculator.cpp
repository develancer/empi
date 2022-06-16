/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "BlockAtomProductCalculator.h"
#include "Extractor.h"
#include "nelder_mead.h"

//////////////////////////////////////////////////////////////////////////////

BlockAtomProductCalculator::BlockAtomProductCalculator(std::shared_ptr<Family> family)
        : family(std::move(family)) {}

double BlockAtomProductCalculator::calculate_squared_product(const BlockAtomParams &p, const BlockAtomParams &q, double q_phase) const {
    IndexRange p_range = family->compute_range(p.position, p.scale);
    IndexRange q_range = family->compute_range(q.position, q.scale);

    IndexRange pq_range = p_range.overlap(q_range);

    double product = 0.0;
    if (pq_range.first_index < pq_range.end_index) {
        IndexRange wide_range(
                std::min(p_range.first_index, q_range.first_index),
                std::max(p_range.end_index, q_range.end_index)
        );
        complex sum_p2 = 0.0;
        double sum_ep2 = 0.0;
        double sum_q2 = 0.0;
        complex sum_pq = 0.0;
        for (index_t i = wide_range.first_index; i < wide_range.end_index; ++i) {
            complex p_value = 0.0;
            double q_value = 0.0;
            if (p_range.includes(i)) {
                double ep_value = family->value((i - p.position) / p.scale);
                p_value = std::polar(ep_value, 2 * M_PI * p.frequency * (i - p.position));
                sum_ep2 += ep_value * ep_value;
                sum_p2 += p_value * p_value;
            }
            if (q_range.includes(i)) {
                double eq_value = family->value((i - q.position) / q.scale);
                q_value = eq_value * cos(2 * M_PI * q.frequency * (i - q.position) + q_phase);
                sum_q2 += q_value * q_value;
            }
            sum_pq += p_value * q_value;
        }

        Corrector corrector(sum_p2 / sum_ep2);
        product = corrector.compute(sum_pq / std::sqrt(sum_q2 * sum_ep2));
    }
    return product;
}

double BlockAtomProductCalculator::calculate_squared_product(const BlockAtomParams &p, const BlockAtomParams &q) const {
    const int GRID = 10;
    double best_phase = 0.0;
    double min_product2 = INFINITY;
    for (int i=0; i<GRID; ++i) {
        double phase = 2 * M_PI * i / GRID;
        double product2 = calculate_squared_product(p, q, phase);
        if (product2 < min_product2) {
            min_product2 = product2;
            best_phase = phase;
        }
    }
    auto result = nelder_mead<double,1>(
            [&](const std::array<double,1>& x) {
                return this->calculate_squared_product(p, q, x[0]);
            },
            std::array<double,1>{best_phase},
            0.01,
            std::array<double,1>{0.1}
    );
    if (result.ifault) {
        throw std::runtime_error("cannot find optimal phase");
    }
    return result.ynewlo;
}

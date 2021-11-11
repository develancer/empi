/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "BlockAtomProductCalculator.h"

//////////////////////////////////////////////////////////////////////////////

BlockAtomProductCalculator::BlockAtomProductCalculator(std::shared_ptr<Family> family)
        : family(std::move(family)) {}

double BlockAtomProductCalculator::calculate_product(const BlockAtomParams &p, const BlockAtomParams &q, double p_phase) const {
    IndexRange p_range = family->compute_range(p.position, p.scale);
    IndexRange q_range = family->compute_range(q.position, q.scale);

    IndexRange pq_range = p_range.overlap(q_range);

    double product = 0.0;
    if (pq_range.first_index < pq_range.end_index) {
        IndexRange wide_range(
                std::min(p_range.first_index, q_range.first_index),
                std::max(p_range.end_index, q_range.end_index)
        );
        index_t mean_i = calculate_mean_index(p, q, pq_range);
        double p_phase_vs_mean_i = p_phase + 2 * M_PI * p.frequency * (mean_i - p.position);
        /*
        const double s02 = 1/( 1/(p.scale*p.scale) + 1/(q.scale*q.scale) );
        const double epsilon = exp(-4*M_PI*p.frequency*q.frequency*s02);
        double p_phase_vs_mean_i = p_phase + 2 * M_PI * p.frequency * (mean_i - p.position);
        double q_phase_vs_mean_i = atan2(
                sin(p_phase_vs_mean_i) * (1 - epsilon),
                cos(p_phase_vs_mean_i) * (1 + epsilon)
        ); */
        double sum_p2 = 0.0;
        double sum_q2 = 0.0;
        double sum_pq = 0.0;
        for (index_t i = wide_range.first_index; i < wide_range.end_index; ++i) {
            double p_value = p_range.includes(i)
                             ? family->value((i - p.position) / p.scale) * cos(2 * M_PI * p.frequency * (i - mean_i) + p_phase_vs_mean_i)
                             : 0.0;
            double q_value = q_range.includes(i)
                             ? family->value((i - q.position) / q.scale) * cos(2 * M_PI * q.frequency * (i - mean_i) + p_phase_vs_mean_i)
                             : 0.0;
            sum_p2 += p_value * p_value;
            sum_q2 += q_value * q_value;
            sum_pq += p_value * q_value;
        }
        product = sum_pq / std::sqrt(sum_p2 * sum_q2);
    }
    return product;
}

index_t BlockAtomProductCalculator::calculate_mean_index(const BlockAtomParams &p, const BlockAtomParams &q, const IndexRange &pq_range) const {
    // first, compute the product of the envelopes and find its mean position
    // so we know how to synchronize phases of atoms
    index_t ref_index = (pq_range.first_index + pq_range.end_index) / 2; // some arbitrary reference index

    double sum_envelopes = 0.0;
    double sum_envelopes_i = 0.0;
    for (index_t i = pq_range.first_index; i < pq_range.end_index; ++i) {
        double pq = family->value((i - p.position) / p.scale) * family->value((i - q.position) / q.scale);
        sum_envelopes += pq;
        sum_envelopes_i += (i - ref_index) * pq; // relative to ref_index, to minimize numerical errors
    }
    return Types::round<index_t>(sum_envelopes_i / sum_envelopes) + ref_index;
}

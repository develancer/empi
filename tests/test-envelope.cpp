/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "Array.h"
#include "GaussianFamily.h"
#include "TriangularFamily.h"
#include "Testing.h"

#define N 1000

void test_integrals(Family& family, const double delta_t = 0.3) {
    const double min_arg = family.min_arg();
    const double max_arg = family.max_arg();
    const double step = (max_arg - min_arg) / N;

    // any values would do fine
    const double delta_log_scale = 0.1;
    const double delta_f = 0.2;

    double sum2 = 0.0, sum2t = 0.0, sum2t2 = 0.0, sumA = 0.0, sumB = 0.0, sumC = 0.0;
    for (int i = 0; i < N; ++i) {
        double t = min_arg + (i + 0.5) * step;
        double value = family.value(t);
        sum2 += value * value;
        sum2t += value * value * t;
        sum2t2 += value * value * t * t;
        sumA += family.value(t * exp(delta_log_scale / 2)) * family.value(t * exp(-delta_log_scale / 2));
        sumB += value * value * std::cos(2 * M_PI * delta_f * t);
        sumC += family.value(t + delta_t / 2) * family.value(t - delta_t / 2);
    }
    ASSERT_APPROX(1.0, sum2 * step, 1.0e-5);
    ASSERT_NEAR_ZERO(sum2t * step);
    ASSERT_NEAR_ZERO(sum2t2 * step - 0.25 / M_PI);
    ASSERT_APPROX(sumA * step, family.scale_integral(delta_log_scale), 1.0e-5);
    ASSERT_APPROX(sumB * step, family.freq_integral(delta_f), 1.0e-5);
    ASSERT_APPROX(sumC * step, family.time_integral(delta_t), 1.0e-5);

    ASSERT_NEAR_ZERO(family.inv_scale_integral(family.scale_integral(delta_log_scale)) - delta_log_scale);
    ASSERT_NEAR_ZERO(family.inv_freq_integral(family.freq_integral(delta_f)) - delta_f);
    ASSERT_NEAR_ZERO(family.inv_time_integral(family.time_integral(delta_t)) - delta_t);
}

void test_generate(Family& family, int expected_sample_count, int expected_first_sample_offset) {
    index_t first_sample_offset;
    index_t sample_count = family.size_for_values(400.0, 100.0, &first_sample_offset);
    ASSERT_EQUALS(expected_sample_count, sample_count);
    ASSERT_EQUALS(expected_first_sample_offset, first_sample_offset);

    Array1D<double> values(sample_count);
    family.generate_values(400.0, 100.0, &first_sample_offset, values.get(), true);
    ASSERT_EQUALS(expected_first_sample_offset, first_sample_offset);

    double sum2 = 0.0, sum2t = 0.0, sum2t2 = 0.0;
    for (int i = 0; i < sample_count; ++i) {
        double t = static_cast<double>(first_sample_offset + i) - 400.0;
        double value = values[i];
        sum2 += value * value;
        sum2t += value * value * t;
        sum2t2 += value * value * t * t;
    }
    ASSERT_NEAR_ZERO(sum2 - 1.0);
    ASSERT_NEAR_ZERO(sum2t);
    ASSERT_APPROX(2500 / M_PI, sum2t2, 0.1);
}

int main() {
    GaussianFamily gauss;
    TriangularFamily triangle;

    test_integrals(gauss);
    test_integrals(triangle);
    test_integrals(triangle, 1.0);
    test_generate(gauss, 601, 100);
    test_generate(triangle, 179, 311);

    puts("OK");
}

#include <cassert>
#include <cstdio>

#include "Array.h"
#include "Family.h"

#define N 100

void test_integrals() {
    GaussianFamily gauss;
    const double min_arg = gauss.min_arg();
    const double max_arg = gauss.max_arg();
    const double step = (max_arg - min_arg) / N;

    // any values would do fine
    const double delta_log_scale = 0.1;
    const double delta_f = 0.2;
    const double delta_t = 0.3;

    double sum2 = 0.0, sum2t = 0.0, sum2t2 = 0.0, sumA = 0.0, sumB = 0.0, sumC = 0.0;
    for (int i=0; i<N; ++i) {
        double t = min_arg + (i+0.5)*step;
        double value = gauss.value(t);
        sum2 += value * value;
        sum2t += value * value * t;
        sum2t2 += value * value * t * t;
        sumA += gauss.value(t * exp(delta_log_scale/2)) * gauss.value(t * exp(-delta_log_scale/2));
        sumB += value * value * std::cos(2*M_PI*delta_f*t);
        sumC += gauss.value(t + delta_t/2) * gauss.value(t - delta_t/2);
    }
    assert(std::abs(sum2*step - 1.0) < 1.0e-12);
    assert(std::abs(sum2t*step) < 1.0e-12);
    assert(std::abs(sum2t2*step - 0.25/M_PI) < 1.0e-12);
    assert(std::abs(sumA*step - gauss.scale_integral(delta_log_scale)) < 1.0e-12);
    assert(std::abs(sumB*step - gauss.freq_integral(delta_f)) < 1.0e-12);
    assert(std::abs(sumC*step - gauss.time_integral(delta_t)) < 1.0e-12);

    assert(std::abs(gauss.inv_scale_integral(gauss.scale_integral(delta_log_scale)) - delta_log_scale) < 1.0e-12);
    assert(std::abs(gauss.inv_freq_integral(gauss.freq_integral(delta_f)) - delta_f) < 1.0e-12);
    assert(std::abs(gauss.inv_time_integral(gauss.time_integral(delta_t)) - delta_t) < 1.0e-12);
}

void test_generate() {
    GaussianFamily gauss;
    index_t first_sample_offset;
    index_t sample_count = gauss.size_for_values(40.0, 10.0, &first_sample_offset);
    assert(sample_count == 61);
    assert(first_sample_offset == 10);

    Array1D<double> values(sample_count);
    gauss.generate_values(40.0, 10.0, &first_sample_offset, values.get(), true);
    assert(first_sample_offset == 10);

    double sum2 = 0.0, sum2t = 0.0, sum2t2 = 0.0;
    for (int i=0; i<sample_count; ++i) {
        double t = static_cast<double>(first_sample_offset + i) - 40.0;
        double value = values[i];
        sum2 += value * value;
        sum2t += value * value * t;
        sum2t2 += value * value * t * t;
    }
    assert(std::abs(sum2 - 1.0) < 1.0e-12);
    assert(std::abs(sum2t) < 1.0e-12);
    assert(std::abs(sum2t2 - 25/M_PI) < 1.0e-12);
}

int main() {
    test_integrals();
    test_generate();

    puts("OK");
}

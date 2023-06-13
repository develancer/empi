/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <random>
#include "Array.h"
#include "BlockAtomObjective.h"
#include "Corrector.h"
#include "Extractor.h"
#include "GaussianFamily.h"
#include "Testing.h"

const int sample_count = 1024;

std::pair<Array2D<double>, double> generate_atom(std::shared_ptr<Family> family, const BlockAtomParams& params, double phase) {
    Array2D<double> atom(1, sample_count);
    atom.fill(0.0);

    index_t offset;
    if (family->size_for_values(params.position, params.scale, &offset) >= sample_count || offset < 0) {
        throw std::runtime_error("test scale is too large");
    }
    double* values = atom[0];
    // generate atom with unit energy
    family->generate_values(params.position, params.scale, nullptr, values+offset, false);

    double sum2 = 0.0;
    for (int i=0; i<sample_count; ++i) {
        values[i] *= std::cos(2*M_PI*params.frequency*(i-params.position) + phase);
        sum2 += values[i] * values[i];
    }
    double norm = 1.0 / std::sqrt(sum2);
    for (int i=0; i<sample_count; ++i) {
        values[i] *= norm;
    }
    double amplitude = family->value(0.0) * norm;
    return std::make_pair(atom, amplitude);
}

void generate_random(double* values) {
    std::generate_n(values, sample_count, std::bind(
            std::normal_distribution<double>(),
            std::default_random_engine()
    ));
}

void test_single_atom(std::shared_ptr<Family> family, double scale, double frequency, double phase) {
    const double position = 0.5 * (sample_count - 1);
    BlockAtomParams params(frequency, position, scale);
    auto [atom, atom_amplitude] = generate_atom(family, params, phase);

    Array2D<double> data(1, sample_count);
    std::copy_n(atom[0], sample_count, data[0]);

    auto converter = std::make_shared<BlockAtomParamsConverter>();
    BlockAtomObjective objective(family, data, extractorVariablePhase, converter);

    ExtraData extra_data;
    double actual_energy = objective.calculate_energy(converter->arrayFromParams(params),nullptr, &extra_data);
    double actual_amplitude = extra_data.amplitude * family->value(0.0) / std::sqrt(scale);

    ASSERT_APPROX(1.0, actual_energy, 1.0e-5);
    ASSERT_APPROX(1.0, extra_data.energy, 1.0e-5);
    ASSERT_APPROX(atom_amplitude, actual_amplitude, 1.0e-5);

    double raw_phase_diff = std::abs(phase - extra_data.phase);
    double correct_phase_diff = std::min(raw_phase_diff, 2*M_PI - raw_phase_diff);
    ASSERT_NEAR_ZERO(correct_phase_diff);
}

void test_random(std::shared_ptr<Family> family, double scale, double frequency) {
    Array2D<double> data(1, sample_count);
    generate_random(data[0]);

    const double position = 0.5 * (sample_count - 1);
    BlockAtomParams params(frequency, position, scale);

    auto converter = std::make_shared<BlockAtomParamsConverter>();
    BlockAtomObjective objective(family, data, extractorVariablePhase, converter);

    ExtraData extra_data;
    double actual_energy = objective.calculate_energy(converter->arrayFromParams(params),nullptr, &extra_data);
    double actual_amplitude = extra_data.amplitude * family->value(0.0) / std::sqrt(scale);

    auto [atom, atom_amplitude] = generate_atom(family, params, extra_data.phase);

    double product = 0.0;
    for (int i=0; i<sample_count; ++i) {
        product += atom[0][i] * data[0][i];
    }
    double expected_energy = product * product;
    double expected_amplitude = product * atom_amplitude;

    ASSERT_APPROX(expected_energy, actual_energy, 1.0e-5);
    ASSERT_APPROX(expected_energy, extra_data.energy, 1.0e-5);
    ASSERT_APPROX(expected_amplitude, actual_amplitude, 1.0e-5);

    // check residual
    for (int i=0; i<sample_count; ++i) {
        data[0][i] -= product * atom[0][i];
    }
    double residual_for_atom = objective.calculate_energy(converter->arrayFromParams(params),nullptr, nullptr);
    ASSERT_NEAR_ZERO(residual_for_atom);
}

double compute_diff(std::shared_ptr<Family> family, double scale, double frequency) {
    Array2D<double> data(1, sample_count);
    generate_random(data[0]);

    const double position = 0.5 * (sample_count - 1);
    BlockAtomParams params(frequency, position, scale);

    auto converter = std::make_shared<BlockAtomParamsConverter>();
    BlockAtomObjective objective(family, data, extractorVariablePhase, converter);

    ExtraData extra_data;
    double actual_energy = objective.calculate_energy(converter->arrayFromParams(params),nullptr, &extra_data);
    double actual_amplitude = extra_data.amplitude * family->value(0.0) / std::sqrt(scale);

    auto [atom, atom_amplitude] = generate_atom(family, params, extra_data.phase);

    double product = 0.0;
    for (int i=0; i<sample_count; ++i) {
        product += atom[0][i] * data[0][i];
    }
    double expected_energy = product * product;

    return std::abs(actual_energy - expected_energy);
}

void test_for_scale_and_frequency(std::shared_ptr<Family> family, double scale, double frequency) {
    if (frequency != 0.5) {
        test_single_atom(family, scale, frequency, 0.0);
        test_single_atom(family, scale, frequency, +M_PI);
        test_single_atom(family, scale, frequency, -M_PI);
    }
    if (frequency != 0.0 && frequency != 0.5) {
        test_single_atom(family, scale, frequency, +M_PI_4);
        test_single_atom(family, scale, frequency, -M_PI_4);
    }
    if (frequency != 0.0) {
        test_single_atom(family, scale, frequency, +M_PI_2);
        test_single_atom(family, scale, frequency, -M_PI_2);
    }
    test_random(family, scale, frequency);
}

void test_for_scale(std::shared_ptr<Family> family, double scale) {
    test_for_scale_and_frequency(family, scale, 0.0);
    test_for_scale_and_frequency(family, scale, 0.001);
    test_for_scale_and_frequency(family, scale, 0.1);
    test_for_scale_and_frequency(family, scale, 0.4);
    test_for_scale_and_frequency(family, scale, 0.499);
    test_for_scale_and_frequency(family, scale, 0.5);
}

int main() {
    auto gauss = std::make_shared<GaussianFamily>();
    test_for_scale(gauss, 10.0);
    test_for_scale(gauss, 50.0);
    puts("OK");
}

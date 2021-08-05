/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "Array.h"
#include "BlockAtom.h"
#include "BlockDictionary.h"
#include "Family.h"
#include "Testing.h"

double compute_residual(PinnedArray2D<double> &original_data, const BlockExtendedAtom &atom, const std::shared_ptr<Family> &family) {
    PinnedArray2D<double> data(original_data.height(), original_data.length());
    for (int c = 0; c < original_data.height(); ++c) {
        std::copy(original_data[c], original_data[c] + original_data.length(), data[c]);
    }
    BlockDictionary dictionary(data, family);
    dictionary.subtract_from_signal(atom);

    double residual_energy = 0.0;
    for (int c = 0; c < 2; ++c) {
        for (int i = 0; i < 100; ++i) {
            residual_energy += data[c][i] * data[c][i];
        }
    }
    return residual_energy;
}

void test_constant_phase(double frequency, double position, double scale, double phase0, double phase1, double amplitude0, double amplitude1) {
    std::shared_ptr<GaussianFamily> family = std::make_shared<GaussianFamily>();
    index_t envelope_offset;
    index_t envelope_length = family->size_for_values(position, scale, nullptr);
    Array1D<double> envelope(envelope_length);
    family->generate_values(position, scale, &envelope_offset, envelope.get(), true);

    PinnedArray2D<double> data(2, 100);
    const double value_max = envelope[envelope.length() / 2];

    double energies[2];
    energies[0] = energies[1] = 0.0;
    for (int i = 0; i < 100; ++i) {
        // a single dual-channel Gabor atom with specified phases and amplitudes
        const int io = i - envelope_offset;
        const double v = (io >= 0 && io < envelope_length) ? envelope[io] / value_max : 0.0;
        const double phi = 2 * M_PI * frequency * (i - position);
        data[0][i] = amplitude0 * v * std::cos(phi + phase0);
        data[1][i] = amplitude1 * v * std::cos(phi + phase1);
        energies[0] += data[0][i] * data[0][i];
        energies[1] += data[1][i] * data[1][i];
    }

    BlockAtom block_atom(data, NAN, family, frequency, position, extractorConstantPhase);
    block_atom.scale = scale;

    std::shared_ptr<BlockExtendedAtom> extended = std::dynamic_pointer_cast<BlockExtendedAtom>(block_atom.extend());
    ASSERT(extended);
    ASSERT_EQUALS(2, extended->extra.length());

    double original_residual = compute_residual(data, *extended, family);
    ASSERT_NEAR_ZERO(original_residual - 2.5);

    extended->scale += 0.01;
    ASSERT(compute_residual(data, *extended, family) > original_residual);
    extended->scale -= 0.02;
    ASSERT(compute_residual(data, *extended, family) > original_residual);
    extended->scale += 0.01;

    extended->frequency += 0.01;
    ASSERT(compute_residual(data, *extended, family) > original_residual);
    extended->frequency -= 0.02;
    ASSERT(compute_residual(data, *extended, family) > original_residual);
    extended->frequency += 0.01;

    extended->position += 1.0;
    ASSERT(compute_residual(data, *extended, family) > original_residual);
    extended->position -= 2.0;
    ASSERT(compute_residual(data, *extended, family) > original_residual);
    extended->position += 1.0;

    for (int c = 0; c < 2; ++c) {
        extended->extra[0].amplitude += 0.01;
        ASSERT(compute_residual(data, *extended, family) > original_residual);
        extended->extra[0].amplitude -= 0.02;
        ASSERT(compute_residual(data, *extended, family) > original_residual);
        extended->extra[0].amplitude += 0.01;
    }

    extended->extra[0].phase += 0.01;
    extended->extra[1].phase += 0.01;
    ASSERT(compute_residual(data, *extended, family) > original_residual);
    extended->extra[0].phase -= 0.02;
    extended->extra[1].phase -= 0.02;
    ASSERT(compute_residual(data, *extended, family) > original_residual);
    extended->extra[0].phase += 0.01;
    extended->extra[1].phase += 0.01;
}

void test_variable_phase(double frequency, double position, double scale, double phase0, double phase1, double amplitude0, double amplitude1) {
    std::shared_ptr<GaussianFamily> family = std::make_shared<GaussianFamily>();
    index_t envelope_offset;
    index_t envelope_length = family->size_for_values(position, scale, nullptr);
    Array1D<double> envelope(envelope_length);
    family->generate_values(position, scale, &envelope_offset, envelope.get(), true);

    PinnedArray2D<double> data(2, 100);
    const double value_max = envelope[envelope.length() / 2];

    double energies[2];
    energies[0] = energies[1] = 0.0;
    for (int i = 0; i < 100; ++i) {
        // a single dual-channel Gabor atom with specified phases and amplitudes
        const int io = i - envelope_offset;
        const double v = (io >= 0 && io < envelope_length) ? envelope[io] / value_max : 0.0;
        const double phi = 2 * M_PI * frequency * (i - position);
        data[0][i] = amplitude0 * v * std::cos(phi + phase0);
        data[1][i] = amplitude1 * v * std::cos(phi + phase1);
        energies[0] += data[0][i] * data[0][i];
        energies[1] += data[1][i] * data[1][i];
    }

    BlockAtom block_atom(data, NAN, family, frequency, position, extractorVariablePhase);
    block_atom.scale = scale;

    std::shared_ptr<BlockExtendedAtom> extended = std::dynamic_pointer_cast<BlockExtendedAtom>(block_atom.extend());
    ASSERT(extended);
    ASSERT_EQUALS(2, extended->extra.length());

    ASSERT_NEAR_ZERO(std::fmod(std::abs(extended->extra[0].phase - phase0), 2 * M_PI));
    ASSERT_NEAR_ZERO(std::fmod(std::abs(extended->extra[1].phase - phase1), 2 * M_PI));

    ASSERT_NEAR_ZERO(extended->extra[0].amplitude - amplitude0);
    ASSERT_NEAR_ZERO(extended->extra[1].amplitude - amplitude1);

    ASSERT_NEAR_ZERO(extended->extra[0].energy - energies[0]);
    ASSERT_NEAR_ZERO(extended->extra[1].energy - energies[1]);

    BlockDictionary dictionary(data, family);
    dictionary.subtract_from_signal(*extended);

    for (int c = 0; c < 2; ++c) {
        for (int i = 0; i < 100; ++i) {
            ASSERT_NEAR_ZERO(data[c][i]);
        }
    }
}

int main() {
    test_variable_phase(
            0.051, 40, 10,
            0.71529, 0.92517,
            1.5, 2.5);
    test_variable_phase(
            0.451, 50, 10,
            0.71529, 0.92517,
            1.5, 2.5);
    test_variable_phase(
            0.25, 60, 0.7,
            0.1, 0.2,
            1.3, 2.8);
    test_variable_phase(
            0.0001, 50, 5,
            0.0, M_PI,
            2.0, 3.0);
    test_variable_phase(
            0.0, 50, 5,
            0.0, M_PI,
            2.0, 3.0);
    test_variable_phase(
            0.5, 50, 5,
            0.0, M_PI,
            2.0, 3.0);
    test_constant_phase(
            0.125, 50, 10,
            0.0, 0.375 * M_PI,
            1.0, sqrt(M_SQRT2));
    puts("OK");
}

/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <err.h>
#include <memory>

#include "BlockAtomProductCalculator.h"
#include "GaussianFamily.h"

class OptimalityCalculator
{
    BlockAtomProductCalculator calculator;
    const double frequency_step;
    const double position_step;
    const double log_scale_step;

public:
    OptimalityCalculator(const std::shared_ptr<Family>& family, double energy_error)
    : calculator(family),
    frequency_step(0.5 * family->inv_freq_integral(1 - energy_error)),
    position_step(0.5 * family->inv_time_integral(1 - energy_error)),
    log_scale_step(0.5 * family->inv_scale_integral(1 - energy_error))
    { }

    double compute_optimality_factor(double scale, double frequency) {
        BlockAtomParams params{frequency, 0.0, scale};
        const double a = std::exp(log_scale_step);
        double result = 1.0;

        for (int i0 = -1; i0 <= 1; i0 += 2) {
            for (int i1 = -1; i1 <= 1; i1 += 2) {
                BlockAtomParams other = params;
                other.scale *= a;
                other.frequency += i0 * frequency_step / other.scale;
                if (other.frequency >= 0 && other.frequency <= 0.5) {
                    other.position += i1 * position_step * other.scale;
                    double y = calculator.calculate_squared_product(other, params);
                    result = std::min(result, y);
                }
            }
        }
        for (int i0 = -1; i0 <= 1; i0 += 2) {
            for (int i1 = -1; i1 <= 1; i1 += 2) {
                BlockAtomParams other = params;
                other.scale /= a;
                other.frequency += i0 * frequency_step / other.scale;
                if (other.frequency >= 0 && other.frequency <= 0.5) {
                    other.position += i1 * position_step * other.scale;
                    double y = calculator.calculate_squared_product(other, params);
                    result = std::min(result, y);
                }
            }
        }
        return result;
    }
};

void run_test() {
    // TODO
}

int main(int argc, char **argv) {
    if (argc == 1) {
        // actual test is run only if no parameters are passed
        run_test();
        puts("OK");
        return 0;
    }

    // otherwise, perform a numerical study for the optimality factor for the given envelope
    if (argc != 5) {
        errx(EXIT_FAILURE, "USAGE: %s family energy_error scale frequency", argv[0]);
    }
    srand(time(nullptr));
    auto it = Family::ALL.find(argv[1]);
    if (it == Family::ALL.end()) {
        errx(EXIT_FAILURE, "ERROR: invalid family %s", argv[1]);
    }

    const std::shared_ptr<Family> family = it->second;
    const double energy_error = atof(argv[2]);
    const double scale = atof(argv[3]);
    const double frequency = atof(argv[4]);

    if (frequency < 0 || frequency > 0.5) {
        errx(EXIT_FAILURE, "frequency %lf is invalid", frequency);
    }

    OptimalityCalculator calculator(family, energy_error);
    const double result = calculator.compute_optimality_factor(scale, frequency);
    printf("%lf %lf %lf  %lf\n", energy_error, scale, frequency, result);
}

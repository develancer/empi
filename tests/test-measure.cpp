/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <err.h>
#include <memory>

#include "BlockAtomProductCalculator.h"
#include "GaussianFamily.h"

int main(int argc, char **argv) {
    if (argc < 3) {
        // success as we don't want this test to fail, as it's not a real test
        errx(EXIT_SUCCESS, "USAGE: %s energy_error scale [ frequency phase ]", argv[0]);
    }
    srand(time(nullptr));
    const double energy_error = atof(argv[1]);
    const double scale = atof(argv[2]);

    auto family = std::make_shared<GaussianFamily>();
    const double frequency_step = 0.5 * family->inv_freq_integral(1 - energy_error) / scale;
    const double position_step = 0.5 * family->inv_time_integral(1 - energy_error) * scale;
    const double log_scale_step = 0.5 * family->inv_scale_integral(1 - energy_error);
    const double frequency = (argc > 3) ? atof(argv[3]) : 0.0;
    const double phase = (argc > 4) ? atof(argv[4]) : 0.0;
    BlockAtomParams params{frequency, 0.0, scale};

    if (frequency < 0 || frequency > 0.5) {
        errx(EXIT_FAILURE, "frequency %lf is invalid", frequency);
    }

    BlockAtomProductCalculator calculator(family);

    const double a = std::exp(log_scale_step);
    double result = 1.0;

    for (int i0 = -1; i0 <= 1; i0 += 2) {
        for (int i1 = -1; i1 <= 1; i1 += 2) {
            BlockAtomParams other = params;
            other.frequency += i0 * frequency_step / a;
            if (other.frequency >= 0 && other.frequency <= 0.5) {
                other.position += i1 * position_step * a;
                other.scale *= a;
                double y = calculator.calculate_squared_product(other, params, phase);
                result = std::min(result, y);
            }
        }
    }
    for (int i0 = -1; i0 <= 1; i0 += 2) {
        for (int i1 = -1; i1 <= 1; i1 += 2) {
            BlockAtomParams other = params;
            other.frequency += i0 * frequency_step * a;
            if (other.frequency >= 0 && other.frequency <= 0.5) {
                other.position += i1 * position_step / a;
                other.scale /= a;
                double y = calculator.calculate_squared_product(other, params, phase);
                result = std::min(result, y);
            }
        }
    }

    printf("%lf %lf %lf  %lf\n", energy_error, scale, frequency, result);
}

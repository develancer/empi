/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "BlockAtomProductCalculator.h"
#include "GaussianFamily.h"
#include "Testing.h"

void test_product() {
    auto family = std::make_shared<GaussianFamily>();
    BlockAtomProductCalculator calculator(family);

    BlockAtomParams x0{0.2, 2.0, 20.0};
    ASSERT_APPROX(1.0, calculator.calculate_squared_product(x0, x0), 1.0e-10);

    const double energy_error = 0.01;
    const double dl = family->inv_scale_integral(1 - energy_error);
    const double dt = family->inv_time_integral(1 - energy_error) * x0.scale;
    const double df = family->inv_freq_integral(1 - energy_error) / x0.scale;

    BlockAtomParams xf(x0), xt(x0), xs(x0);
    xf.frequency += df;
    xt.position += dt;
    xs.scale *= std::exp(dl);

    const double expected = std::pow(1 - energy_error, 2);
    ASSERT_APPROX(expected, calculator.calculate_squared_product(x0, xf), 1.0e-10);
    ASSERT_APPROX(expected, calculator.calculate_squared_product(x0, xt), 1.0e-10);
    ASSERT_APPROX(expected, calculator.calculate_squared_product(x0, xs), 1.0e-10);
}

int main() {
    test_product();

    puts("OK");
}

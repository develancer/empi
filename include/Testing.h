/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_TESTING_H
#define EMPI_TESTING_H

#include <cmath>
#include <err.h>

#define ASSERT(ACTUAL) if (!(ACTUAL)) errx(EXIT_FAILURE, #ACTUAL " evaluates to false")
#define ASSERT_EQUALS(EXPECTED, ACTUAL) if ((EXPECTED) != (ACTUAL)) errx(EXIT_FAILURE, #ACTUAL " does not equal expected " #EXPECTED)
#define ASSERT_NEAR_ZERO(ACTUAL) if (std::abs(ACTUAL) > 1.0e-9) errx(EXIT_FAILURE, #ACTUAL " (%le) is significantly nonzero", (ACTUAL))
#define ASSERT_APPROX(EXPECTED, ACTUAL, EPSILON) if (std::abs((EXPECTED)-(ACTUAL)) > (EPSILON)) \
  errx(EXIT_FAILURE, #ACTUAL " (%lf) is significantly (" #EPSILON ") different than " #EXPECTED " (%lf)", (ACTUAL), (EXPECTED))
#define ASSERT_SAME_PHASE(EXPECTED, ACTUAL, EPSILON) if (std::min(std::abs((EXPECTED)-(ACTUAL)), 2*M_PI-std::abs((EXPECTED)-(ACTUAL))) > (EPSILON)) \
  errx(EXIT_FAILURE, #ACTUAL " (%lf) and " #EXPECTED " (%lf) differ significantly (" #EPSILON ") in phase", (ACTUAL), (EXPECTED))

#endif //EMPI_TESTING_H

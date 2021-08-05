/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_TESTING_H
#define EMPI_TESTING_H

#include <cmath>
#include <stdexcept>

#define ASSERT(ACTUAL) if (!(ACTUAL)) throw std::runtime_error(#ACTUAL " evaluates to false");
#define ASSERT_EQUALS(EXPECTED, ACTUAL) if ((EXPECTED) != (ACTUAL)) throw std::runtime_error(#ACTUAL " does not equal expected " #EXPECTED)
#define ASSERT_NEAR_ZERO(ACTUAL) if (std::abs(ACTUAL) > 1.0e-10) throw std::runtime_error(#ACTUAL " is significantly nonzero")

#endif //EMPI_TESTING_H

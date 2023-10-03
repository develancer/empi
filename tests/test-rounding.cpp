/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "Testing.h"
#include "Types.h"

template<typename T>
void assert_equals(double input, T expected, T (*fun)(double)) {
    T actual = fun(input);
    ASSERT_EQUALS(expected, actual);
}

template<typename T>
void test_all_cases() {
    assert_equals<T>(-4.0, -4, Types::round);
    assert_equals<T>(-3.3, -3, Types::round);
    assert_equals<T>(-2.6, -3, Types::round);
    assert_equals<T>(-1.9, -2, Types::round);
    assert_equals<T>(-1.2, -1, Types::round);
    assert_equals<T>(-0.5, 0, Types::round);
    assert_equals<T>(0.2, 0, Types::round);
    assert_equals<T>(0.9, 1, Types::round);
    assert_equals<T>(1.6, 2, Types::round);
    assert_equals<T>(2.3, 2, Types::round);
    assert_equals<T>(3.0, 3, Types::round);

    assert_equals<T>(-4.0, -4, Types::ceil);
    assert_equals<T>(-3.3, -3, Types::ceil);
    assert_equals<T>(-2.6, -2, Types::ceil);
    assert_equals<T>(-1.9, -1, Types::ceil);
    assert_equals<T>(-1.2, -1, Types::ceil);
    assert_equals<T>(-0.5, 0, Types::ceil);
    assert_equals<T>(0.2, 1, Types::ceil);
    assert_equals<T>(0.9, 1, Types::ceil);
    assert_equals<T>(1.6, 2, Types::ceil);
    assert_equals<T>(2.3, 3, Types::ceil);
    assert_equals<T>(3.0, 3, Types::ceil);

    assert_equals<T>(-4.0, -4, Types::floor);
    assert_equals<T>(-3.3, -4, Types::floor);
    assert_equals<T>(-2.6, -3, Types::floor);
    assert_equals<T>(-1.9, -2, Types::floor);
    assert_equals<T>(-1.2, -2, Types::floor);
    assert_equals<T>(-0.5, -1, Types::floor);
    assert_equals<T>(0.2, 0, Types::floor);
    assert_equals<T>(0.9, 0, Types::floor);
    assert_equals<T>(1.6, 1, Types::floor);
    assert_equals<T>(2.3, 2, Types::floor);
    assert_equals<T>(3.0, 3, Types::floor);
}

int main() {
    test_all_cases<int>();
    test_all_cases<long>();
    test_all_cases<long long>();

    puts("OK");
}

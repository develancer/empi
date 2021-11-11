/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_MINIMIZER_H
#define EMPI_MINIMIZER_H

#include <array>
#include <atomic>
#include <functional>
#include <memory>
#include <gsl/gsl_multimin.h>

#define GSL_UNIQUE_POINTER(GSL_TYPE) std::unique_ptr<GSL_TYPE, decltype(&GSL_TYPE ## _free)>

template<int N>
class Minimizer {
    using minimizer_ptr = GSL_UNIQUE_POINTER(gsl_multimin_fminimizer);
    using vector_ptr = GSL_UNIQUE_POINTER(gsl_vector);

    minimizer_ptr minimizer;

    template<typename T>
    static double internal_target(const gsl_vector *x, void *p) {
        T *target = reinterpret_cast<T *>(p);
        std::array<double, N> array;
        for (int i = 0; i < N; ++i) {
            array[i] = gsl_vector_get(x, i);
        }
        return (*target)(array);
    }

public:
    double initial_step = 0.2;
    double tolerance = 1.0e-3;
    int max_iterations = 1000;

    Minimizer() : minimizer(
            gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2rand, N),
            &gsl_multimin_fminimizer_free
    ) {}

    template<typename T>
    int minimize(T *target, std::array<double, N> &array, const std::atomic<bool> *running = nullptr) {
        vector_ptr x(gsl_vector_alloc(N), &gsl_vector_free);
        for (int i = 0; i < N; ++i) {
            gsl_vector_set(x.get(), i, array[i]);
        }

        vector_ptr s(gsl_vector_alloc(N), &gsl_vector_free);
        gsl_vector_set_all(s.get(), initial_step);

        gsl_multimin_function fun;
        fun.n = N;
        fun.f = &internal_target<T>;
        fun.params = reinterpret_cast<void *>(target);
        gsl_multimin_fminimizer_set(minimizer.get(), &fun, x.get(), s.get());

        int iter = 0;
        int status;
        do {
            iter++;
            status = gsl_multimin_fminimizer_iterate(minimizer.get());
            if (status) break;

            status = gsl_multimin_test_size(gsl_multimin_fminimizer_size(minimizer.get()), tolerance);
            if (status == GSL_SUCCESS) {
                for (int i = 0; i < N; ++i) {
                    array[i] = gsl_vector_get(minimizer->x, i);
                }
                return iter;
            }
        } while ((!running || *running) && status == GSL_CONTINUE && iter < max_iterations);

        return 0;
    }
};

#undef GSL_UNIQUE_POINTER

#endif //EMPI_MINIMIZER_H

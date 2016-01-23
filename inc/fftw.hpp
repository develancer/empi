/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_FFTW_HPP
#define	EMPI_FFTW_HPP

#include <complex>
#include <cstring>
#include <fftw3.h>

/**
 * RAII-style wrapper for arrays allocated for FFTW.
 * @param T  type of values stored in array
 */
template<typename T>
class fftwArray {
	const int length;
	T* pointer = nullptr;

	fftwArray(const fftwArray&) =delete;

public:
	/**
	 * Allocate array of given length.
	 * The contents of the array will remain uninitialized.
     * @param length  array's requested length (number of elements)
     */
	fftwArray(int length)
	: length(length) {
		pointer = static_cast<T*>(fftw_malloc(sizeof(T)*length));
		if (!pointer) {
			throw std::bad_alloc();
		}
	}

	/**
	 * Destroy array by calling fftw_free on internal pointer.
     */
	~fftwArray() {
		if (pointer) fftw_free(pointer);
	}

	/**
	 * Return an internal pointer to the underlying data.
     * @return  C-style pointer to data
     */
	inline T* operator&() const {
		return pointer;
	}

	/**
	 * Return i-th element of an array (counting from i=0).
	 * No range-checking is performed.
     * @param i  index of requested element
     * @return reference to i-th value
     */
	inline T& operator[](int i) {
		return pointer[i];
	}

	/**
	 * Return i-th element of an array (counting from i=0).
	 * No range-checking is performed.
     * @param i  index of requested element
     * @return constant reference to i-th value
     */
	inline const T& operator[](int i) const {
		return pointer[i];
	}

	/**
     * @return array's length (number of elements)
     */
	inline int size() const {
		return length;
	}

	/**
	 * Fill contents of the array with zeroes.
     */
	inline void zero() {
		memset(pointer, 0, sizeof(T)*length);
	}
};

// convenience names for arrays allocated for FFTW
typedef fftwArray<double> fftwDouble;
typedef fftwArray<std::complex<double>> fftwComplex;

/**
 * RAII-style wrapper for FFTW plan structure.
 */
class fftwPlan {
	fftw_plan plan = nullptr;

	fftwPlan(const fftwPlan&) =delete;
	void operator=(const fftwPlan&) =delete;

public:
	/**
	 * Create plan of given length for complex-to-complex transform.
	 * All parameters are analogous to fftw_plan_dft_1d function.
	 */
	inline fftwPlan(int Nfft, fftwComplex& input, fftwComplex& output, int sign, unsigned flags) {
		#ifdef _OPENMP
		#pragma omp critical
		#endif
		{
			plan = fftw_plan_dft_1d(Nfft, reinterpret_cast<fftw_complex*>(&input), reinterpret_cast<fftw_complex*>(&output), sign, flags);
		}
	}

	/**
	 * Create plan of given length for real-to-complex transform.
	 * All parameters are analogous to fftw_plan_dft_r2c_1d specification.
	 */
	inline fftwPlan(int Nfft, fftwDouble& input, fftwComplex& output, unsigned flags) {
		#ifdef _OPENMP
		#pragma omp critical
		#endif
		{
			plan = fftw_plan_dft_r2c_1d(Nfft, &input, reinterpret_cast<fftw_complex*>(&output), flags);
		}
	}

	/**
	 * Destroy plan by calling fftw_destroy_plan.
     */
	inline ~fftwPlan() {
		#ifdef _OPENMP
		#pragma omp critical
		#endif
		{
			if (plan) fftw_destroy_plan(plan);
		}
	}

	/**
     * @return internal C-style type for plan
     */
	inline fftw_plan operator&() const {
		return plan;
	}

	/**
	 * Execute plan.
     */
	inline void execute() const {
		fftw_execute(plan);
	}
};

/**
 * Return smallest integer power of 2 greater or equal than given value.
 * @param x  numeric value of any type
 * @return  smallest 2^n >= x in the same type as x
 */
template<typename T>
inline T fftwRound(T x) {
	T y = 2;
	while (y < x) {
		y *= 2;
	}
	return y;
}

#endif	/* EMPI_FFTW_HPP */

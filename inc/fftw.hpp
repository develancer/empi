/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#ifndef EMPI_FFTW_HPP
#define	EMPI_FFTW_HPP

#include <complex>
#include <cstring>
#include <fftw3.h>

template<typename T>
class fftwArray {
	const long length;
	T* pointer = 0;

	fftwArray(const fftwArray&); // forbidden

public:
	fftwArray(size_t length)
	: length(length) {
		pointer = static_cast<T*>(fftw_malloc(sizeof(T)*length));
		if (!pointer) {
			throw std::bad_alloc();
		}
	}

	~fftwArray() {
		if (pointer) fftw_free(pointer);
	}

	T* operator&() const {
		return pointer;
	}

	T& operator[](size_t i) {
		return pointer[i];
	}

	const T& operator[](size_t i) const {
		return pointer[i];
	}

	fftwArray& operator=(const fftwArray& source) {
		memcpy(pointer, source.pointer, sizeof(T)*length);
		return *this;
	}

	size_t size() const {
		return length;
	}

	void zero() {
		memset(pointer, 0, sizeof(T)*length);
	}
};

class fftwPlan {
	fftw_plan plan;

	fftwPlan(const fftwPlan&); // forbidden
	void operator=(const fftwPlan&); // forbidden

public:
	inline fftwPlan(fftw_plan plan) : plan(plan) { }

	inline ~fftwPlan() {
		fftw_destroy_plan(plan);
	}

	inline fftw_plan operator&() const {
		return plan;
	}

	inline void execute() const {
		fftw_execute(plan);
	}
};

template<typename T>
inline T fftwRound(T x) {
	T y = 2;
	while (y < x) {
		y *= 2;
	}
	return y;
}

typedef fftwArray<double> fftwDouble;
typedef fftwArray<std::complex<double>> fftwComplex;

#endif	/* EMPI_FFTW_HPP */

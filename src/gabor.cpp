/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#include "fftw.hpp"
#include "gabor.hpp"

/**
 * If the scalar product between two normalized Gabor atoms is smaller
 * than this value, the atoms are classified as orthogonal.
 */
const double ORTHOGONALITY = 1.0e-8;

/**
 * Multiply given array/vector/buffer by discretely sampled Gabor function,
 * i.e. transform x[t] := x[t] g[t], where g[t] ~ exp(−π(t-t₀)²/s²).
 * The function is L²-normalized: Σ g[t]² Δt = 1.
 *
 * @param buffer  array values x[0], x[Δt], x[2Δt], ...
 * @param step    sampling step Δt (inverse of sampling frequency)
 * @param center  center t₀ of Gabor function
 * @param width   width s of Gabor function
 */
template<typename T>
static void gaussize(T& buffer, double step, double center, double width) {
	const double norm = sqrt(M_SQRT2 / width);
	size_t N = buffer.size();
	for (size_t i=0; i<N; ++i) {
		double x_width = (i * step - center) / width;
		buffer[i] *= norm * exp(-M_PI * x_width * x_width);
	}
}

//------------------------------------------------------------------------------

void GaborMapGenerator0::generate(TimeFreqMap<complex>& map, const SingleSignal& signal) {
	for (size_t ti=0; ti<map.tCount; ++ti) {
		map.t(ti) = ti * tMax / (map.tCount - 1);
	}
	for (size_t fi=0; fi<map.fCount; ++fi) {
		map.f(fi) = (fi + 1) * freqSampling / Nfft;
	}

	const size_t N = signal.samples.size();
	const double norm = sqrt(M_SQRT2 / map.s) / freqSampling;
	for (size_t ti=0; ti<map.tCount; ++ti) {
		for (size_t fi=0; fi<map.fCount; ++fi) {
			complex sum = 0.0;
			for (size_t i=0; i<N; ++i) {
				double x = i/freqSampling - map.t(ti);
				double x_width = x / map.s;
				sum += signal.samples[i] * exp(-M_PI * x_width * x_width) * std::polar(1.0, 2*M_PI*map.f(fi)*x);
			}
			map.value(fi, ti) = norm * sum;
		}
	}
}

void GaborMapGenerator1::generate(TimeFreqMap<complex>& map, const SingleSignal& signal) {
	for (size_t ti=0; ti<map.tCount; ++ti) {
		map.t(ti) = ti * tMax / (map.tCount - 1);
	}
	for (size_t fi=0; fi<map.fCount; ++fi) {
		map.f(fi) = (fi + 1) * freqSampling / Nfft;
	}

	fftwDouble input(Nfft);
	fftwComplex output(map.fCount + 1);
	fftwPlan plan(fftw_plan_dft_r2c_1d(Nfft, &input, reinterpret_cast<fftw_complex*>(&output), FFTW_ESTIMATE | FFTW_DESTROY_INPUT));

	const long N = signal.samples.size();
	for (size_t ti=0; ti<map.tCount; ++ti) {
		const double t0 = map.t(ti);
		long iL = std::max(0L, lrint((t0 - hwGabor) * freqSampling));
		long iR = std::min(N-1, lrint((t0 + hwGabor) * freqSampling));
		const double t0fixed = t0 - iL / freqSampling;

		input.zero();
		for (long i=iL; i<=iR; ++i) {
			input[i-iL] = signal.samples[i];
		}
		gaussize(input, 1.0/freqSampling, t0fixed, map.s);
		plan.execute();

		// kopiujemy wartości output[1..fCount]
		const double norm = 1.0 / freqSampling;
		const double mult = 2 * M_PI * t0fixed;
		for (size_t fi=0; fi<map.fCount; ++fi) {
			map.value(fi, ti) = norm * output[fi+1] * std::polar(1.0, map.f(fi) * mult);
		}
	}
}

void GaborMapGenerator2::generate(TimeFreqMap<complex>& map, const SingleSignal& signal) {
	for (size_t ti=0; ti<map.tCount; ++ti) {
		map.t(ti) = ti / freqSampling;
	}
	for (size_t fi=0; fi<map.fCount; ++fi) {
		map.f(fi) = (fi + 1) * 0.5*freqSampling / map.fCount;
	}

	fftwComplex transform(Nfft);
	fftwComplex input(Nfft);
	fftwComplex output(Nfft);
	fftwPlan planForward(fftw_plan_dft_1d(Nfft, reinterpret_cast<fftw_complex*>(&input), reinterpret_cast<fftw_complex*>(&transform), FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT));
	fftwPlan planBackward(fftw_plan_dft_1d(Nfft, reinterpret_cast<fftw_complex*>(&input), reinterpret_cast<fftw_complex*>(&output), FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT));

	input.zero();
	for (size_t i=0; i<signal.samples.size(); ++i) {
		input[i] = signal.samples[i];
	}
	planForward.execute();

	for (size_t fi=0; fi<map.fCount; ++fi) {
		const double f0 = map.f(fi);

		input = transform;
		gaussize(input, freqSampling/Nfft, f0, 1/map.s);
		planBackward.execute();

		const double norm = 1.0 / Nfft;
		for (size_t ti=0; ti<map.tCount; ++ti) {
			map.value(fi, ti) = norm * output[ti];
		}
	}
}

//------------------------------------------------------------------------------

struct GaborProductBase {
	const double s1, s12, s2, s22, S, A, s3, s32, product;

	GaborProductBase(double s1, double s2) :
	s1(s1), s12(s1*s1),
	s2(s2), s22(s2*s2),
	S(s12*s22), A(s12 + s22),
	s3(sqrt(S / A)), s32(s3*s3),
	product(M_SQRT2 * s3 / sqrt(s1 * s2))
	{ }
};

struct GaborProductCalculator {
	const GaborProductBase* const base;
	const double t1, t2, B, C, product;

	GaborProductCalculator(const GaborProductBase* base, double t1, double t2) :
	base(base),
	t1(t1), t2(t2),
	B(base->s12*t2 + base->s22*t1),
	C(base->s12*t2*t2 + base->s22*t1*t1),
	product(base->product * exp(-M_PI/base->S*(C-B*B/base->A)))
	{ }

	complex calculate(double f1, double f2) const {
		double f3 = f1 - f2;
		double phi3 = (t2 - t1) * (f1*base->s12 + f2*base->s22) / base->A;
		return product * exp(-M_PI*f3*f3*base->s32) * std::polar(1.0, 2*M_PI*phi3);
	}
};

Atom GaborWorkspace::findBestMatch() const {
	double modulusTotal = 0;
	Atom atomTotal;
	for (auto map : maps) {
		const double s = map->s;
		const double amplFactor = 2.0 * sqrt(M_SQRT2 / s);
		#ifdef _OPENMP
		#pragma omp parallel
		#endif
		{
			double modulusBest = 0;
			Atom atomBest;

			#ifdef _OPENMP
			#pragma omp for
			#endif
			for (size_t fIndex=0; fIndex<map->fCount; ++fIndex) {
				const double f = map->f(fIndex);
				const double expFactor = exp(-2*M_PI*s*s*f*f);
				for (size_t tIndex=0; tIndex<map->tCount; ++tIndex) {
					complex value = map->value(fIndex, tIndex);
					double phase = std::arg(value);
					double realFactor = 1.0 + cos(2*phase) * expFactor;
					double modulus = std::abs(value) * M_SQRT2 * sqrt(freqSampling / realFactor);
					if (modulus > modulusBest) {
						atomBest.type = ATOM_GABOR;
						atomBest.params.resize(6);
						atomBest.params[0] = modulus;
						atomBest.params[1] = std::abs(value) * amplFactor / realFactor;
						atomBest.params[2] = map->t(tIndex) * freqSampling;
						atomBest.params[3] = s * freqSampling;
						atomBest.params[4] = f * 2.0 / freqSampling;
						atomBest.params[5] = phase;
						modulusBest = modulus;
					}
				}
			}
			#ifdef _OPENMP
			#pragma omp critical
			#endif
			{
				if (modulusBest > modulusTotal) {
					atomTotal = atomBest;
					modulusTotal = modulusBest;
				}
			}
		}
	}
	return atomTotal;
}

size_t GaborWorkspace::getAtomCount(void) const {
	size_t atomCount = 0;
	for (auto map : maps) {
		atomCount += map->fCount * map->tCount;
	}
	return atomCount;
}

void GaborWorkspace::subtractAtom(const Atom& atom) {
	const double sA = atom.params[3] / freqSampling;
	const double fA = atom.params[4] * freqSampling / 2.0;
	const double tA = atom.params[2] / freqSampling;
	complex V = std::polar(atom.params[1] * sqrt(sA / M_SQRT2), atom.params[5]);
	for (auto map : maps) {
		GaborProductBase base(sA, map->s);
		if (base.product >= ORTHOGONALITY) {
			#ifdef _OPENMP
			#pragma omp parallel for schedule(static,1)
			#endif
			for (size_t tIndex=0; tIndex<map->tCount; ++tIndex) {
				GaborProductCalculator calc(&base, tA, map->t(tIndex));
				if (calc.product >= ORTHOGONALITY) {
					for (size_t fIndex=0; fIndex<map->fCount; ++fIndex) {
						const double f0 = map->f(fIndex);
						complex prodP = calc.calculate(fA, f0);
						complex prodN = calc.calculate(-fA, f0);
						map->value(fIndex, tIndex) -= 0.5 * (V * prodP + std::conj(V) * prodN);
					}
				}
			}
		}
	}
}

//------------------------------------------------------------------------------

Workspace* GaborWorkspaceBuilder::buildWorkspace(const SingleSignal& signal) const {
	double scaleMax = signal.samples.size() / signal.freqSampling;
	double root = sqrt(-2.0/M_PI * log(1.0-energyError));
	double aDenomSqrt = 1.0 - energyError;
	double aNominPart = energyError*(2.0-energyError)*(energyError*energyError-2*energyError+2);
	double a = (1.0 + sqrt(aNominPart)) / (aDenomSqrt * aDenomSqrt);

	const long N = signal.samples.size();
	const double tMax = (N-1) / signal.freqSampling;

	std::vector<std::shared_ptr<TimeFreqMap<complex>>> maps;
	std::vector<std::shared_ptr<GaborMapGenerator>> generators;
	int count = 0;
	for (double s=scaleMin; s<=scaleMax; s*=a) {
		const double dt = root * s;
		const double df = root / s;
		const double hwGabor = 3.0 * s; // TODO

		const long Ngauss = lrint(2.0 * hwGabor * signal.freqSampling - 0.5) + 1;

		const long Nfft1 = fftwRound( std::max(Ngauss, lrint(signal.freqSampling/df + 0.5)) );
		const long tCount1 = lrint(tMax/dt + 0.5) + 1;
		const long fCount1 = Nfft1 / 2;

		/* GaborMapGenerator2 is optimized for short (non-oscillating) atoms
		   but is quite buggy at the moment. */
//		const long Nfft2 = fftwRound(N + Ngauss);
//		const long tCount2 = N;
//		const long fCount2 = lrint(0.5 * signal.freqSampling / df + 0.5);

		TimeFreqMap<complex>* map;
		GaborMapGenerator* generator;
//		if (tCount1 * fCount1 <= tCount2 * fCount2) {
			map = new TimeFreqMap<complex>(s, fCount1, tCount1);
			generator = new GaborMapGenerator1(Ngauss, Nfft1, hwGabor, tMax, signal.freqSampling);
//		} else {
//			map = new TimeFreqMap<complex>(s, fCount2, tCount2);
//			generator = new GaborMapGenerator2(Ngauss, Nfft2, hwGabor, tMax, signal.freqSampling);
//		}
		maps.push_back(std::shared_ptr<TimeFreqMap<complex>>(map));
		generators.push_back(std::shared_ptr<GaborMapGenerator>(generator));
		++count;
	}

	for (int i=0; i<count; ++i) {
		generators[i]->generate(*maps[i], signal);
	}
	return new GaborWorkspace(signal.freqSampling, std::move(maps));
}

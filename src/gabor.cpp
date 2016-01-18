/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "gabor.hpp"

/**
 * Minimum scale parameter for Gabor atoms, in number of samples.
 */
const int MIN_SCALE_IN_SAMPLES = 50;

/**
 * If the scalar product between two normalized Gabor atoms is smaller
 * than this value, the atoms are classified as orthogonal.
 */
const double ORTHOGONALITY = 1.0e-8;

/**
 * Multiply given array by discretely sampled Gabor function,
 * i.e. transform x[t] := x[t] g[t], where g[t] ~ exp(−π(t-t₀)²/s²).
 * The function is L²-normalized: Σ g[t]² Δt = 1.
 *
 * @param buffer  array values x[0], x[Δt], x[2Δt], ...
 * @param N       length of the array
 * @param step    sampling step Δt (inverse of sampling frequency)
 * @param center  center t₀ of Gabor function
 * @param width   width s of Gabor function
 */
template<typename T>
static void gaussize(T* const __restrict buffer, const int N, double step, double center, double width) {
	const double norm = sqrt(M_SQRT2 / width);
	const double w = sqrt(M_PI) / width;
	const double a = w * step;
	const double b = w * center;
	for (int i=0; i<N; ++i) {
		const double x_width = a * i - b;
		buffer[i] *= norm * exp(-x_width * x_width);
	}
}

//------------------------------------------------------------------------------

GaborProductEstimator::GaborProductEstimator(double s1, double s2) :
s12(s1*s1), s22(s2*s2), A(s12 + s22),
part(M_SQRT2 * sqrt(s1 * s2 / A))
{ }

double GaborProductEstimator::estimate(double t1, double t2) const {
	double dt = t1 - t2;
	return part * exp(-M_PI/A * dt * dt);
}

//------------------------------------------------------------------------------

GaborWorkspaceMap::GaborWorkspaceMap(double s, size_t fCount, size_t tCount, double freqSampling, double tMax)
:	TimeFreqMap<complex>(fCount, tCount),
	Nfft((fCount-1)*2), freqSampling(freqSampling),
	input(Nfft), output(fCount),
	plan(Nfft, input, output, FFTW_ESTIMATE | FFTW_DESTROY_INPUT),
	s(s), atomCount(fCount * tCount)
{
	for (size_t ti=0; ti<tCount; ++ti) {
		t(ti) = ti * tMax / (tCount - 1);
	}
	for (size_t fi=0; fi<fCount; ++fi) {
		f(fi) = fi * freqSampling / Nfft;
	}
}

void GaborWorkspaceMap::compute(const SingleSignal& signal) {
	for (size_t tIndex=0; tIndex<tCount; ++tIndex) {
		compute(signal, tIndex);
	}
}

void GaborWorkspaceMap::compute(const SingleSignal& signal, size_t tIndex) {
	const long N = signal.samples.size();
	const double t0 = t(tIndex);
	const double hwGabor = 3.0 * s; // TODO
	long iL = std::max(0L, lrint((t0 - hwGabor) * signal.freqSampling));
	long iR = std::min(N-1, lrint((t0 + hwGabor) * signal.freqSampling));
	const double t0fixed = t0 - iL / signal.freqSampling;

	input.zero();
	for (long i=iL; i<=iR; ++i) {
		input[i-iL] = signal.samples[i];
	}
	gaussize(&input, iR-iL, 1.0/signal.freqSampling, t0fixed, s);
	plan.execute();

	const double norm = 1.0 / signal.freqSampling;
	const double mult = 2 * M_PI * t0fixed;
	for (size_t fIndex=0; fIndex<fCount; ++fIndex) {
		value(fIndex, tIndex) = norm * output[fIndex] * std::polar(1.0, f(fIndex) * mult);
	}
}

//------------------------------------------------------------------------------

Atom GaborWorkspace::findBestMatch() const {
	const double freqNyquist = 0.5 * freqSampling;
	double modulusTotal = 0;
	Atom atomTotal;
	size_t testedCount = 0;
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
				const bool isHighFreq = (f > 0.5 * freqNyquist);
				// for high frequencies, one can write
				// cos(2π(fN-Δf)(t-t₀)+φ) = (−1)^n cos(2πΔf(t-t₀)+(2πfNt₀−φ))
				const double f4Norm = isHighFreq ? (freqNyquist - f) : f;
				const double expFactor = exp(-2*M_PI*s*s*f4Norm*f4Norm);
				for (size_t tIndex=0; tIndex<map->tCount; ++tIndex) {
					complex value = map->value(fIndex, tIndex);
					double phase = std::arg(value);
					double phase4Norm = isHighFreq ? (2*M_PI*freqNyquist*map->t(tIndex) - phase) : phase;
					double realFactor = 1.0 + cos(2*phase4Norm) * expFactor;
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
		std::cout << "TESTED" << '\t' << (testedCount += map->atomCount) << std::endl;
	}

	return atomTotal;
}

size_t GaborWorkspace::getAtomCount(void) const {
	size_t atomCount = 0;
	for (auto map : maps) {
		atomCount += map->atomCount;
	}
	return atomCount;
}

void GaborWorkspace::subtractAtom(const Atom& atom, SingleSignal& signal) {
	const double sA = atom.params[3] / freqSampling;
	const double fA = atom.params[4] * freqSampling / 2.0;
	const double tA = atom.params[2] / freqSampling;

	const double amplitude = atom.params[1];
	const double omega = 2.0 * M_PI * fA;
	const double phase = atom.params[5];

	const size_t N = signal.samples.size();
	for (size_t i=0; i<N; ++i) {
		double t = i / freqSampling - tA;
		double t_width = t / sA;
		signal.samples[i] -= amplitude * exp(-M_PI * t_width * t_width) * cos(omega * t + phase);
	}

	const int count = maps.size();
	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for (int i=0; i<count; ++i) {
		GaborWorkspaceMap& map = *maps[i];
		GaborProductEstimator estimator(sA, map.s);
		if (estimator.part >= ORTHOGONALITY) {
			for (size_t tIndex=0; tIndex<map.tCount; ++tIndex) {
				if (estimator.estimate(tA, map.t(tIndex)) >= ORTHOGONALITY) {
					map.compute(signal, tIndex);
				}
			}
		}
	}
}

//------------------------------------------------------------------------------

Workspace* GaborWorkspaceBuilder::buildWorkspace(const SingleSignal& signal) const {
	double scaleMin = MIN_SCALE_IN_SAMPLES / signal.freqSampling;
	double scaleMax = signal.samples.size() / signal.freqSampling;
	double root = sqrt(-2.0/M_PI * log(1.0-energyError));
	double aDenomSqrt = 1.0 - energyError;
	double aNominPart = energyError*(2.0-energyError)*(energyError*energyError-2*energyError+2);
	double a = (1.0 + sqrt(aNominPart)) / (aDenomSqrt * aDenomSqrt);

	const long N = signal.samples.size();
	const double tMax = (N-1) / signal.freqSampling;

	std::vector<std::shared_ptr<GaborWorkspaceMap>> maps;
	int count = 0;
	for (double s=scaleMin; s<=scaleMax; s*=a) {
		const double dt = root * s;
		const double df = root / s;
		const double hwGabor = 3.0 * s; // TODO

		const long Ngauss = lrint(2.0 * hwGabor * signal.freqSampling - 0.5) + 1;

		const long Nfft = fftwRound( std::max(Ngauss, lrint(signal.freqSampling/df + 0.5)) );
		const long tCount = lrint(tMax/dt + 0.5) + 1;
		const long fCount = Nfft / 2 + 1;

		maps.push_back(std::make_shared<GaborWorkspaceMap>(s, fCount, tCount, signal.freqSampling, tMax));
		++count;
	}

	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for (int i=0; i<count; ++i) {
		maps[i]->compute(signal);
	}
	return new GaborWorkspace(signal.freqSampling, std::move(maps));
}

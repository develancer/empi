/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include "gabor.hpp"

/**
 * Minimum scale parameter for Gabor atoms, in number of samples.
 */
const double MIN_SCALE_IN_SAMPLES = 2.0;

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
static void gaussize(T* const __restrict buffer, const int N, double step, double center, double width, bool normalize) {
	const double norm = normalize ? sqrt(M_SQRT2 / width) : 1.0;
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

GaborWorkspaceMap::GaborWorkspaceMap(double s, int fCount, int tCount, double freqSampling, double tMax, int channelCount)
:	TimeFreqMap<complex>(channelCount, fCount, tCount),
	Nfft((fCount-1)*2), freqSampling(freqSampling),
	input(Nfft), output(fCount),
	plan(Nfft, input, output, FFTW_ESTIMATE | FFTW_DESTROY_INPUT),
	s(s), atomCount(fCount * tCount)
{
	for (int ti=0; ti<tCount; ++ti) {
		t(ti) = ti * tMax / (tCount - 1);
	}
	for (int fi=0; fi<fCount; ++fi) {
		f(fi) = fi * freqSampling / Nfft;
	}
}

void GaborWorkspaceMap::compute(const MultiSignal& signal) {
	for (int tIndex=0; tIndex<tCount; ++tIndex) {
		for (int cIndex=0; cIndex<cCount; ++cIndex) {
			compute(signal.channels[cIndex], cIndex, tIndex);
		}
	}
}

void GaborWorkspaceMap::compute(const SingleSignal& signal, int cIndex, int tIndex) {
	const long N = signal.samples.size();
	const double t0 = t(tIndex);
	const double hwGabor = 3.0 * s; // TODO
	// type casts are safe since we know signal length fits in int
	int iL = static_cast<int>( std::max(0L, lrint((t0 - hwGabor) * signal.freqSampling)) );
	int iR = static_cast<int>( std::min(N-1, lrint((t0 + hwGabor) * signal.freqSampling)) );
	const double t0fixed = t0 - iL / signal.freqSampling;

	input.zero();
	for (int i=iL; i<=iR; ++i) {
		input[i-iL] = signal.samples[i];
	}
	gaussize(&input, iR-iL+1, 1.0/signal.freqSampling, t0fixed, s, true);
	plan.execute();

	const double norm = 1.0 / signal.freqSampling;
	const double mult = 2 * M_PI * t0fixed;
	for (int fIndex=0; fIndex<fCount; ++fIndex) {
		value(cIndex, fIndex, tIndex) = norm * output[fIndex] * std::polar(1.0, f(fIndex) * mult);
	}
}

//------------------------------------------------------------------------------

void GaborWorkspace::compute(const MultiSignal& signal) {
	const int count = static_cast<int>(maps.size());
	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for (int i=0; i<count; ++i) {
		maps[i]->compute(signal);
	}
}

Atoms GaborWorkspace::findBestMatch(MultichannelConstraint constraint) const {
	const double freqNyquist = 0.5 * freqSampling;
	double squareTotal = 0;
	Atoms atomsTotal;
	size_t testedCount = 0;
	for (auto map : maps) {
		const double s = map->s;
		const double amplFactor = 2.0 * sqrt(M_SQRT2 / s);
		#ifdef _OPENMP
		#pragma omp parallel
		#endif
		{
			double squareBest = 0;
			Atoms atomsBest(map->cCount);
			std::vector<complex> values(map->cCount);
			std::vector<double> amplitudes(map->cCount), moduli(map->cCount), phases(map->cCount);

			#ifdef _OPENMP
			#pragma omp for
			#endif
			for (int fIndex=0; fIndex<map->fCount; ++fIndex) {
				const double f = map->f(fIndex);
				const bool isHighFreq = (f > 0.5 * freqNyquist);
				// for high frequencies, one can write
				// cos(2π(fN-Δf)(t-t₀)+φ) = (−1)^n cos(2πΔf(t-t₀)+(2πfNt₀−φ))
				const double f4Norm = isHighFreq ? (freqNyquist - f) : f;
				const double expFactor = exp(-2*M_PI*s*s*f4Norm*f4Norm);
				for (int tIndex=0; tIndex<map->tCount; ++tIndex) {
					for (int cIndex=0; cIndex<map->cCount; ++cIndex) {
						values[cIndex] = map->value(cIndex, fIndex, tIndex);
					}
					if (constraint && !isHighFreq) {
						(*constraint)(values);
					}
					double square = 0.0;
					for (int cIndex=0; cIndex<map->cCount; ++cIndex) {
						phases[cIndex] = std::arg(values[cIndex]);
						double phase4Norm = isHighFreq ? (2*M_PI*freqNyquist*map->t(tIndex) - phases[cIndex]) : phases[cIndex];
						double realFactor = 1.0 + cos(2*phase4Norm) * expFactor;
						moduli[cIndex] = std::abs(values[cIndex]) * M_SQRT2 * sqrt(freqSampling / realFactor);
						amplitudes[cIndex] = std::abs(values[cIndex]) * amplFactor / realFactor;
						square += moduli[cIndex] * moduli[cIndex];
					}
					if (square > squareBest) {
						for (int cIndex=0; cIndex<map->cCount; ++cIndex) {
							atomsBest[cIndex].type = ATOM_GABOR;
							atomsBest[cIndex].params.resize(6);
							atomsBest[cIndex].params[0] = moduli[cIndex];
							atomsBest[cIndex].params[1] = amplitudes[cIndex];
							atomsBest[cIndex].params[2] = map->t(tIndex) * freqSampling;
							atomsBest[cIndex].params[3] = s * freqSampling;
							atomsBest[cIndex].params[4] = f * 2.0 / freqSampling;
							atomsBest[cIndex].params[5] = phases[cIndex];
						}
						squareBest = square;
					}
				}
			}
			#ifdef _OPENMP
			#pragma omp critical
			#endif
			{
				if (squareBest > squareTotal) {
					atomsTotal = atomsBest;
					squareTotal = squareBest;
				}
			}
		}
		std::cout << "TESTED" << '\t' << (testedCount += map->atomCount) << std::endl;
	}

	return atomsTotal;
}

size_t GaborWorkspace::getAtomCount(void) const {
	size_t atomCount = 0;
	for (auto map : maps) {
		atomCount += map->atomCount;
	}
	return atomCount;
}

void GaborWorkspace::subtractAtomFromSignal(Atom& atom, SingleSignal& signal, bool fit) {
	const double sA = atom.params[3] / signal.freqSampling;
	const double fA = atom.params[4] * signal.freqSampling / 2.0;
	const double tA = atom.params[2] / signal.freqSampling;

	const double amplitude = fit ? 1.0 : atom.params[1];
	const double omega = 2.0 * M_PI * fA;
	const double phase = atom.params[5];

	const long N = signal.samples.size();
	const double hwGabor = 3.0 * sA; // TODO
	// type casts are safe since we know signal length fits in int
	int iL = static_cast<int>( std::max(0L, lrint((tA - hwGabor) * signal.freqSampling)) );
	int iR = static_cast<int>( std::min(N-1, lrint((tA + hwGabor) * signal.freqSampling)) );
	const int Nw = 1 + iR - iL;
	std::vector<double> waveform(Nw);
	for (int i=0; i<Nw; ++i) {
		double t = (i+iL) / signal.freqSampling - tA;
		waveform[i] = amplitude * cos(omega * t + phase);
	}
	gaussize(waveform.data(), Nw, 1.0/signal.freqSampling, tA-iL/signal.freqSampling, sA, false);

	if (fit) {
		double norm2 = 0.0;
		for (int i=0; i<Nw; ++i) {
			norm2 += waveform[i] * waveform[i];
		}
		double norm = sqrt(norm2);
		for (int i=0; i<Nw; ++i) {
			waveform[i] /= norm;
		}
		double modulus = 0.0;
		for (int i=0; i<Nw; ++i) {
			modulus += waveform[i] * signal.samples[iL+i];
		}
		for (int i=0; i<Nw; ++i) {
			waveform[i] *= modulus;
		}
		atom.params[0] = modulus;
		atom.params[1] = modulus / norm;
	}
	for (int i=0; i<Nw; ++i) {
		signal.samples[iL+i] -= waveform[i];
	}
}

void GaborWorkspace::subtractAtom(const Atom& atom, SingleSignal& signal, int channel) {
	const double sA = atom.params[3] / signal.freqSampling;
	const double tA = atom.params[2] / signal.freqSampling;
	GaborWorkspace::subtractAtomFromSignal(const_cast<Atom&>(atom), signal, false);

	const int count = maps.size();
	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for (int i=0; i<count; ++i) {
		GaborWorkspaceMap& map = *maps[i];
		GaborProductEstimator estimator(sA, map.s);
		if (estimator.part >= ORTHOGONALITY) {
			for (int tIndex=0; tIndex<map.tCount; ++tIndex) {
				if (estimator.estimate(tA, map.t(tIndex)) >= ORTHOGONALITY) {
					map.compute(signal, channel, tIndex);
				}
			}
		}
	}
}

//------------------------------------------------------------------------------

Workspace* GaborWorkspaceBuilder::prepareWorkspace(double freqSampling, int channelCount, int sampleCount) const {
	double scaleMin = MIN_SCALE_IN_SAMPLES / freqSampling;
	double scaleMax = sampleCount / freqSampling;
	double root = sqrt(-2.0/M_PI * log(1.0-energyError));
	double aDenomSqrt = 1.0 - energyError;
	double aNominPart = energyError*(2.0-energyError)*(energyError*energyError-2*energyError+2);
	double a = (1.0 + sqrt(aNominPart)) / (aDenomSqrt * aDenomSqrt);

	const double tMax = (sampleCount-1) / freqSampling;

	std::vector<std::shared_ptr<GaborWorkspaceMap>> maps;
	int count = 0;
	for (double s=scaleMin; s<=scaleMax; s*=a) {
		const double dt = root * s;
		const double df = root / s;
		const double hwGabor = 3.0 * s; // TODO

		const long Ngauss = lrint(2.0 * hwGabor * freqSampling - 0.5) + 1;

		const long Nfft = fftwRound( std::max(Ngauss, lrint(freqSampling/df + 0.5)) );
		const long tCount = lrint(tMax/dt + 0.5) + 1;
		if (std::max(Nfft, tCount) > static_cast<long>(std::numeric_limits<int>::max())) {
			throw Exception("signalFileIsTooLongForThisDecomposition");
		}
		const long fCount = Nfft / 2 + 1;

		maps.push_back(std::make_shared<GaborWorkspaceMap>(s, static_cast<int>(fCount), static_cast<int>(tCount), freqSampling, tMax, channelCount));
		++count;
	}
	return new GaborWorkspace(freqSampling, std::move(maps));
}

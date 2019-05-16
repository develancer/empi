/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include "gabor.hpp"
#include "timer.hpp"

/**
 * Reconstructed half-width of the triangular envelope,
 * in units of scale parameter (s).
 */
const double TRI_HALF_WIDTH = sqrt(2.5/M_PI);

/**
 * Minimum scale parameter for Gabor atoms, in number of samples.
 */
const double MIN_SCALE_IN_SAMPLES = 2.0;

/**
 * Multiply given array by discretely sampled Gabor function,
 * i.e. transform x[t] := x[t] g[t], where g[t] = exp(−π(t-t₀)²/s²).
 *
 * @param buffer  array values x[0], x[Δt], x[2Δt], ...
 * @param N       length of the array
 * @param step    sampling step Δt (inverse of sampling frequency)
 * @param center  center t₀ of Gabor function
 * @param width   width s of Gabor function
 */
template<typename T>
static void gaussize(T* const __restrict buffer, const int N, double step, double center, double width) {
	const double hw = sqrt(2.5/M_PI) * width;
	for (int i=0; i<N; ++i) {
		const double t = step * i - center;
		buffer[i] *= std::max(0.0, 1.0 - fabs(t)/hw);
	}
}

//------------------------------------------------------------------------------

static double I(double dl) {
    double x = exp(-dl/2);
    return (3 - x*x)*x/2;
}

static double J(double x) {
    x *= sqrt(10*M_PI);
    return 6/(x*x) * (1 - sin(x)/x);
}

static double K(double x) {
    x = fabs(x);
    return x*x <= 2.5/M_PI ? 1 - 0.6*M_PI*x*x + sqrt(9*M_PI*M_PI*M_PI/250)*x*x*x : (
            x*x <= 10/M_PI ? 2*pow(1 - x*sqrt(M_PI/10), 3) : 0
    );
}

//------------------------------------------------------------------------------

GaborComputer::GaborComputer(const GaborWorkspaceMap* map, int fIndex) :
	map(map), fIndex(fIndex),
	normComplex(std::sqrt(std::sqrt(0.9*M_PI) / map->s)),
	freqNyquist(0.5 * map->freqSampling),
	f(map->f(fIndex)),
	isHighFreq(f > 0.5 * freqNyquist),
	f4Norm(isHighFreq ? (freqNyquist - f) : f),
	expFactor(f4Norm>0 ? J(2*f4Norm*map->s) : 1.0)
{ }

double GaborComputer::compute(int tIndex, std::vector<complex>& buffer) {
	const double t = map->t(tIndex);
	for (int cIndex=0; cIndex<map->cCount; ++cIndex) {
		buffer[cIndex] = map->value(cIndex, fIndex, tIndex);
	}
	if (map->constraint && !isHighFreq) {
		(*map->constraint)(buffer);
	}
	double square = 0.0;
	for (int cIndex=0; cIndex<map->cCount; ++cIndex) {
		double phase = std::arg(buffer[cIndex]);
		double phase4Norm = isHighFreq ? (2*M_PI*freqNyquist*t - phase) : phase;
		double normReal = M_SQRT2 * normComplex / std::sqrt(1.0 + std::cos(2*phase4Norm) * expFactor);
		double moduli = std::abs(buffer[cIndex]) * normReal * std::sqrt(map->freqSampling);
		square += moduli * moduli;
	}
	return square;
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

Atoms GaborWorkspace::findBestMatch(void) const {
	TimeFreqMapValue best{nullptr, 0, 0, 0.0};

	const int mapCount = maps.size();
	#pragma omp parallel
	{
		TimeFreqMapValue bestPerThread{nullptr, 0, 0, 0.0};

		#pragma omp for schedule(dynamic,1)
		for (int i=0; i<mapCount; ++i) {
			TimeFreqValue peak = maps[i]->max();
			if (peak.value > bestPerThread.value) {
				bestPerThread = TimeFreqMapValue{maps[i], peak.fIndex, peak.tIndex, peak.value};
			}
		}

		#pragma omp critical
		{
			if (bestPerThread.value > best.value) {
				best = bestPerThread;
			}
		}
	}
	if (!best.map) {
		throw Exception("energyMapIsEmpty");
	}
	GaborComputer updater(best.map.get(), best.fIndex);
	std::vector<complex> values(best.map->cCount);
	updater.compute(best.tIndex, values);

	Atoms atomsBest(best.map->cCount);
	for (int cIndex=0; cIndex<best.map->cCount; ++cIndex) {
		double phase = std::arg(values[cIndex]);
		double phase4Norm = updater.isHighFreq ? (2*M_PI*updater.freqNyquist*best.map->t(best.tIndex) - phase) : phase;
		double normReal = M_SQRT2 * updater.normComplex / std::sqrt(1.0 + std::cos(2*phase4Norm) * updater.expFactor);
		double modulus = std::abs(values[cIndex]) * normReal * std::sqrt(freqSampling);
		double amplitude = std::abs(values[cIndex]) * normReal * normReal;

		atomsBest[cIndex].type = ATOM_GABOR;
		atomsBest[cIndex].params.resize(6);
		atomsBest[cIndex].params[0] = modulus;
		atomsBest[cIndex].params[1] = amplitude;
		atomsBest[cIndex].params[2] = best.map->t(best.tIndex) * freqSampling;
		atomsBest[cIndex].params[3] = best.map->s * freqSampling;
		atomsBest[cIndex].params[4] = updater.f * 2.0 / freqSampling;
		atomsBest[cIndex].params[5] = phase;
	}
	return atomsBest;
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
	const double hwGabor = TRI_HALF_WIDTH * sA;
	// type casts are safe since we know signal length fits in int
	int iL = static_cast<int>( std::max(0L, std::lrint((tA - hwGabor) * signal.freqSampling)) );
	int iR = static_cast<int>( std::min(N-1, std::lrint((tA + hwGabor) * signal.freqSampling)) );
	const int Nw = 1 + iR - iL;
	std::vector<double> waveform(Nw);
	for (int i=0; i<Nw; ++i) {
		double t = (i+iL) / signal.freqSampling - tA;
		waveform[i] = amplitude * std::cos(omega * t + phase);
	}
	gaussize(waveform.data(), Nw, 1.0/signal.freqSampling, tA-iL/signal.freqSampling, sA);

	if (fit) {
		double norm2 = 0.0;
		for (int i=0; i<Nw; ++i) {
			norm2 += waveform[i] * waveform[i];
		}
		double norm = std::sqrt(norm2);
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

	TIMER_START(subtractAtomFromSignal);
	GaborWorkspace::subtractAtomFromSignal(const_cast<Atom&>(atom), signal, false);
	TIMER_STOP(subtractAtomFromSignal);

	const int count = maps.size();
	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for (int i=0; i<count; ++i) {
		GaborWorkspaceMap& map = *maps[i];
		for (int tIndex=0; tIndex<map.tCount; ++tIndex) {
			if (fabs(map.t(tIndex) - tA) <= (map.s + sA) * TRI_HALF_WIDTH) {
				map.compute(signal, channel, tIndex);
			}
		}
	}
}

//------------------------------------------------------------------------------

static double solve(double (*F)(double), double result) {
    double Fx = F(result);
    double xL = result, xP = result;
    if (Fx > result) {
        while (Fx > result) {
            xL = xP;
            xP *= 2;
            Fx = F(xP);
        }
    } else if (Fx < result) {
        while (Fx < result) {
            xP = xL;
            xL /= 2;
            Fx = F(xL);
        }
    } else {
        return result; // lucky first guess
    }

    // szukamy pomiędzy xL a xP
    double x = (xP + xL) / 2;
    while (xL < x && x < xP) {
        Fx = F(x);
        if (Fx > result) {
            xL = x;
        } else if (Fx < result) {
            xP = x;
        } else {
            return x;
        }
        x = (xP + xL) / 2;
    }
    return x;
}

Workspace* GaborWorkspaceBuilder::prepareWorkspace(double freqSampling, int channelCount, int sampleCount, MultichannelConstraint constraint) const {
	double scaleMin = std::max(this->scaleMin, MIN_SCALE_IN_SAMPLES / freqSampling);
	double scaleMax = std::min(this->scaleMax, sampleCount / freqSampling);

	double dl = solve(I, 1.0 - energyError);
	double dfs = solve(J, 1.0 - energyError);
	double dt_s = solve(K, 1.0 - energyError);

	const double tMax = (sampleCount-1) / freqSampling;
	const double lMin = log(scaleMin);
	const double lMax = log(scaleMax);
	long lCount = std::lrint((lMax - lMin)/dl + 0.5) + 1;

	std::vector<std::shared_ptr<GaborWorkspaceMap>> maps;
	int count = 0;
	for (int i=0; i<lCount; ++i) {
		const double s = exp(lMin + i*(lMax - lMin)/(lCount - 1));
		const double dt = dt_s * s;
		const double df = dfs / s;
		const double hwGabor = TRI_HALF_WIDTH * s;

		const long Ngauss = 2 * std::lrint(hwGabor * freqSampling - 0.5) + 1;

		const long Nfft = fftwRound( std::max(Ngauss, std::lrint(freqSampling/df + 0.5)) );
		const long tCount = std::lrint(tMax/dt + 0.5) + 1;
		if (std::max(Nfft, tCount) > static_cast<long>(std::numeric_limits<int>::max())) {
			throw Exception("signalFileIsTooLongForThisDecomposition");
		}
		long fCount = Nfft / 2 + 1;
		if (std::isfinite(freqMax) && freqMax > 0) {
			fCount = std::min(fCount, std::lrint(freqMax/freqSampling*Nfft + 0.5) + 1);
		}

		maps.push_back(std::make_shared<GaborWorkspaceMap>(s, static_cast<int>(Nfft), static_cast<int>(fCount), static_cast<int>(tCount), freqSampling, tMax, channelCount, constraint));
		++count;
	}

	return new GaborWorkspace(freqSampling, std::move(maps));
}

//------------------------------------------------------------------------------

GaborWorkspaceMap::GaborWorkspaceMap(double s, int Nfft, int fCount, int tCount, double freqSampling, double tMax, int channelCount, MultichannelConstraint constraint)
:	TimeFreqMap<complex>(channelCount, fCount, tCount),
	Nfft(Nfft), input(Nfft), output(Nfft/2+1),
	plan(Nfft, input, output, FFTW_ESTIMATE | FFTW_DESTROY_INPUT),
	s(s), freqSampling(freqSampling), atomCount(fCount * tCount), constraint(constraint)
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
	index.reset( new GaborWorkspaceIndex(this) );
}

void GaborWorkspaceMap::compute(const SingleSignal& signal, int cIndex, int tIndex) {
	const long N = signal.samples.size();
	const double t0 = t(tIndex);
	const double hwGabor = TRI_HALF_WIDTH * s;
	// type casts are safe since we know signal length fits in int
	int iL = static_cast<int>( std::max(0L, std::lrint((t0 - hwGabor) * signal.freqSampling)) );
	int iR = static_cast<int>( std::min(N-1, std::lrint((t0 + hwGabor) * signal.freqSampling)) );
	const double t0fixed = t0 - iL / signal.freqSampling;

	input.zero();
	for (int i=iL; i<=iR; ++i) {
		input[i-iL] = signal.samples[i];
	}
	gaussize(&input, iR-iL+1, 1.0/signal.freqSampling, t0fixed, s);
	plan.execute();

	const double norm = 1.0 / signal.freqSampling;
	const double mult = 2 * M_PI * t0fixed;
	for (int fIndex=0; fIndex<fCount; ++fIndex) {
		value(cIndex, fIndex, tIndex) = norm * output[fIndex] * std::polar(1.0, f(fIndex) * mult);
	}
	if (index) {
		index->mark(tIndex);
	}
}

TimeFreqValue GaborWorkspaceMap::max(void) {
	if (!index) {
		throw Exception("internalLogicError");
	}
	return index->max();
}

//------------------------------------------------------------------------------

void GaborWorkspaceIndex::update(void) {
	for (int fIndex=0; fIndex<map->fCount; ++fIndex) {
		GaborComputer updater(map, fIndex);
		for (int tIndex=0; tIndex<map->tCount; ++tIndex) if (tIndexNeedUpdate[tIndex]) {
			size_t key = fIndex * static_cast<size_t>(map->tCount) + tIndex;
			double value = updater.compute(tIndex, buffer);
			heap->update(key, value);
		}
	}
	std::fill(tIndexNeedUpdate.begin(), tIndexNeedUpdate.end(), false);
	tIndexNeedUpdateFlag = false;
}

GaborWorkspaceIndex::GaborWorkspaceIndex(const GaborWorkspaceMap* map)
: map(map), buffer(map->cCount) {
	tIndexNeedUpdate.resize(map->tCount);
	std::vector<double> data(map->tCount * map->fCount);
	size_t key = 0;
	for (int fIndex=0; fIndex<map->fCount; ++fIndex) {
		GaborComputer updater(map, fIndex);
		for (int tIndex=0; tIndex<map->tCount; ++tIndex) {
			data[key++] = updater.compute(tIndex, buffer);
		}
	}
	this->heap.reset(new Heap<double>(std::move(data)));
}

void GaborWorkspaceIndex::mark(int tIndex) {
	tIndexNeedUpdate[tIndex] = true;
	tIndexNeedUpdateFlag = true;
}

TimeFreqValue GaborWorkspaceIndex::max(void) {
	if (tIndexNeedUpdateFlag) {
		update();
	}
	HeapItem<double> peak = heap->peek();
	int fIndex = peak.key / map->tCount;
	int tIndex = peak.key % map->tCount;
	return TimeFreqValue{fIndex, tIndex, peak.value};
}

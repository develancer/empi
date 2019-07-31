/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_GABOR_HPP
#define	EMPI_GABOR_HPP

#include <memory>
#include "classes.hpp"
#include "fftw.hpp"
#include "heap.hpp"

class GaborWorkspace;
class GaborWorkspaceBuilder;
class GaborWorkspaceIndex;
class GaborWorkspaceMap;

struct TimeFreqValue {
	int fIndex, tIndex;
	double value;
};

struct TimeFreqMapValue {
	std::shared_ptr<GaborWorkspaceMap> map;
	int fIndex, tIndex;
	double value;
};

struct GaborComputer {
	const GaborWorkspaceMap* map;
	const int fIndex;
	const double normComplex;
	const double freqNyquist;
	const double f;
	const bool isHighFreq;
	const double f4Norm;
	const double expFactor;

	GaborComputer(const GaborWorkspaceMap* map, int fIndex);
	double compute(int tIndex, std::vector<complex>& buffer);
};

//------------------------------------------------------------------------------

class GaborWorkspace : public Workspace {
	const double freqSampling;
	std::vector< std::shared_ptr<GaborWorkspaceMap> > maps;

public:
	static void subtractAtomFromSignal(Atom& atom, SingleSignal& signal, bool fit);

	GaborWorkspace(double freqSampling, std::vector< std::shared_ptr<GaborWorkspaceMap> >&& maps)
	: freqSampling(freqSampling), maps(std::move(maps)) { }

	void compute(const MultiSignal& signal);

	Atoms findBestMatch(void) const;

	size_t getAtomCount(void) const;

	void subtractAtom(const Atom& atom, SingleSignal& signal, int channel);
};

//------------------------------------------------------------------------------

class GaborWorkspaceBuilder : public WorkspaceBuilder {
	const double energyError;
	const double scaleMin;
	const double scaleMax;
	const double freqMax;

public:
	GaborWorkspaceBuilder(double energyError, double scaleMin, double scaleMax, double freqMax)
	: energyError(energyError), scaleMin(scaleMin), scaleMax(scaleMax), freqMax(freqMax) { }

	Workspace* prepareWorkspace(double freqSampling, int channelCount, int sampleCount, MultichannelConstraint constraint) const;
};

//------------------------------------------------------------------------------

class GaborWorkspaceIndex {

	const GaborWorkspaceMap* const map;
	std::vector<complex> buffer;
	std::unique_ptr< Heap<double> > heap;

	std::vector<bool> tIndexNeedUpdate;
	bool tIndexNeedUpdateFlag;

	void update(void);

public:
	GaborWorkspaceIndex(const GaborWorkspaceMap* map);
	void mark(int tIndex);
	TimeFreqValue max(void);
};

//------------------------------------------------------------------------------

class GaborWorkspaceMap : public TimeFreqMap<complex> {
	const int Nfft;
	fftwDouble input;
	fftwComplex output;
	fftwPlan plan;
	std::unique_ptr<GaborWorkspaceIndex> index;

public:
	const double s;
	const double freqSampling;
	const size_t atomCount;
	const MultichannelConstraint constraint;

	std::vector<double> expFactors;

	GaborWorkspaceMap(double s, int Nfft, int fCount, int tCount, double freqSampling, double tMax, int channelCount, MultichannelConstraint constraint);

	void compute(const MultiSignal& signal);

	void compute(const SingleSignal& signal, int cIndex, int tIndex);

	TimeFreqValue max(void);

	GaborWorkspaceMap(const GaborWorkspaceMap&) =delete;
	void operator=(const GaborWorkspaceMap&) =delete;
};

//------------------------------------------------------------------------------

#endif	/* EMPI_GABOR_HPP */

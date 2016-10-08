/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_CLASSES_HPP
#define	EMPI_CLASSES_HPP

#include "base.hpp"

template<typename T>
class TimeFreqMap {
	T* pointer = 0;
	double* fValues = 0;
	double* tValues = 0;

	TimeFreqMap(const TimeFreqMap&); // forbidden
	void operator=(const TimeFreqMap&); // forbidden

public:	
	const int cCount, fCount, tCount;

	TimeFreqMap(int cCount, int fCount, int tCount)
	: cCount(cCount), fCount(fCount), tCount(tCount) {
		pointer = static_cast<T*>(calloc(fCount * tCount * cCount, sizeof(T)));
		fValues = static_cast<double*>(calloc(fCount, sizeof(double)));
		tValues = static_cast<double*>(calloc(tCount, sizeof(double)));
	}

	~TimeFreqMap() {
		free(pointer);
		free(fValues);
		free(tValues);
	}

	double& f(int fIndex) {
		return fValues[fIndex];
	}

	const double& f(int fIndex) const {
		return fValues[fIndex];
	}

	double& t(int tIndex) {
		return tValues[tIndex];
	}

	const double& t(int tIndex) const {
		return tValues[tIndex];
	}

	T& value(int cIndex, int fIndex, int tIndex) {
		return pointer[(cIndex * fCount + fIndex) * tCount + tIndex];
	}

	const T& value(int cIndex, int fIndex, int tIndex) const {
		return pointer[(cIndex * fCount + fIndex) * tCount + tIndex];
	}
};

//------------------------------------------------------------------------------

typedef void (*MultichannelConstraint)(std::vector<complex>&);

class Workspace {
public:
	static void subtractAtomFromSignal(Atom& atom, SingleSignal& signal, bool fit);

	virtual ~Workspace() =default;
	virtual void compute(const MultiSignal& signal) =0;
	virtual Atoms findBestMatch(MultichannelConstraint constraint = nullptr) const =0;
	virtual size_t getAtomCount(void) const =0;
	virtual void subtractAtom(const Atom& atom, SingleSignal& signal, int channel) =0;
};

class WorkspaceBuilder {
public:
	virtual ~WorkspaceBuilder() =default;
	virtual Workspace* prepareWorkspace(double freqSampling, int channelCount, int sampleCount) const =0;
};

//------------------------------------------------------------------------------

struct DecompositionSettings {
	int iterationMax;
	double residualEnergy;
};

class Decomposition {
	const MultichannelConstraint constraint;

protected:
	Decomposition(MultichannelConstraint constraint)
	: constraint(constraint) { }

public:
	virtual ~Decomposition() =default;
	virtual MultiChannelResult compute(const DecompositionSettings& settings, Workspace* workspace, const MultiSignal& signal);
};

//------------------------------------------------------------------------------

class SmpDecomposition : public Decomposition {
public:
	SmpDecomposition(void) : Decomposition(nullptr) { }

	MultiChannelResult compute(const DecompositionSettings& settings, Workspace* workspace, const MultiSignal& signal);
};

#endif /* EMPI_CLASSES_HPP */

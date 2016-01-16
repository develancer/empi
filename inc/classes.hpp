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
	const size_t fCount, tCount;

	TimeFreqMap(size_t fCount, size_t tCount)
	: fCount(fCount), tCount(tCount) {
		pointer = static_cast<T*>(calloc(fCount * tCount, sizeof(T)));
		fValues = static_cast<double*>(calloc(fCount, sizeof(double)));
		tValues = static_cast<double*>(calloc(tCount, sizeof(double)));
	}

	~TimeFreqMap() {
		free(pointer);
		free(fValues);
		free(tValues);
	}

	double& f(size_t fIndex) {
		return fValues[fIndex];
	}

	const double& f(size_t fIndex) const {
		return fValues[fIndex];
	}

	double& t(size_t tIndex) {
		return tValues[tIndex];
	}

	const double& t(size_t tIndex) const {
		return tValues[tIndex];
	}

	T& value(size_t fIndex, size_t tIndex) {
		return pointer[fIndex * tCount + tIndex];
	}

	const T& value(size_t fIndex, size_t tIndex) const {
		return pointer[fIndex * tCount + tIndex];
	}
};

//------------------------------------------------------------------------------

class Workspace {
public:
	virtual ~Workspace() =default;
	virtual Atom findBestMatch() const =0;
	virtual size_t getAtomCount(void) const =0;
	virtual void subtractAtom(const Atom& atom, SingleSignal& signal) =0;
};

class WorkspaceBuilder {
public:
	virtual ~WorkspaceBuilder() =default;
	virtual Workspace* buildWorkspace(const SingleSignal&) const =0;
};

//------------------------------------------------------------------------------

struct DecompositionSettings {
	int iterationMax;
	double residualEnergy;
};

class Decomposition {
protected:
	static SingleChannelResult compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const SingleSignal& signal);

public:
	virtual ~Decomposition() =default;
	virtual MultiChannelResult compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const MultiSignal& signal) =0;
};

//------------------------------------------------------------------------------

class SmpDecomposition : public Decomposition {
public:
	MultiChannelResult compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const MultiSignal& signal);
};

#endif /* EMPI_CLASSES_HPP */

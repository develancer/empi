/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#ifndef EMPI_GABOR_HPP
#define	EMPI_GABOR_HPP

#include <memory>
#include "classes.hpp"
#include "fftw.hpp"

struct GaborProductEstimator {
	const double s12, s22, A, part;

	GaborProductEstimator(double s1, double s2);

	double estimate(double t1, double t2) const;
};

class GaborWorkspaceMap : public TimeFreqMap<complex> {
	const int Nfft;
	const double freqSampling;
	fftwDouble input;
	fftwComplex output;
	fftwPlan plan;

	GaborWorkspaceMap(const GaborWorkspaceMap&); // forbidden
	void operator=(const GaborWorkspaceMap&); // forbidden

public:
	const double s;

	GaborWorkspaceMap(double s, size_t fCount, size_t tCount, double freqSampling, double tMax);

	void compute(const SingleSignal& signal);

	void compute(const SingleSignal& signal, size_t tIndex);
};

class GaborWorkspace : public Workspace {
	const double freqSampling;
	std::vector<std::shared_ptr<GaborWorkspaceMap>> maps;
	std::vector<double> buffer;

public:
	GaborWorkspace(double freqSampling, std::vector<std::shared_ptr<GaborWorkspaceMap>>&& maps, int N)
	: freqSampling(freqSampling), maps(maps), buffer(N) { }

	Atom findBestMatch() const;

	size_t getAtomCount(void) const;

	void subtractAtom(const Atom& atom, SingleSignal& signal);
};

class GaborWorkspaceBuilder : public WorkspaceBuilder {
	double energyError;
	double scaleMin;

public:
	GaborWorkspaceBuilder(double energyError, double scaleMin)
	: energyError(energyError), scaleMin(scaleMin) { }

	Workspace* buildWorkspace(const SingleSignal&) const;

	void writeAtom(const Atom&, FILE*) const;
};

#endif	/* EMPI_GABOR_HPP */

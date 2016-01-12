/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#ifndef EMPI_GABOR_HPP
#define	EMPI_GABOR_HPP

#include <memory>
#include "classes.hpp"

class GaborWorkspace : public Workspace {
	const double freqSampling;
	std::vector<std::shared_ptr<TimeFreqMap<complex>>> maps;

public:
	static complex computeProduct(double s1, double t1, double f1, double s2, double t2, double f2);

	GaborWorkspace(double freqSampling, std::vector<std::shared_ptr<TimeFreqMap<complex>>>&& maps)
	: freqSampling(freqSampling), maps(maps) { }

	Atom findBestMatch() const;

	size_t getAtomCount(void) const;

	void subtractAtom(const Atom&);
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

class GaborMapGenerator {
protected:
	const long Ngauss, Nfft;
	const double hwGabor, tMax, freqSampling;

public:
	GaborMapGenerator(long Ngauss, long Nfft, double hwGabor, double tMax, double freqSampling)
	: Ngauss(Ngauss), Nfft(Nfft), hwGabor(hwGabor), tMax(tMax), freqSampling(freqSampling) { }

	virtual ~GaborMapGenerator() { }

	virtual void generate(TimeFreqMap<complex>& map, const SingleSignal& signal) =0;
};

class GaborMapGenerator0 : public GaborMapGenerator {
public:
	GaborMapGenerator0(long Ngauss, long Nfft, double hwGabor, double tMax, double freqSampling)
	: GaborMapGenerator(Ngauss, Nfft, hwGabor, tMax, freqSampling) { }

	void generate(TimeFreqMap<complex>& map, const SingleSignal& signal);
};

class GaborMapGenerator1 : public GaborMapGenerator {
public:
	GaborMapGenerator1(long Ngauss, long Nfft, double hwGabor, double tMax, double freqSampling)
	: GaborMapGenerator(Ngauss, Nfft, hwGabor, tMax, freqSampling) { }

	void generate(TimeFreqMap<complex>& map, const SingleSignal& signal);
};

class GaborMapGenerator2 : public GaborMapGenerator {
public:
	GaborMapGenerator2(long Ngauss, long Nfft, double hwGabor, double tMax, double freqSampling)
	: GaborMapGenerator(Ngauss, Nfft, hwGabor, tMax, freqSampling) { }

	void generate(TimeFreqMap<complex>& map, const SingleSignal& signal);
};

#endif	/* EMPI_GABOR_HPP */

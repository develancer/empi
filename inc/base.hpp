/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#ifndef EMPI_BASE_HPP
#define	EMPI_BASE_HPP

#include <complex>
#include <vector>

enum AtomType { ATOM_GABOR = 13 };

struct Atom {
	AtomType type;
	std::vector<double> params;
};

typedef std::complex<double> complex;
typedef std::vector<double> Samples;
typedef std::vector<Atom> SingleChannelResult;
typedef std::vector<SingleChannelResult> MultiChannelResult;

struct SingleSignal {
	double freqSampling;
	Samples samples;
};

struct MultiSignal {
	std::vector<SingleSignal> channels;
};

#endif	/* EMPI_BASE_HPP */

/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#ifndef EMPI_BASE_HPP
#define	EMPI_BASE_HPP

#include <complex>
#include <iostream>
#include <string>
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

	double computeEnergy() const {
		const size_t N = samples.size();
		double sum = 0.0;
		for (size_t i=0; i<N; ++i) {
			sum += samples[i] * samples[i];
		}
		return sum / freqSampling;
	}
};

struct MultiSignal {
	std::vector<SingleSignal> channels;
};

class Exception : public std::runtime_error {
public:
	explicit Exception(const std::string& camelCaseMessage)
	: std::runtime_error(camelCaseMessage) { }
};

#endif	/* EMPI_BASE_HPP */

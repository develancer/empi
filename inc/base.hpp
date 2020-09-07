/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BASE_HPP
#define	EMPI_BASE_HPP

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

enum AtomType { ATOM_GABOR = 13 };

struct Atom {
	AtomType type;
	std::vector<double> params;
};

typedef std::complex<double> complex;
typedef std::vector<Atom> Atoms;
typedef std::vector<Atoms> MultiChannelResult;

struct SingleSignal {
	double freqSampling;
	std::vector<double> samples;

	double computeEnergy() const {
		const int N = samples.size();
		double sum = 0.0;
		for (int i=0; i<N; ++i) {
			sum += samples[i] * samples[i];
		}
		return sum / freqSampling;
	}
};

struct MultiSignal {
	long segmentOffset; // in samples
	std::vector<SingleSignal> channels;

	double computeEnergy() const {
		double sum = 0.0;
		for (const SingleSignal& channel : channels) {
			sum += channel.computeEnergy();
		}
		return sum;
	}

	double getFreqSampling() const {
		return channels.empty() ? NAN : channels[0].freqSampling;
	}

	int getSampleCount() const {
		return channels.empty() ? 0 : channels[0].samples.size();
	}
};

class Exception : public std::runtime_error {
public:
	explicit Exception(const std::string& camelCaseMessage)
	: std::runtime_error(camelCaseMessage) { }
};

#endif	/* EMPI_BASE_HPP */

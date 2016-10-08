/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "gabor.hpp"
#include "mmp.hpp"

void Mmp1Decomposition::constraint(std::vector<complex>& values) {
	complex direction = 0;
	for (complex value : values) {
		direction += value * value;
	}
	double angle = 0.5 * std::arg(direction);
//	commented version does not work, should be equivalent
//	direction = std::polar(1.0, angle);
	for (complex& value : values) {
		value = std::polar(std::abs(value) * cos(std::arg(value) - angle), angle);
//		value = std::abs(value) * std::real(value * direction) * direction;
	}
}

MultiChannelResult Mmp2Decomposition::compute(const DecompositionSettings& settings, Workspace* workspace, const MultiSignal& signal) {
	MultiSignal sum;
	const int N = signal.getSampleCount();
	const int channelCount = signal.channels.size();
	sum.channels.push_back(signal.channels[0]);
	for (int c=1; c<channelCount; ++c) {
		for (int i=0; i<N; ++i) {
			sum.channels[0].samples[i] += signal.channels[c].samples[i];
		}
	}

	Atoms sumResult = Decomposition::compute(settings, workspace, sum)[0];

	MultiSignal residue(signal);
	MultiChannelResult result(channelCount);
	for (int c=0; c<channelCount; ++c) {
		for (Atom atom : sumResult) {
			Workspace::subtractAtomFromSignal(atom, residue.channels[c], true);
			result[c].push_back(atom);
		}
	}
	return result;
}

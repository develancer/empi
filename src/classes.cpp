/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <cstdio>
#include <memory>
#include "classes.hpp"
#include "gabor.hpp"

//------------------------------------------------------------------------------

void Workspace::subtractAtomFromSignal(Atom& atom, SingleSignal& signal, bool fit) {
	switch (atom.type) {
		case ATOM_GABOR:
			return GaborWorkspace::subtractAtomFromSignal(atom, signal, fit);
		default:
			throw Exception("invalidAtomGenerated");
	}
}

//------------------------------------------------------------------------------

MultiChannelResult Decomposition::compute(const DecompositionSettings& settings, Workspace* workspace, const MultiSignal& signal) {
	workspace->compute(signal);
	const int channelCount = signal.channels.size();
	MultiSignal residue(signal);
	MultiChannelResult result(channelCount);
	size_t atomCount = workspace->getAtomCount();
	const double totalEnergy = signal.computeEnergy();
	double residueEnergy = totalEnergy;
	if (totalEnergy == 0) {
		throw Exception("signalIsEmpty");
	}
	for (int iteration=1; iteration<=settings.iterationMax; ++iteration) {
		double progress = std::max(100.0 * (iteration - 1) / settings.iterationMax, 100.0 * (1.0 - residueEnergy / totalEnergy));
		std::cout << "ATOM\t" << (iteration - 1) << '\t' << atomCount << '\t' << progress << '\t' << progress << std::endl;
		Atoms bestMatches = workspace->findBestMatch(constraint);
		for (int c=0; c<channelCount; ++c) {
			result[c].push_back(bestMatches[c]);
		}
		if (iteration == settings.iterationMax) {
			break;
		}
		for (int c=0; c<channelCount; ++c) {
			workspace->subtractAtom(bestMatches[c], residue.channels[c], c);
		}
		residueEnergy = residue.computeEnergy();
		if (residueEnergy / totalEnergy <= settings.residualEnergy) {
			break;
		}
	}
	return result;
}

//------------------------------------------------------------------------------

MultiChannelResult SmpDecomposition::compute(const DecompositionSettings& settings, Workspace* workspace, const MultiSignal& signal) {
	MultiChannelResult result;
	int channelNumber = 0;
	for (const auto& channel : signal.channels) {
		MultiSignal wrapper;
		wrapper.channels.push_back(channel);
		std::cout << "CHANNEL\t" << channelNumber++ << std::endl;
		result.push_back(Decomposition::compute(settings, workspace, wrapper)[0]);
	}
	return result;
}

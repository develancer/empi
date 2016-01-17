/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <memory>
#include "classes.hpp"

//------------------------------------------------------------------------------

SingleChannelResult Decomposition::compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const SingleSignal& signal) {
	std::unique_ptr<Workspace> workspace( builder.buildWorkspace(signal) );
	SingleSignal residue(signal);
	SingleChannelResult result;
	size_t atomCount = workspace->getAtomCount();
	const double totalEnergy = signal.computeEnergy();
	double residueEnergy = totalEnergy;
	if (totalEnergy == 0) {
		throw Exception("signalIsEmpty");
	}
	for (int iteration=1; iteration<=settings.iterationMax; ++iteration) {
		double progress = std::max(100.0 * (iteration-1) / settings.iterationMax, 100.0 * (1.0 - residueEnergy / totalEnergy));
		std::cout << "ATOM" << '\t' << (iteration-1) << '\t' << atomCount << '\t' << progress << '\t' << progress << std::endl;
		Atom best = workspace->findBestMatch();
		result.push_back(best);
		if (iteration == settings.iterationMax) {
			break;
		}
		workspace->subtractAtom(best, residue);
		residueEnergy = residue.computeEnergy();
		if (residueEnergy / totalEnergy <= settings.residualEnergy) {
			break;
		}
	}
	return result;
}

//------------------------------------------------------------------------------

MultiChannelResult SmpDecomposition::compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const MultiSignal& signal) {
	MultiChannelResult result;
	int channelNumber = 0;
	for (const auto& channel : signal.channels) {
		std::cout << "CHANNEL" << '\t' << channelNumber++ << std::endl;
		fflush(stdout);
		result.push_back(Decomposition::compute(settings, builder, channel));
	}
	return result;
}

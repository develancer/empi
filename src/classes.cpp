/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#include <cstdio>
#include <memory>
#include "classes.hpp"

//------------------------------------------------------------------------------

SingleChannelResult Decomposition::compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const SingleSignal& signal) {
	std::unique_ptr<Workspace> workspace( builder.buildWorkspace(signal) );
	SingleChannelResult result;
	size_t atomCount = workspace->getAtomCount();
	for (int iteration=1; iteration<=settings.iterationMax; ++iteration) {
		printf("ATOM\t%d\t%zd\t%.2f\t%.2f\n", iteration-1, atomCount, 0.0, 0.0);
		Atom best = workspace->findBestMatch();
		result.push_back(best);
		if (iteration == settings.iterationMax) {
			break;
		}
		workspace->subtractAtom(best);
	}
	return result;
}

//------------------------------------------------------------------------------

MultiChannelResult SmpDecomposition::compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const MultiSignal& signal) {
	MultiChannelResult result;
	int channelNumber = 0;
	for (const auto& channel : signal.channels) {
		printf("CHANNEL\t%d\n", channelNumber++);
		result.push_back(Decomposition::compute(settings, builder, channel));
	}
	return result;
}

/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_MMP_HPP
#define	EMPI_MMP_HPP

#include "classes.hpp"

class Mmp1Decomposition : public Decomposition {
public:
	MultiChannelResult compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const MultiSignal& signal);
};

class Mmp2Decomposition : public Decomposition {
public:
	MultiChannelResult compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const MultiSignal& signal);
};

class Mmp3Decomposition : public Decomposition {
public:
	MultiChannelResult compute(const DecompositionSettings& settings, const WorkspaceBuilder& builder, const MultiSignal& signal);
};

#endif	/* EMPI_MMP_HPP */

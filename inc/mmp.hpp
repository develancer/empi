/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_MMP_HPP
#define	EMPI_MMP_HPP

#include "classes.hpp"

class Mmp1Decomposition : public Decomposition {
	static void constraint(std::vector<complex>& values);

public:
	Mmp1Decomposition(void) : Decomposition(&constraint)
	{ }
};

class Mmp2Decomposition : public Decomposition {
public:
	Mmp2Decomposition(void) : Decomposition(nullptr)
	{ }

	MultiChannelResult compute(const DecompositionSettings& settings, Workspace* workspace, const MultiSignal& signal);
};

class Mmp3Decomposition : public Decomposition {
public:
	Mmp3Decomposition(void) : Decomposition(nullptr)
	{ }
};

#endif	/* EMPI_MMP_HPP */

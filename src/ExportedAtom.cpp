/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cmath>
#include "ExportedAtom.h"

//////////////////////////////////////////////////////////////////////////////

ExportedAtom::ExportedAtom(double energy)
        : amplitude(NAN), energy(energy), frequency(NAN), phase(NAN), scale(NAN), position(NAN) {}

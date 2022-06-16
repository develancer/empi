/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "Family.h"
#include "GaussianFamily.h"
#include "TriangularFamily.h"

const std::map<std::string, std::shared_ptr<Family>> Family::ALL = {
        { "gauss", std::make_shared<GaussianFamily>() },
        { "triangular", std::make_shared<TriangularFamily>() },
};

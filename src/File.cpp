/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <stdexcept>
#include "File.h"

//////////////////////////////////////////////////////////////////////////////

File::File(const char *path, const char *mode)
        : std::shared_ptr<FILE>(fopen(path, mode), [](FILE *f) { f && fclose(f); }) {
    if (!get()) {
        throw std::runtime_error("failed to open file");
    }
}

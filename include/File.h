/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_FILE_H
#define EMPI_FILE_H

#include <cstdio>
#include <memory>

/**
 * Simple smart pointer to FILE based on std::shared_ptr.
 */
class File : public std::shared_ptr<FILE> {
public:
    /**
     * Create a new file handle.
     * If fopen() call fails, throw an exception.
     *
     * @param path path to the file
     * @param mode as in fopen()
     */
    File(const char *path, const char *mode);
};

#endif //EMPI_FILE_H

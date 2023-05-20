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
     * Create a null file handle.
     */
    File() = default;

    /**
     * Create a new file handle.
     * If fopen() call fails, throw an exception.
     *
     * @param path path to the file
     * @param mode as in fopen()
     */
    File(const char *path, const char *mode);
};

/**
 * Smart pointer to binary file in read mode.
 * Includes a method to get size of the file.
 */
class FileToRead : public File {
    const size_t file_size;

    static size_t read_file_size(FILE* file);

public:
    /**
     * Create a new file handle opened in read-only binary mode.
     * If fopen() call fails, throw an exception.
     *
     * @param path path to the file
     */
    explicit FileToRead(const char* path);

    /**
     * @return file size in bytes or 0 if file size could not be queried
     */
    [[nodiscard]] size_t get_file_size() const;
};

#endif //EMPI_FILE_H

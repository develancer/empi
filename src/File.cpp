/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
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

//////////////////////////////////////////////////////////////////////////////

FileToRead::FileToRead(const char *path)
        : File(path, "rb"), file_size(read_file_size(get())) {}

size_t FileToRead::read_file_size(FILE *file) {
    off_t end_position;
    if (fseeko(file, 0, SEEK_END)
        || (end_position = ftello(file)) < 0) {
        end_position = 0; // default value when file size could not be read
    }
    rewind(file);
    return end_position;
}

size_t FileToRead::get_file_size() const {
    return file_size;
}

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <ctime>
#include "Logger.h"

//////////////////////////////////////////////////////////////////////////////

void Logger::header(const char* type) {
    time_t now = time(nullptr);
    char ctime_buffer[26];
    ctime_r(&now, ctime_buffer);
    ctime_buffer[24] = 0;

    fprintf(stderr, "[%s] %s: ", ctime_buffer, type);
}

void Logger::message(const char* type, const char* string) {
    header(type);
    fputs(string, stderr);
    fputc('\n', stderr);
}
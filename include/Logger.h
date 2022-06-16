/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_LOGGER_H
#define EMPI_LOGGER_H

#include <cstdio>

class Logger
{
    static void header(const char* type);

    static void message(const char* type, const char* string);

    template<typename... VALUES>
    static void message(const char* type, const char* format, VALUES... parameters) {
        header(type);
        fprintf(stderr, format, parameters...);
        fputc('\n', stderr);
    }

public:
    Logger() = delete;

    template<typename... VALUES>
    static void error(const char* format, VALUES... parameters) {
        message("ERROR", format, parameters...);
    }

    template<typename... VALUES>
    static void info(const char* format, VALUES... parameters) {
        message("INFO", format, parameters...);
    }

    template<typename... VALUES>
    static void internal_error(const char* format, VALUES... parameters) {
        message("INTERNAL ERROR", format, parameters...);
    }
};

#endif //EMPI_LOGGER_H

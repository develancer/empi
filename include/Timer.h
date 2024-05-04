/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_TIMER_H
#define EMPI_TIMER_H

#include <chrono>
#include <memory>

class ElapsingTimer {
    const bool hasMessage_;

    std::chrono::time_point<std::chrono::steady_clock> start_;

    ElapsingTimer(const ElapsingTimer &);

    ElapsingTimer &operator=(const ElapsingTimer &);

public:
    ElapsingTimer();

    ElapsingTimer(const char *msg, va_list ap);

    ~ElapsingTimer();

    [[nodiscard]] float time() const;
};

class Timer {
    std::unique_ptr<ElapsingTimer> timer_;

public:
    Timer() =default;

    void start();

    void start(const char *msg, ...)
    #ifdef __GNUC__
    __attribute__ (( format (printf, 2, 3)))
    #endif
    ;

    [[nodiscard]] float time() const;

    void stop();
};

#endif //EMPI_TIMER_H

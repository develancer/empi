/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <cstdarg>
#include "Timer.h"

//////////////////////////////////////////////////////////////////////////////

ElapsingTimer::ElapsingTimer()
        : hasMessage_(false) {
    start_ = std::chrono::steady_clock::now();
}

ElapsingTimer::ElapsingTimer(const char *msg, va_list ap)
        : hasMessage_(msg) {
    if (hasMessage_) {
        vfprintf(stderr, msg, ap);
        fputs("... ", stderr);
        fflush(stderr);
    }

    start_ = std::chrono::steady_clock::now();
}

ElapsingTimer::~ElapsingTimer() {
    if (hasMessage_) {
        fprintf(stderr, "(%.3f s)\n", time());
        fflush(stderr);
    }
}

float ElapsingTimer::time() const {
    auto now = std::chrono::steady_clock::now();
    return static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(now - start_).count()) * 1.0e-6f;
}

//////////////////////////////////////////////////////////////////////////////

void Timer::start() {
    timer_ = std::make_unique<ElapsingTimer>();
}

void Timer::start(const char *msg, ...) {
    va_list ap;
    va_start(ap, msg);
    timer_ = std::make_unique<ElapsingTimer>(msg, ap);
    va_end(ap);
}

float Timer::time() const {
    return timer_ ? timer_->time() : 0.0f;
}

void Timer::stop() {
    timer_.reset();
}

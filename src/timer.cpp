/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <string>
#include <sys/time.h>

#include "timer.hpp"

//----------------------------------------------------------------------

SingleTimer::SingleTimer(void) {
	gettimeofday(&start_, 0);
}

double SingleTimer::time(void) const {
	struct timeval now;
	gettimeofday(&now,0);
	return (now.tv_usec - start_.tv_usec)*1.0e-6 + (now.tv_sec - start_.tv_sec);
}

//----------------------------------------------------------------------

Timer::Timer(void)
{ }

void Timer::start(void) {
	timer_.reset( new SingleTimer() );
}

double Timer::time(void) const {
	return time_;
}

void Timer::stop(void) {
	if (timer_) {
		time_ += timer_->time();
	}
	timer_.reset();
}

//----------------------------------------------------------------------

std::map<std::string, Timer> timers;

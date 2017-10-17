/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_TIMER_HPP
#define EMPI_TIMER_HPP

#include <map>
#include <memory>

class SingleTimer {
	struct timeval start_;

	SingleTimer(const SingleTimer&) =delete;
	SingleTimer& operator=(const SingleTimer&) =delete;

public:
	SingleTimer(void);

	double time(void) const;
};

class Timer {
	double time_;
	std::unique_ptr<SingleTimer> timer_;

public:
	Timer(void);

	void start(void);

	void stop(void);

	double time(void) const;
};

extern std::map<std::string, Timer> timers;

#ifdef NDEBUG
#define TIMER_START(NAME)
#define TIMER_STOP(NAME)
#define PRINT_TIMERS
#else
#define TIMER_START(NAME) timers[#NAME].start()
#define TIMER_STOP(NAME) timers[#NAME].stop()
#define PRINT_TIMERS for (const auto& stats: timers) std::cerr << stats.first << '\t' << stats.second.time() << std::endl
#endif

#endif  /* EMPI_TIMER_HPP */

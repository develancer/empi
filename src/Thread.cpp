/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <pthread.h>
#include "Thread.h"

void Thread::affix_to_cpu() {
    static unsigned cpu_index = 0;
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(cpu_index, &cpuset);
    pthread_setaffinity_np(native_handle(), sizeof(cpu_set_t), &cpuset);
    cpu_index = (cpu_index + 1) % std::thread::hardware_concurrency();
}

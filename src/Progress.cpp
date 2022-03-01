/**********************************************************
 * Piotr T. Różański (c) 2015-2022                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <iostream>
#include "Progress.h"
#include "Types.h"

Progress::Progress(int total_epoch_count)
        : total_epoch_count(total_epoch_count), epochs_completed(0), last_progress(0) {}

void Progress::print_progress() {
    double completed = epochs_completed;
    for (const auto &pair: progress_map) {
        completed += pair.second;
    }
    int progress = Types::floor<int>(100.0 * completed / total_epoch_count);
    if (progress > last_progress) {
        printf("%d%% completed (finished %d out of %d segments)\n", progress, epochs_completed, total_epoch_count);
        fflush(stdout);
        last_progress = progress;
    }
}

void Progress::epoch_started(EpochIndex index) {
    std::lock_guard<std::mutex> lock(mutex);
    progress_map[index] = 0;
}

void Progress::epoch_progress(EpochIndex index, double progress) {
    progress_map[index] = progress;
    std::lock_guard<std::mutex> lock(mutex);
    print_progress();
}

void Progress::epoch_finished(EpochIndex index) {
    std::lock_guard<std::mutex> lock(mutex);
    progress_map.erase(index);
    ++epochs_completed;
    print_progress();
}

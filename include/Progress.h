/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_PROGRESS_H
#define EMPI_PROGRESS_H

#include <map>
#include <mutex>
#include "EpochIndex.h"

class Progress {
    const int total_epoch_count;
    int epochs_completed;
    int last_progress;

    std::map<EpochIndex, double> progress_map;
    std::mutex mutex;

    void print_progress();

public:
    explicit Progress(int total_epoch_count);

    void epoch_started(EpochIndex index);

    void epoch_progress(EpochIndex index, double progress);

    void epoch_finished(EpochIndex index);
};

#endif //EMPI_PROGRESS_H

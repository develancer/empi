/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_SPECTROGRAM_LOOP_H
#define EMPI_SPECTROGRAM_LOOP_H

#include <memory>
#include "SpectrogramCalculator.h"
#include "SpectrogramRequest.h"
#include "TaskQueue.h"

/**
 * Runnable object used with SpectrogramCalculator instances to be run in separate consumer threads.
 */
class SpectrogramLoop {
    std::shared_ptr<TaskQueue<SpectrogramRequest>> task_queue;
    std::unique_ptr<SpectrogramCalculator> calculator;
    bool reverse_order;

public:
    SpectrogramLoop(std::shared_ptr<TaskQueue<SpectrogramRequest>> task_queue, std::unique_ptr<SpectrogramCalculator> calculator, bool reverse_order = false)
            : task_queue(std::move(task_queue)), calculator(std::move(calculator)), reverse_order(reverse_order) {}

    void operator()(bool wait = true) {
        SpectrogramRequest task;
        while (task_queue->get(task, wait, reverse_order)) {
            calculator->compute(task);
            task.interface->notify();
            task_queue->notify();
        }
    }
};

#endif //EMPI_SPECTROGRAM_LOOP_H

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_WORKER_LOOP_H
#define EMPI_WORKER_LOOP_H

#include <memory>
#include "SpectrogramRequest.h"
#include "TaskQueue.h"
#include "Worker.h"

/**
 * Runnable object used with Worker instances to be run in separate consumer threads.
 */
class WorkerLoop {
    std::shared_ptr<TaskQueue<SpectrogramRequest>> task_queue;
    std::unique_ptr<Worker> worker;
    bool reverse_order;

public:
    WorkerLoop(std::shared_ptr<TaskQueue<SpectrogramRequest>> task_queue, std::unique_ptr<Worker> worker, bool reverse_order = false)
            : task_queue(task_queue), worker(std::move(worker)), reverse_order(reverse_order) {}

    void operator()() {
        SpectrogramRequest task;
        while (task_queue->get(task, reverse_order)) {
            worker->compute(task);
            task.interface->notify();
            task_queue->notify();
        }
    }
};

#endif //EMPI_WORKER_LOOP_H

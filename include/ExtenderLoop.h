/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_EXTENDER_LOOP_H
#define EMPI_EXTENDER_LOOP_H

#include <memory>
#include "Atom.h"
#include "TaskQueue.h"

/**
 * Runnable object used internally in Computer to form separate consumer threads for BlockAtom::extend() calculations.
 */
class ExtenderLoop {
    std::shared_ptr<TaskQueue<BasicAtomPointer>> task_queue;

public:
    explicit ExtenderLoop(std::shared_ptr<TaskQueue<BasicAtomPointer>> task_queue)
            : task_queue(std::move(task_queue)) {}

    void operator()(bool wait = true) {
        BasicAtomPointer atom;
        while (task_queue->get(atom, wait)) {
            atom->extend(true); // ExtendedAtomPointer will be put into cache
            task_queue->notify();
        }
    }
};

#endif //EMPI_EXTENDER_LOOP_H

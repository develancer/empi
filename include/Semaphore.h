//
// Created by piotr on 09.02.2022.
//

#ifndef EMPI_SEMAPHORE_H
#define EMPI_SEMAPHORE_H

#include <mutex>
#include <condition_variable>

class Semaphore {
    std::mutex mutex;
    std::condition_variable condition;
    unsigned long state = 0;

public:
    explicit Semaphore(int state = 0) : state(state) {}

    Semaphore(const Semaphore &) = delete;

    Semaphore(Semaphore &&) = delete;

    void operator=(const Semaphore &) = delete;

    void release(int increment = 1) {
        std::lock_guard<std::mutex> lock(mutex);
        state += increment;
        if (increment == 1) {
            condition.notify_one();
        } else if (increment > 1) {
            condition.notify_all();
        }
    }

    void acquire(int decrement = 1) {
        std::unique_lock<std::mutex> lock(mutex);
        while (state < decrement) {
            condition.wait(lock);
        }
        state -= decrement;
    }
};

#endif //EMPI_SEMAPHORE_H

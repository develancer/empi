/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_TASK_QUEUE_H
#define EMPI_TASK_QUEUE_H

#include <atomic>
#include <condition_variable>
#include <list>
#include <mutex>

/**
 * Synchronized (thread-safe) task queue for producer-consumer pattern with multiple consumers.
 * All public methods are thread-safe.
 *
 * @tparam T type of elements (tasks) that will be stored in queue.
 */
template<typename T>
class TaskQueue {
    std::mutex mutex;
    std::list<T> tasks;
    bool terminated;
    size_t unfinished_count;

    std::condition_variable task_count_up_or_terminated, unfinished_count_down_or_terminated;

public:
    /**
     * Create a new queue with no tasks.
     */
    TaskQueue() : terminated(false), unfinished_count(0) {}

    /**
     * Add a task to the queue and wait until all tasks are finished.
     * This method should be called by the producer thread.
     */
    void put(const T &task) {
        std::unique_lock lock(mutex);
        unfinished_count++;
        tasks.push_back(task);
        lock.unlock();

        task_count_up_or_terminated.notify_one();
    }

    /**
     * Add all tasks from the given list to the queue.
     * This method should be called by the producer thread.
     */
    void put(std::list<T> &&list_of_tasks) {
        std::unique_lock lock(mutex);
        unfinished_count += list_of_tasks.size();
        tasks.splice(tasks.end(), list_of_tasks);
        lock.unlock();

        task_count_up_or_terminated.notify_all();
    }

    /**
     * Get a next task from the queue. If no tasks are available, wait until one became available and then return.
     * This method should be called by consumer threads.
     * After processing the task, each consumer should call notify() prior to get-ting the next one.
     *
     * @param task reference to store the task taken from the queue
     * @param wait if true, and the queue is empty, wait until the task is available or the queue terminates;
     * if false, return false immediately
     * @return true if task was taken from the queue, false if the queue was terminated instead
     * (or not available and wait=false)
     */
    bool get(T &task, bool wait = true, bool from_back = false) {
        return from_back
               ? custom_get(task, wait, &decltype(tasks)::back, &decltype(tasks)::pop_back)
               : custom_get(task, wait, &decltype(tasks)::front, &decltype(tasks)::pop_front);
    }

    /**
     * Notify the producers that the task previously obtained with get() is now completed.
     */
    void notify() {
        std::unique_lock lock(mutex);
        unfinished_count--;
        lock.unlock();
        unfinished_count_down_or_terminated.notify_all();
    }

    /**
     * Terminate the queue. All consumers and producers will exit.
     */
    void terminate() {
        std::unique_lock lock(mutex);
        terminated = true;
        lock.unlock();
        task_count_up_or_terminated.notify_all();
        unfinished_count_down_or_terminated.notify_all();
    }

    /**
     * Wait until all tasks are finished.
     * This method should be called by the producer thread.
     */
    void wait_for_tasks() {
        std::unique_lock lock(mutex);
        unfinished_count_down_or_terminated.wait(lock, [&] { return terminated || !unfinished_count; });
        if (unfinished_count) {
            throw std::logic_error("TaskQueue was terminated while waiting for tasks");
        }
        lock.unlock();
    }

private:
    bool custom_get(T &task, bool wait, T &(std::list<T>::*ref)(), void (std::list<T>::*pop)()) {
        std::unique_lock lock(mutex);
        if (wait) {
            task_count_up_or_terminated.wait(lock, [&] { return terminated || !tasks.empty(); });
        } else if (tasks.empty()) {
            return false;
        }
        bool result = false;
        if (!terminated) {
            task = (tasks.*ref)();
            (tasks.*pop)();
            result = true;
        }

        return result;
    }
};

#endif //EMPI_TASK_QUEUE_H

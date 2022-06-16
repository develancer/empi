/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_WORKER_H
#define EMPI_WORKER_H

#include <list>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>

#include "Array.h"
#include "BlockInterface.h"
#include "Dictionary.h"
#include "IndexRange.h"
#include "OptimizationMode.h"
#include "SpectrogramCalculator.h"
#include "SpectrogramLoop.h"
#include "SpectrogramRequest.h"
#include "TaskQueue.h"
#include "Thread.h"
#include "Types.h"

/**
 * The central point of the computation process.
 * Worker instances facilitate communication between Dictionary objects
 * (which have an internal cache pointing to the current best atom)
 * and SpectrogramCalculator objects performing the main part of calculations and allowing dictionaries to update.
 * Resulting decomposition can be accessed by repeated calls to get_next_atom().
 */
class Worker {
    Array2D<real> data;
    const OptimizationMode mode;
    IndexRange updated_index_range;

    std::shared_ptr<TaskQueue<BasicAtomPointer>> atom_queue;
    std::shared_ptr<TaskQueue<SpectrogramRequest>> spectrogram_queue;
    std::unique_ptr<SpectrogramLoop> primary_calculator_loop;
    std::vector<Thread> threads;

    std::list<std::unique_ptr<Dictionary>> dictionaries;

public:
    /**
     * Create a new Worker to analyze given multi-channel signal.
     * Dictionaries and workers must be added manually after constructing the Worker instance.
     *
     * @param data reference to multi-channel data of the analysed signal
     * @param cpu_threads number of CPU threads to use for global optimization
     */
    explicit Worker(Array2D<real> data, unsigned cpu_threads, OptimizationMode mode = OPTIMIZATION_DISABLED);

    Worker(Worker&& source) = default;
    Worker(const Worker&) = delete;
    void operator=(const Worker&) = delete;

    /**
     * Destroy the Worker instance by terminating the task queue and closing all worker threads.
     */
    ~Worker();

    /**
     * Associate a SpectrogramCalculator object with this Worker instance.
     * It will be used during following calls to get_next_atom().
     * Each calculator should only be added once.
     *
     * @param calculator smart pointer to SpectrogramCalculator instance
     * @param prefer_long_fft should be set to true in case of GPU calculators, false otherwise
     */
    void add_calculator(std::unique_ptr<SpectrogramCalculator> calculator, bool prefer_long_fft = false);

    /**
     * Associate a Dictionary object with this Worker instance.
     * It will be used during following calls to get_next_atom().
     * Each dictionary should only be added once.
     *
     * @param dictionary smart pointer to Dictionary instance
     */
    void add_dictionary(std::unique_ptr<Dictionary> dictionary);

    /**
     * Compute and return the atom that was a best match for the analyzed signal.
     * This atom will be subtracted from signal before returning,
     * so consecutive calls to get_next_atom() will provide consecutive atoms for the decomposition.
     *
     * @return best matching atom
     */
    ExtendedAtomPointer get_next_atom();

    /**
     * Create a list of all request templates that can be requested by this worker's dictionaries.
     */
    std::list<ProtoRequest> get_proto_requests();

    /**
     * Reset the internal state of all dictionaries.
     * This method should be called whenever the input signal has been replaced with a new signal segment
     * and computation should be started from scratch.
     */
    void reset();
};

#endif //EMPI_WORKER_H

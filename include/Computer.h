/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_COMPUTER_H
#define EMPI_COMPUTER_H

#include <list>
#include <memory>
#include <mutex>
#include <queue>

#include "Array.h"
#include "Dictionary.h"
#include "IndexRange.h"
#include "BlockInterface.h"
#include "SpectrogramRequest.h"
#include "Types.h"
#include "Worker.h"

/**
 * The central point of the computation process.
 * Computer instances facilitate communication between Dictionary objects
 * (which have an internal cache pointing to the current best atom)
 * and Worker objects performing the main part of calculations and allowing dictionaries to update.
 * Resulting decomposition can be accessed by repeated calls to get_next_atom().
 */
class Computer {
    Array2D<real> data;
    IndexRange updated_index_range;

    std::mutex requests_mutex;
    std::list<SpectrogramRequest> requests;

    std::list<std::unique_ptr<Dictionary>> dictionaries;
    std::list<std::unique_ptr<Worker>> calculators_cpu;
    std::list<std::unique_ptr<Worker>> calculators_gpu;

    void work(Worker *calculator);

    void work_gpu(Worker *calculator);

    bool next_cpu(SpectrogramRequest &request);

    bool next_gpu(SpectrogramRequest &request);

public:
    /**
     * Create a new Computer to analyze given multi-channel signal.
     * Dictionaries and workers must be added manually after constructing the Computer instance.
     *
     * @param data reference to multi-channel data of the analysed signal
     */
    explicit Computer(Array2D<real> data) : data(std::move(data)) {
        reset();
    }

    /**
     * Associate a Worker object with this Computer instance.
     * It will be used during following calls to get_next_atom().
     * Each worker should only be added once.
     *
     * @param calculator smart pointer to Worker instance
     * @param prefer_long_fft should be set to true in case of GPU workers, false otherwise
     */
    void add_calculator(std::unique_ptr<Worker> calculator, bool prefer_long_fft = false);

    /**
     * Associate a Dictionary object with this Computer instance.
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
     * Create a list of all request templates that can be requested by this computer's dictionaries.
     */
    std::list<ProtoRequest> get_proto_requests();

    /**
     * Reset the internal state of all dictionaries.
     * This method should be called whenever the input signal has been replaced with a new signal segment
     * and computation should be started from scratch.
     */
    void reset();
};

#endif //EMPI_COMPUTER_H
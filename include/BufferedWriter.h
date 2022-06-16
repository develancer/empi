/**********************************************************
 * Piotr T. Różański (c) 2015–2022                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BUFFERED_WRITER_H
#define EMPI_BUFFERED_WRITER_H

#include <list>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>
#include "Array.h"
#include "BookWriter.h"
#include "Semaphore.h"

struct BufferedResult {
    std::mutex mutex;
    EpochIndex index;
    Array2D<double> data;
    std::vector<std::list<ExportedAtom>> atoms;
};

/**
 * Wrapper for book writer implementations to allow a thread-safe access to it.
 * Instances of this class can be used by multiple threads and epochs may be passed to write in any sequence.
 * Moreover, in case of the single-channel decomposition, individual channels can be passed as well.
 */
class BufferedWriter : public BookWriter {
    const int total_channel_count;
    const int real_epoch_count;
    std::unique_ptr<BookWriter> actual_writer;

    std::vector<Semaphore> semaphores;
    std::vector<BufferedResult> results;

    BufferedResult& get_result_ref(int epoch_counter, index_t sample_count);

public:
    BufferedWriter(int total_channel_count, int real_epoch_count, std::unique_ptr<BookWriter> actual_writer);

    void finalize() final;

    void write(Array2D<double> data, EpochIndex index, const std::vector<std::list<ExportedAtom>> &atoms) final;
};

#endif //EMPI_BUFFERED_WRITER_H

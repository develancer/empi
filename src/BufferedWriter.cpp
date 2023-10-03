/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "BufferedWriter.h"

//////////////////////////////////////////////////////////////////////////////

BufferedWriter::BufferedWriter(int total_channel_count, int real_epoch_count, std::unique_ptr<BookWriter> actual_writer)
        : BookWriter(*actual_writer), total_channel_count(total_channel_count), real_epoch_count(real_epoch_count),
        actual_writer(std::move(actual_writer)),
        semaphores(real_epoch_count), results(real_epoch_count) {}

BufferedResult& BufferedWriter::get_result_ref(int epoch_counter, index_t sample_count) {
    BufferedResult& result = results[epoch_counter];
    std::lock_guard<std::mutex> lock(result.mutex);
    if (!result.data) {
        result.data = Array2D<double>(total_channel_count, sample_count);
    }
    if (result.atoms.empty()) {
        result.atoms.resize(total_channel_count);
    }
    return result;
}

void BufferedWriter::finalize() {
    for (int epoch=0; epoch<real_epoch_count; ++epoch) {
        semaphores[epoch].acquire(total_channel_count);
        BufferedResult& result = results[epoch];
        actual_writer->write(result.data, result.index, result.atoms);
        result.data.reset();
        result.atoms.clear();
    }
    actual_writer->finalize();
}

void BufferedWriter::write(Array2D<double> data, EpochIndex index, const std::vector<std::list<ExportedAtom>> &atoms) {
    const index_t sample_count = data.length();
    BufferedResult& result = get_result_ref(index.epoch_counter, sample_count);
    result.index = EpochIndex{index.epoch_counter, index.epoch_offset, 0};
    for (int c=0; c<data.height(); ++c) {
        std::copy(data[c], data[c] + sample_count, result.data[index.channel_offset + c]);
        result.atoms[index.channel_offset + c] = atoms[c];
    }
    semaphores[index.epoch_counter].release(data.height());
}

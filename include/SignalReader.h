/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_SIGNAL_READER_H
#define EMPI_SIGNAL_READER_H

#include <cstdio>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <vector>
#include "Array.h"
#include "EpochIndex.h"
#include "File.h"

class SignalReader {
public:
    virtual int get_epoch_channel_count() const = 0;

    virtual size_t get_epoch_count() const = 0;

    virtual index_t get_epoch_sample_count() const = 0;

    virtual std::optional<EpochIndex> read(Array2D<double> buffer) = 0;

    virtual ~SignalReader() = default;

protected:
    std::mutex mutex;
};

//------------------------------------------------------------------------------

template<typename REAL>
class SignalReaderBase : public SignalReader {
    const std::vector<int> selected_channels;

protected:
    const int channel_count;
    FileToRead file;

    SignalReaderBase(const char *signal_file_path, int channel_count, std::vector<int> selected_channels)
            : selected_channels(std::move(selected_channels)), channel_count(channel_count), file(signal_file_path) {
        if (this->selected_channels.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("too many selected channels");
        }
    }

    index_t read_into_buffer(Array2D<double>& buffer) {
        std::vector<REAL> sample(channel_count);
        const int selected_channel_count = selected_channels.size();
        if (buffer.height() != selected_channel_count) {
            throw std::logic_error("buffer has invalid height");
        }

        index_t samples_read = 0;
        for (index_t i = 0; i < buffer.length(); ++i) {
            if (fread(sample.data(), sizeof(REAL), channel_count, file.get()) != static_cast<size_t>(channel_count)) {
                break;
            }
            for (int c = 0; c < selected_channel_count; ++c) {
                buffer[c][i] = sample[selected_channels[c] - 1];
            }
            ++samples_read;
        }
        return samples_read;
    }

public:
    int get_epoch_channel_count() const final {
        return selected_channels.size();
    }
};

//------------------------------------------------------------------------------

template<typename REAL>
class SignalReaderForAllEpochs : public SignalReaderBase<REAL> {
    int epochs_read;

protected:
    const index_t epoch_sample_count;

    std::optional<EpochIndex> read_unlocked(Array2D<double> buffer) {
        if (buffer.length() != epoch_sample_count) {
            throw std::logic_error("buffer has invalid length");
        }
        index_t samples_read = SignalReaderBase<REAL>::read_into_buffer(buffer);
        if (!samples_read) {
            return std::nullopt;
        }
        if (samples_read < epoch_sample_count) {
            // the last epoch in a signal file may be truncated
            for (int c = 0; c < buffer.height(); ++c) {
                std::fill(buffer[c] + samples_read, buffer[c] + epoch_sample_count, 0);
            }
        }
        int epoch_offset = epochs_read++;
        return EpochIndex{epoch_offset, epoch_offset};
    }

public:
    SignalReaderForAllEpochs(const char *signal_file_path, int channel_count, std::vector<int> selected_channels, index_t epoch_sample_count)
            : SignalReaderBase<REAL>(signal_file_path, channel_count, std::move(selected_channels)),
              epochs_read(0),
              epoch_sample_count(epoch_sample_count) {
        if (epoch_sample_count <= 0) {
            throw std::logic_error("invalid epoch size");
        }
    }

    size_t get_epoch_count() const override {
        size_t file_size_in_bytes = this->file.get_file_size();
        if (!file_size_in_bytes) {
            return 0;
        }
        // last epoch may be truncated
        return (file_size_in_bytes - 1) / (sizeof(REAL) * epoch_sample_count * this->channel_count) + 1;
    }

    index_t get_epoch_sample_count() const final {
        return epoch_sample_count;
    }

    std::optional<EpochIndex> read(Array2D<double> buffer) override {
        std::lock_guard<std::mutex> lock(this->mutex);
        return read_unlocked(buffer);
    }
};

//------------------------------------------------------------------------------

template<typename REAL>
class SignalReaderForSelectedEpochs : public SignalReaderForAllEpochs<REAL> {
    std::vector<int> epochs;
    std::vector<int>::size_type next_epoch_index;

public:
    SignalReaderForSelectedEpochs(const char *signal_file_path, int channel_count, std::vector<int> selected_channels, index_t epoch_size,
                                  std::vector<int> &&selected_epochs)
            : SignalReaderForAllEpochs<REAL>(signal_file_path, channel_count, std::move(selected_channels), epoch_size), epochs(selected_epochs), next_epoch_index(0) { }

    size_t get_epoch_count() const final {
        return epochs.size();
    }

    std::optional<EpochIndex> read(Array2D<double> buffer) final {
        std::lock_guard<std::mutex> lock(this->mutex);
        if (next_epoch_index >= epochs.size()) {
            return std::nullopt;
        }
        int epoch_counter = next_epoch_index++;
        int epoch_offset = epochs[epoch_counter] - 1;
        size_t seek_offset = sizeof(REAL) * static_cast<size_t>(epoch_offset) * static_cast<size_t>(this->epoch_sample_count) * static_cast<size_t>(this->channel_count);
        if (seek_offset > static_cast<size_t>(std::numeric_limits<long>::max()) || fseek(this->file.get(), seek_offset, SEEK_SET)) {
            throw std::runtime_error("could not seek signal file");
        }
        if (!SignalReaderForAllEpochs<REAL>::read_unlocked(buffer)) {
            return std::nullopt;
        }
        return EpochIndex{epoch_counter, epoch_offset};
    }
};

//------------------------------------------------------------------------------

template<typename REAL>
class SignalReaderForWholeSignal : public SignalReaderBase<REAL> {
    index_t signal_sample_count;

public:
    SignalReaderForWholeSignal(const char *signal_file_path, int channel_count, std::vector<int> selected_channels)
            : SignalReaderBase<REAL>(signal_file_path, channel_count, std::move(selected_channels)) {
        signal_sample_count = this->file.get_file_size() / (sizeof(REAL) * channel_count);
    }

    size_t get_epoch_count() const final {
        return 1;
    }

    index_t get_epoch_sample_count() const final {
        return signal_sample_count;
    }

    std::optional<EpochIndex> read(Array2D<double> buffer) final {
        if (buffer.length() != signal_sample_count) {
            throw std::logic_error("buffer has invalid length");
        }
        std::lock_guard<std::mutex> lock(this->mutex);
        index_t samples_read = SignalReaderBase<REAL>::read_into_buffer(buffer);
        if (!samples_read) {
            return std::nullopt;
        }
        if (samples_read != signal_sample_count) {
            throw std::runtime_error("could not read from file");
        }
        return EpochIndex{0, 0};
    }
};

//------------------------------------------------------------------------------

class SignalReaderSingleChannel : public SignalReader {
    std::shared_ptr<SignalReader> source;
    Array2D<double> epoch;
    std::optional<EpochIndex> last_epoch;

public:
    explicit SignalReaderSingleChannel(std::shared_ptr<SignalReader> source_);

    int get_epoch_channel_count() const final;

    size_t get_epoch_count() const final;

    index_t get_epoch_sample_count() const final;

    std::optional<EpochIndex> read(Array2D<double> buffer) final;
};

//------------------------------------------------------------------------------

#endif //EMPI_SIGNAL_READER_H

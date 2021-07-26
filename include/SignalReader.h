/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_SIGNALREADER_H
#define	EMPI_SIGNALREADER_H

#include <cstdio>
#include <memory>
#include <optional>
#include <vector>

#include "Array.h"
#include "File.h"

struct EpochIndex {
    int epoch;
    int channel_offset;

    EpochIndex(int epoch, int channel_offset = 0) : epoch(epoch), channel_offset(channel_offset) { }
};

//------------------------------------------------------------------------------

class SignalReader {
public:
    virtual int get_epoch_channel_count() const = 0;

    virtual index_t get_epoch_sample_count() const = 0;

    virtual std::optional<EpochIndex> read(Array2D<double> buffer) = 0;

    virtual ~SignalReader() = default;
};

//------------------------------------------------------------------------------

class SignalReaderBase : public SignalReader {
    const std::vector<int> selected_channels;

protected:
    const int channel_count;
    File file;

    SignalReaderBase(const char* signal_file_path, int channel_count, std::vector<int>&& selected_channels);

    index_t read_into_buffer(Array2D<double> buffer);

public:
    int get_epoch_channel_count() const final;
};

//------------------------------------------------------------------------------

class SignalReaderForAllEpochs : public SignalReaderBase {
    int epochs_read;

protected:
	const index_t epoch_sample_count;

public:
	SignalReaderForAllEpochs(const char* signal_file_path, int channel_count, std::vector<int>&& selected_channels, index_t epoch_size);

    index_t get_epoch_sample_count() const final;

    std::optional<EpochIndex> read(Array2D<double> buffer) override;
};

//------------------------------------------------------------------------------

class SignalReaderForSelectedEpochs : public SignalReaderForAllEpochs {
	std::vector<int> epochs;
	std::vector<int>::size_type next_epoch_index;

public:
	SignalReaderForSelectedEpochs(const char* signal_file_path, int channel_count, std::vector<int>&& selected_channels, index_t epoch_size, std::vector<int>&& selected_epochs);

    std::optional<EpochIndex> read(Array2D<double> buffer) final;
};

//------------------------------------------------------------------------------

class SignalReaderForWholeSignal : public SignalReaderBase {
    index_t signal_sample_count;

public:
	SignalReaderForWholeSignal(const char* signal_file_path, int channel_count, std::vector<int>&& selected_channels);

    index_t get_epoch_sample_count() const final;

    std::optional<EpochIndex> read(Array2D<double> buffer) final;
};

//------------------------------------------------------------------------------

class SignalReaderSingleChannel : public SignalReader {
    std::unique_ptr<SignalReader> source;
    Array2D<double> epoch;
    std::optional<EpochIndex> last_epoch;

public:
    explicit SignalReaderSingleChannel(std::unique_ptr<SignalReader> &&source_);

    int get_epoch_channel_count() const final;

    index_t get_epoch_sample_count() const final;

    std::optional<EpochIndex> read(Array2D<double> buffer) final;
};

//------------------------------------------------------------------------------

#endif //EMPI_SIGNALREADER_H

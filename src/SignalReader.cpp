/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <limits>
#include "SignalReader.h"

//------------------------------------------------------------------------------

SignalReaderBase::SignalReaderBase(const char* signal_file_path, int channel_count, std::vector<int>&& selected_channels)
: selected_channels(selected_channels), channel_count(channel_count), file(signal_file_path, "rb") {
	if (this->selected_channels.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
	    throw std::runtime_error("too many selected channels");
	}
}

int SignalReaderBase::get_epoch_channel_count() const {
    return selected_channels.size();
}

index_t SignalReaderBase::read_into_buffer(Array2D<double> buffer) {
	std::vector<float> sample(channel_count);
	const int selected_channel_count = selected_channels.size();
    if (buffer.height() != selected_channel_count) {
        throw std::logic_error("buffer has invalid height");
    }

	index_t samples_read = 0;
	for (index_t i=0; i<buffer.length(); ++i) {
        if (fread(sample.data(), sizeof(float), channel_count, file.get()) != static_cast<size_t>(channel_count)) {
            break;
        }
        for (int c=0; c<selected_channel_count; ++c) {
            buffer[c][i] = sample[selected_channels[c]-1];
        }
        ++samples_read;
	}
	return samples_read;
}

//------------------------------------------------------------------------------

SignalReaderForAllEpochs::SignalReaderForAllEpochs(const char* signal_file_path, int channel_count, std::vector<int>&& selected_channels, index_t epoch_sample_count)
: SignalReaderBase(signal_file_path, channel_count, std::move(selected_channels)), epochs_read(0), epoch_sample_count(epoch_sample_count) {
    if (epoch_sample_count <= 0) {
        throw std::logic_error("invalid epoch size");
    }
}

index_t SignalReaderForAllEpochs::get_epoch_sample_count() const {
    return epoch_sample_count;
}

std::optional<EpochIndex> SignalReaderForAllEpochs::read(Array2D<double> buffer) {
    if (buffer.length() != epoch_sample_count) {
        throw std::logic_error("buffer has invalid length");
    }
    index_t samples_read = read_into_buffer(buffer);
    if (!samples_read) {
        return std::nullopt;
    }
    if (samples_read < epoch_sample_count) {
        // the last epoch in a signal file may be truncated
        for (int c=0; c<buffer.height(); ++c) {
            std::fill(buffer[c]+samples_read, buffer[c]+epoch_sample_count, 0);
        }
    }
    return epochs_read++;
}

//------------------------------------------------------------------------------

SignalReaderForSelectedEpochs::SignalReaderForSelectedEpochs(const char* signal_file_path, int channel_count, std::vector<int>&& selected_channels, index_t epoch_size, std::vector<int>&& selected_epochs)
: SignalReaderForAllEpochs(signal_file_path, channel_count, std::move(selected_channels), epoch_size), epochs(selected_epochs), next_epoch_index(0) { }

std::optional<EpochIndex> SignalReaderForSelectedEpochs::read(Array2D<double> buffer) {
    if (next_epoch_index >= epochs.size()) {
        return std::nullopt;
    }
    int actual_epoch = epochs[next_epoch_index++] - 1;
    size_t seek_offset = sizeof(float) * static_cast<size_t>(actual_epoch) * static_cast<size_t>(epoch_sample_count) * static_cast<size_t>(channel_count);
    if (seek_offset > static_cast<size_t>(std::numeric_limits<long>::max()) || fseek(file.get(), seek_offset, SEEK_SET)) {
		throw std::runtime_error("could not seek signal file");
	}
    if (!SignalReaderForAllEpochs::read(buffer)) {
        return std::nullopt;
    }
    return actual_epoch;
}

//------------------------------------------------------------------------------

SignalReaderForWholeSignal::SignalReaderForWholeSignal(const char* signal_file_path, int channel_count, std::vector<int>&& selected_channels)
: SignalReaderBase(signal_file_path, channel_count, std::move(selected_channels)) {
    // TODO support for files larger than 2 GB on 32-bit systems
    fseek(file.get(), 0, SEEK_END);
    signal_sample_count = ftell(file.get()) / sizeof(float) / channel_count;
    rewind(file.get());
}

index_t SignalReaderForWholeSignal::get_epoch_sample_count() const {
    return signal_sample_count;
}

std::optional<EpochIndex> SignalReaderForWholeSignal::read(Array2D<double> buffer) {
    if (buffer.length() != signal_sample_count) {
        throw std::logic_error("buffer has invalid length");
    }
    index_t samples_read = read_into_buffer(buffer);
    if (!samples_read) {
        return std::nullopt;
    }
    if (samples_read != signal_sample_count) {
        throw std::runtime_error("could not read from file");
    }
    return 0; // first and only epoch
}

//------------------------------------------------------------------------------

SignalReaderSingleChannel::SignalReaderSingleChannel(std::unique_ptr<SignalReader> &&source_)
: source(std::move(source_)), epoch(source->get_epoch_channel_count(), source->get_epoch_sample_count())
{ }

int SignalReaderSingleChannel::get_epoch_channel_count() const {
    return 1;
}

index_t SignalReaderSingleChannel::get_epoch_sample_count() const {
    return source->get_epoch_sample_count();
}

std::optional<EpochIndex> SignalReaderSingleChannel::read(Array2D<double> buffer) {
    if (!last_epoch || ++last_epoch->channel_offset >= epoch.height()) {
        last_epoch = source->read(epoch);
        if (!last_epoch) {
            return std::nullopt;
        }
    }
    std::copy(epoch[last_epoch->channel_offset], epoch[last_epoch->channel_offset] + epoch.length(), buffer[0]);

    return last_epoch;
}

//------------------------------------------------------------------------------

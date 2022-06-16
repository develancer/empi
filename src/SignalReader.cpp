/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "SignalReader.h"

//------------------------------------------------------------------------------

SignalReaderSingleChannel::SignalReaderSingleChannel(std::shared_ptr<SignalReader> source_)
: source(std::move(source_)), epoch(source->get_epoch_channel_count(), source->get_epoch_sample_count())
{ }

size_t SignalReaderSingleChannel::get_epoch_count() const {
    return source->get_epoch_count() * source->get_epoch_channel_count();
}

int SignalReaderSingleChannel::get_epoch_channel_count() const {
    return 1;
}

index_t SignalReaderSingleChannel::get_epoch_sample_count() const {
    return source->get_epoch_sample_count();
}

std::optional<EpochIndex> SignalReaderSingleChannel::read(Array2D<double> buffer) {
    std::lock_guard<std::mutex> lock(this->mutex);
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

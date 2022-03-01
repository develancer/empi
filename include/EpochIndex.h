/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_EPOCH_INDEX_H
#define EMPI_EPOCH_INDEX_H

struct EpochIndex {
    int epoch_counter; // continuous, counting from 0
    int epoch_offset; // epoch offset (not necessarily continuous), the first epoch being 0
    int channel_offset;

    EpochIndex()
    : epoch_counter(0), epoch_offset(0), channel_offset(0) {}

    EpochIndex(int epoch_counter, int epoch_offset, int channel_offset = 0)
    : epoch_counter(epoch_counter), epoch_offset(epoch_offset), channel_offset(channel_offset) {}

    bool operator<(const EpochIndex& other) const {
        if (epoch_counter < other.epoch_counter) {
            return true;
        }
        if (epoch_counter > other.epoch_counter) {
            return false;
        }
        return (channel_offset < other.channel_offset);
    }
};

//------------------------------------------------------------------------------

#endif //EMPI_EPOCH_INDEX_H

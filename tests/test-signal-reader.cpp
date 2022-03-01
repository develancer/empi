/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "SignalReader.h"
#include "Testing.h"

const char *const tmp_name = ".test-signal-reader.tmp";

void prepare_signal(int total_sample_count) {
    FILE *file = fopen(tmp_name, "w+b");
    float data[total_sample_count];
    for (int i = 0; i < total_sample_count; ++i) {
        data[i] = static_cast<float>(i);
    }
    fwrite(data, sizeof(float), total_sample_count, file);
    fclose(file);
}

void test_empty_file_all_epochs() {
    prepare_signal(0);

    std::vector<int> selected_channels = {1, 2};
    SignalReaderForAllEpochs<float> reader(tmp_name, 3, std::move(selected_channels), 8);

    ASSERT_EQUALS(2, reader.get_epoch_channel_count());
    ASSERT_EQUALS(8, reader.get_epoch_sample_count());

    Array2D<double> buffer(2, 8);
    auto result = reader.read(buffer);
    ASSERT(!result);
}

void test_empty_file_whole_signal() {
    prepare_signal(0);

    std::vector<int> selected_channels = {1, 2};
    SignalReaderForWholeSignal<float> reader(tmp_name, 3, std::move(selected_channels));

    ASSERT_EQUALS(2, reader.get_epoch_channel_count());
    ASSERT_EQUALS(0, reader.get_epoch_sample_count());

    Array2D<double> buffer(2, 0);
    auto result = reader.read(buffer);
    ASSERT(!result);
}

void test_all_epochs() {
    prepare_signal(36);

    std::vector<int> selected_channels = {1, 2};
    SignalReaderForAllEpochs<float> reader(tmp_name, 3, std::move(selected_channels), 8);

    ASSERT_EQUALS(2, reader.get_epoch_channel_count());
    ASSERT_EQUALS(8, reader.get_epoch_sample_count());

    Array2D<double> buffer(2, 8);
    auto result = reader.read(buffer);
    ASSERT(result);
    ASSERT_EQUALS(0, result->channel_offset);
    ASSERT_EQUALS(0, result->epoch_offset);
    ASSERT_EQUALS(0, result->epoch_counter);
    for (int i = 0; i < 8; ++i) {
        ASSERT_EQUALS(3 * i, buffer[0][i]);
        ASSERT_EQUALS(3 * i + 1, buffer[1][i]);
    }
    result = reader.read(buffer);
    ASSERT(result);
    ASSERT_EQUALS(0, result->channel_offset);
    ASSERT_EQUALS(1, result->epoch_offset);
    ASSERT_EQUALS(1, result->epoch_counter);
    for (int i = 0; i < 4; ++i) {
        ASSERT_EQUALS(3 * (i + 8), buffer[0][i]);
        ASSERT_EQUALS(3 * (i + 8) + 1, buffer[1][i]);
    }
    for (int i = 4; i < 8; ++i) {
        ASSERT_EQUALS(0, buffer[0][i]);
        ASSERT_EQUALS(0, buffer[1][i]);
    }
    result = reader.read(buffer);
    ASSERT(!result);
}

void test_whole_signal() {
    prepare_signal(36);

    std::vector<int> selected_channels = {1, 2};
    SignalReaderForWholeSignal<float> reader(tmp_name, 3, std::move(selected_channels));

    ASSERT_EQUALS(2, reader.get_epoch_channel_count());
    ASSERT_EQUALS(12, reader.get_epoch_sample_count());

    Array2D<double> buffer(2, 12);
    auto result = reader.read(buffer);
    ASSERT(result);
    ASSERT_EQUALS(0, result->channel_offset);
    ASSERT_EQUALS(0, result->epoch_offset);
    ASSERT_EQUALS(0, result->epoch_counter);
    for (int i = 0; i < 12; ++i) {
        ASSERT_EQUALS(3 * i, buffer[0][i]);
        ASSERT_EQUALS(3 * i + 1, buffer[1][i]);
    }
    result = reader.read(buffer);
    ASSERT(!result);
}

int main() {
    test_empty_file_all_epochs();
    test_empty_file_whole_signal();
    test_all_epochs();
    test_whole_signal();

    remove(tmp_name);
    puts("OK");
}

#include <cassert>
#include <cstdio>

#include "SignalReader.h"

const char* const tmp_name = ".test-signal-reader.tmp";

void prepare_signal(int total_sample_count) {
    FILE* file = fopen(tmp_name, "w+b");
    float data[total_sample_count];
    for (int i=0; i<total_sample_count; ++i) {
        data[i] = static_cast<float>(i);
    }
    fwrite(data, sizeof(float), total_sample_count, file);
    fclose(file);
}

void test_empty_file_all_epochs() {
    prepare_signal(0);

    std::vector<int> selected_channels = {1, 2};
    SignalReaderForAllEpochs reader(tmp_name, 3, std::move(selected_channels), 8);

    assert(reader.get_epoch_channel_count() == 2);
    assert(reader.get_epoch_sample_count() == 8);

    Array2D<double> buffer(2, 8);
    auto result = reader.read(buffer);
    assert(!result);
}

void test_empty_file_whole_signal() {
    prepare_signal(0);

    std::vector<int> selected_channels = {1, 2};
    SignalReaderForWholeSignal reader(tmp_name, 3, std::move(selected_channels));

    assert(reader.get_epoch_channel_count() == 2);
    assert(reader.get_epoch_sample_count() == 0);

    Array2D<double> buffer(2, 0);
    auto result = reader.read(buffer);
    assert(!result);
}

void test_all_epochs() {
    prepare_signal(36);

    std::vector<int> selected_channels = {1, 2};
    SignalReaderForAllEpochs reader(tmp_name, 3, std::move(selected_channels), 8);

    assert(reader.get_epoch_channel_count() == 2);
    assert(reader.get_epoch_sample_count() == 8);

    Array2D<double> buffer(2, 8);
    auto result = reader.read(buffer);
    assert(result);
    assert(result->channel_offset == 0);
    assert(result->epoch == 0);
    for (int i=0; i<8; ++i) {
        assert(buffer[0][i] == 3*i);
        assert(buffer[1][i] == 3*i+1);
    }
    result = reader.read(buffer);
    assert(result);
    assert(result->channel_offset == 0);
    assert(result->epoch == 1);
    for (int i=0; i<4; ++i) {
        assert(buffer[0][i] == 3*(i+8));
        assert(buffer[1][i] == 3*(i+8)+1);
    }
    for (int i=4; i<8; ++i) {
        assert(buffer[0][i] == 0);
        assert(buffer[1][i] == 0);
    }
    result = reader.read(buffer);
    assert(!result);
}

void test_whole_signal() {
    prepare_signal(36);

    std::vector<int> selected_channels = {1, 2};
    SignalReaderForWholeSignal reader(tmp_name, 3, std::move(selected_channels));

    assert(reader.get_epoch_channel_count() == 2);
    assert(reader.get_epoch_sample_count() == 12);

    Array2D<double> buffer(2, 12);
    auto result = reader.read(buffer);
    assert(result);
    assert(result->channel_offset == 0);
    assert(result->epoch == 0);
    for (int i=0; i<12; ++i) {
        assert(buffer[0][i] == 3*i);
        assert(buffer[1][i] == 3*i+1);
    }
    result = reader.read(buffer);
    assert(!result);
}

int main() {
    test_empty_file_all_epochs();
    test_empty_file_whole_signal();
    test_all_epochs();
    test_whole_signal();

    remove(tmp_name);
    puts("OK");
}

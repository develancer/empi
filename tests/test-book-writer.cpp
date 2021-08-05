/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <cstring>
#include "BookWriter.h"
#include "Testing.h"

const char *const tmp_name = ".test-book-writer.tmp";

const char *const expected = "{\n\
\"channel_count\": 2,\n\
\"sampling_frequency_Hz\": 16,\n\
\"segments\": [{\n\
\"sample_count\": 4,\n\
\"segment_length_s\": 0.25,\n\
\"segment_offset_s\": 2,\n\
\"channels\": [{\n\
\"atoms\": [{\n\
\"amplitude\": 20,\n\
\"energy\": 10,\n\
\"envelope\": \"gauss\",\n\
\"f_Hz\": 8,\n\
\"phase\": 3,\n\
\"scale_s\": 3,\n\
\"t0_s\": 0.25,\n\
\"t0_abs_s\": 2.25\n\
}],\n\
\"samples\": [\n\
0,\n\
0,\n\
0,\n\
0\n\
]\n\
},{\n\
\"atoms\": [{\n\
\"amplitude\": 40,\n\
\"energy\": 20,\n\
\"envelope\": \"gauss\",\n\
\"f_Hz\": 8,\n\
\"phase\": 3,\n\
\"scale_s\": 3,\n\
\"t0_s\": 0.25,\n\
\"t0_abs_s\": 2.25\n\
}],\n\
\"samples\": [\n\
0,\n\
0,\n\
0,\n\
0\n\
]\n\
}]\n\
}],\n\
\"segment_count\": 1\n\
}";

ExportedAtom prepare_gabor(double energy) {
    ExportedAtom result(energy);
    result.envelope = "gauss";
    result.amplitude = energy / 8;
    result.frequency = 0.5;
    result.phase = 3.0;
    result.position = 4.0;
    result.scale = 48.0;
    return result;
}

std::vector<std::list<ExportedAtom>> prepare_atoms() {
    std::vector<std::list<ExportedAtom>> atoms(2);
    atoms[0].push_back(prepare_gabor(160.0));
    atoms[1].push_back(prepare_gabor(320.0));
    return atoms;
}

void test_json_writer() {
    Array2D<double> data(2, 4);
    auto atoms = prepare_atoms();
    JsonBookWriter writer(tmp_name);
    writer.write(data, 32, 16, atoms);
    writer.finalize();

    size_t length = strlen(expected);
    char buffer[length];
    FILE *f = fopen(tmp_name, "rb");
    ASSERT_EQUALS(length, fread(buffer, 1, length, f));
    ASSERT(!memcmp(buffer, expected, length));
    fclose(f);
}

int main() {
    test_json_writer();

    remove(tmp_name);
    puts("OK");
}

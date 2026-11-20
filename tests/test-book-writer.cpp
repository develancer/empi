/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <cstring>
#include <string>
#include "BookWriter.h"
#include "Testing.h"

const char *const tmp_name = ".test-book-writer.tmp";
const char *const tmp_db_name = ".test-book-writer.db.tmp";

const BookWriter::Metadata metadata = {{"variant", "MMP1"}, {"signal_info", "a \"quoted\" \\ text"}};

const char *const expected = "{\n\
\"version\": \"" APP_VERSION "\",\n\
\"channel_count\": 2,\n\
\"sampling_frequency_Hz\": 16,\n\
\"variant\": \"MMP1\",\n\
\"signal_info\": \"a \\\"quoted\\\" \\\\ text\",\n\
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
    data.fill(0);
    auto atoms = prepare_atoms();
    JsonBookWriter writer(16, 4, tmp_name, metadata);
    writer.write(data, EpochIndex{0, 8}, atoms);
    writer.finalize();

    size_t length = strlen(expected);
    char buffer[length];
    FILE *f = fopen(tmp_name, "rb");
    ASSERT_EQUALS(length, fread(buffer, 1, length, f));
    ASSERT(!memcmp(buffer, expected, length));
    fclose(f);
}

std::string read_metadata(sqlite3 *db, const char *param) {
    sqlite3_stmt *stmt = nullptr;
    ASSERT_EQUALS(SQLITE_OK, sqlite3_prepare_v2(db, "SELECT value FROM metadata WHERE param = ?", -1, &stmt, nullptr));
    sqlite3_bind_text(stmt, 1, param, -1, SQLITE_STATIC);
    std::string value;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        value = reinterpret_cast<const char *>(sqlite3_column_text(stmt, 0));
    }
    sqlite3_finalize(stmt);
    return value;
}

void test_sqlite_writer() {
    Array2D<double> data(2, 4);
    data.fill(0);
    {
        SQLiteBookWriter writer(16, 4, tmp_db_name, metadata);
        writer.write(data, EpochIndex{0, 8}, prepare_atoms());
        writer.finalize();
    }

    sqlite3 *db = nullptr;
    ASSERT_EQUALS(SQLITE_OK, sqlite3_open_v2(tmp_db_name, &db, SQLITE_OPEN_READONLY, nullptr));
    ASSERT(read_metadata(db, "channel_count") == "2");
    ASSERT(read_metadata(db, "segment_count") == "1");
    ASSERT(read_metadata(db, "variant") == "MMP1");
    ASSERT(read_metadata(db, "signal_info") == "a \"quoted\" \\ text");
    sqlite3_close(db);
}

int main() {
    test_json_writer();
    test_sqlite_writer();

    remove(tmp_name);
    remove(tmp_db_name);
    puts("OK");
}

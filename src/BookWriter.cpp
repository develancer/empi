/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <memory>
#include <stdexcept>
#include <unistd.h>
#include "BookWriter.h"

#if defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN

template<typename X, typename T>
void setBE(X& z, T t) {
    z = static_cast<X>(t);
}

#else

template<typename T>
void setBE(uint32_t &z, T t) {
    const auto x = static_cast<uint32_t>(t);
    z = (x >> 24u) | ((x >> 8u) & 0x0000ff00u) | ((x << 8u) & 0x00ff0000u) | (x << 24u);
}

template<typename T>
void setBE(float &z, T t) {
    const auto x = static_cast<float>(t);
    const auto *p = reinterpret_cast<const uint32_t *>(&x);
    setBE(*reinterpret_cast<uint32_t *>(&z), *p);
}

#endif

//////////////////////////////////////////////////////////////////////////////

BookWriter::BookWriter(double freq_sampling, index_t epoch_sample_count)
        : freq_sampling(freq_sampling), epoch_sample_count(epoch_sample_count) {}

//////////////////////////////////////////////////////////////////////////////

FileBackedBookWriter::FileBackedBookWriter(double freq_sampling, index_t epoch_sample_count, const std::string &path_to_book_file)
: BookWriter(freq_sampling, epoch_sample_count), file(path_to_book_file.c_str(), "wb") {}

void FileBackedBookWriter::finalize() {
    fflush(file.get());
}

//////////////////////////////////////////////////////////////////////////////

void JsonBookWriter::finalize() {
    fprintf(file.get(), "}],\n\"segment_count\": %d\n}", total_segments_written);
    FileBackedBookWriter::finalize();
}

void JsonBookWriter::write(Array2D<double> data, EpochIndex epoch, const std::vector<std::list<ExportedAtom>> &atoms) {
    const int C = data.height();
    const index_t N = data.length();
    const index_t segment_offset = epoch.epoch_offset * epoch_sample_count;
    const double segment_offset_s = static_cast<double>(segment_offset) / freq_sampling;

    if (!total_segments_written) {
        fputs("{\n", file.get());
        fprintf(file.get(), "\"version\": \"%s\",\n", APP_VERSION);
        fprintf(file.get(), "\"channel_count\": %d,\n", C);
        fprintf(file.get(), "\"sampling_frequency_Hz\": %.8lg,\n", freq_sampling);
        fputs("\"segments\": [{\n", file.get());
    } else {
        fputs("},{\n\n", file.get());
    }
    fprintf(file.get(), "\"sample_count\": %ld,\n", N);
    fprintf(file.get(), "\"segment_length_s\": %.15lg,\n", N / freq_sampling);
    fprintf(file.get(), "\"segment_offset_s\": %.15lg,\n", segment_offset_s);
    fputs("\"channels\": [{\n", file.get());

    for (int c = 0; c < C; ++c) {
        if (c > 0) {
            fputs("},{\n", file.get());
        }
        fputs("\"atoms\": [{\n", file.get());

        bool first_atom_written = false;
        for (const ExportedAtom &atom : atoms[c]) {
            if (first_atom_written) {
                fputs("},{\n", file.get());
            }

            fprintf(file.get(), "\"amplitude\": %.8lg,\n", atom.amplitude);
            fprintf(file.get(), "\"energy\": %.8lg,\n", atom.energy / freq_sampling);
            fprintf(file.get(), "\"envelope\": \"%s\",\n", atom.envelope.c_str());
            fprintf(file.get(), "\"f_Hz\": %.8lg,\n", atom.frequency * freq_sampling);
            fprintf(file.get(), "\"phase\": %.8lg,\n", atom.phase);
            fprintf(file.get(), "\"scale_s\": %.8lg,\n", atom.scale / freq_sampling);
            fprintf(file.get(), "\"t0_s\": %.8lg,\n", atom.position / freq_sampling);
            fprintf(file.get(), "\"t0_abs_s\": %.8lg\n", segment_offset_s + atom.position / freq_sampling);

            first_atom_written = true;
        }

        fputs("}],\n\"samples\": [\n", file.get());

        fprintf(file.get(), "%.8lg", data[c][0]);
        for (int i = 1; i < N; ++i) {
            fprintf(file.get(), ",\n%.8lg", data[c][i]);
        }

        fputs("\n]\n", file.get());
    }

    fputs("}]\n", file.get());

    total_segments_written++;
}

//////////////////////////////////////////////////////////////////////////////

SQLiteBookWriter::SQLiteBookWriter(double freq_sampling, index_t epoch_sample_count, const std::string &path_to_book_file)
: BookWriter(freq_sampling, epoch_sample_count) {
    db = create_database(path_to_book_file.c_str());
    execute("PRAGMA journal_mode = OFF");

    create_tables();

    stmt_insert_atom = prepare("INSERT INTO atoms VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
    stmt_insert_metadata = prepare("INSERT INTO metadata VALUES (?, ?)");
    stmt_insert_samples = prepare("INSERT INTO samples VALUES (?, ?, ?)");
    stmt_insert_segment = prepare("INSERT INTO segments VALUES (?, ?, ?, ?)");
}

void SQLiteBookWriter::finalize() {
    insert_metadata("segment_count", total_segments_written, sqlite3_bind_int);
}

void SQLiteBookWriter::write(Array2D<double> data, EpochIndex epoch, const std::vector<std::list<ExportedAtom>> &atoms) {
    const int C = data.height();
    const index_t N = data.length();
    const index_t segment_offset = epoch.epoch_offset * epoch_sample_count;

    execute("BEGIN TRANSACTION");
    if (!total_segments_written) {
        insert_metadata("version", APP_VERSION, [](auto ...params) {
            return sqlite3_bind_text(params..., -1, SQLITE_STATIC);
        });
        insert_metadata("channel_count", C, sqlite3_bind_int);
        insert_metadata("sampling_frequency_Hz", freq_sampling, sqlite3_bind_double);
    }
    insert_segment(total_segments_written, N, N / freq_sampling, segment_offset / freq_sampling);

    std::vector<float> samples(N);
    for (int c = 0; c < C; ++c) {
        int atom_index = 0;
        for (const ExportedAtom &atom : atoms[c]) {
            insert_atom(total_segments_written, c, atom_index++,
                        atom.amplitude,
                        atom.energy / freq_sampling,
                        atom.envelope.c_str(),
                        atom.frequency * freq_sampling,
                        atom.phase,
                        atom.scale / freq_sampling,
                        atom.position / freq_sampling,
                        (segment_offset + atom.position) / freq_sampling
            );
        }
        const double *channel = data[c];
        for (int i = 0; i < N; ++i) {
            setBE(samples[i], channel[i]);
        }
        insert_samples(total_segments_written, c, samples);
    }
    execute("COMMIT TRANSACTION");
    total_segments_written++;
}

std::shared_ptr<sqlite3> SQLiteBookWriter::create_database(const char *path) {
    sqlite3 *db = nullptr;
    unlink(path);
    int result = sqlite3_open_v2(path, &db, SQLITE_OPEN_CREATE | SQLITE_OPEN_READWRITE, nullptr);
    if (result != SQLITE_OK || db == nullptr) {
        throw std::runtime_error("cannot open SQLite file");
    }
    return {db, sqlite3_close_v2};
}

std::shared_ptr<sqlite3_stmt> SQLiteBookWriter::prepare(const char *sql) {
    sqlite3_stmt *stmt = nullptr;
    int result = sqlite3_prepare_v2(db.get(), sql, -1, &stmt, nullptr);
    if (result != SQLITE_OK || stmt == nullptr) {
        throw std::runtime_error("cannot prepare SQLite query");
    }
    return {stmt, sqlite3_finalize};
}

void SQLiteBookWriter::execute(const char *sql) {
    int result = sqlite3_exec(db.get(), sql, nullptr, nullptr, nullptr);
    if (result != SQLITE_OK) {
        throw std::runtime_error("cannot execute SQLite query");
    }
}

void SQLiteBookWriter::create_tables() {
    execute(
            "CREATE TABLE metadata ("
            "param TEXT PRIMARY KEY,"
            "value TEXT NOT NULL"
            ");"
            "CREATE TABLE atoms ("
            "segment_id INTEGER NOT NULL,"
            "channel_id INTEGER NOT NULL,"
            "iteration UNSIGNED INTEGER NOT NULL,"
            "amplitude REAL NOT NULL,"
            "energy REAL NOT NULL,"
            "envelope TEXT NOT NULL,"
            "f_Hz REAL,"
            "phase REAL,"
            "scale_s REAL,"
            "t0_s REAL,"
            "t0_abs_s REAL,"
            "PRIMARY KEY (segment_id, channel_id, iteration)"
            ");"
            "CREATE TABLE segments ("
            "segment_id INTEGER PRIMARY KEY,"
            "sample_count UNSIGNED INTEGER NOT NULL,"
            "segment_length_s REAL NOT NULL,"
            "segment_offset_s REAL NOT NULL"
            ");"
            "CREATE TABLE samples ("
            "segment_id INTEGER NOT NULL,"
            "channel_id INTEGER NOT NULL,"
            "samples_float32 BLOB NOT NULL,"
            "PRIMARY KEY (segment_id, channel_id)"
            ");"
    );
}

void SQLiteBookWriter::insert_atom(int segment_id, int channel_id, int iteration, double amplitude, double energy,
                                   const char *envelope, double f, double phase, double scale, double t0, double t0_abs) {
    if (sqlite3_reset(stmt_insert_atom.get()) != SQLITE_OK
        || sqlite3_bind_int(stmt_insert_atom.get(), 1, segment_id) != SQLITE_OK
        || sqlite3_bind_int(stmt_insert_atom.get(), 2, channel_id) != SQLITE_OK
        || sqlite3_bind_int(stmt_insert_atom.get(), 3, iteration) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_atom.get(), 4, amplitude) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_atom.get(), 5, energy) != SQLITE_OK
        || sqlite3_bind_text(stmt_insert_atom.get(), 6, envelope, -1, SQLITE_STATIC) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_atom.get(), 7, f) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_atom.get(), 8, phase) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_atom.get(), 9, scale) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_atom.get(), 10, t0) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_atom.get(), 11, t0_abs) != SQLITE_OK
        || sqlite3_step(stmt_insert_atom.get()) != SQLITE_DONE
            ) {
        throw std::runtime_error("cannot insert atoms into SQLite file");
    }
}

template<typename B, typename V>
void SQLiteBookWriter::insert_metadata(const char *param, V value, B sqlite3_bind) {
    if (sqlite3_reset(stmt_insert_metadata.get()) != SQLITE_OK
        || sqlite3_bind_text(stmt_insert_metadata.get(), 1, param, -1, SQLITE_STATIC) != SQLITE_OK
        || sqlite3_bind(stmt_insert_metadata.get(), 2, value) != SQLITE_OK
        || sqlite3_step(stmt_insert_metadata.get()) != SQLITE_DONE
            ) {
        throw std::runtime_error("cannot insert metadata into SQLite file");
    }
}

void SQLiteBookWriter::insert_samples(int segment_id, int channel_id, const std::vector<float> &samples) {
    const sqlite3_uint64 segment_bytes = static_cast<sqlite3_uint64>(samples.size()) * sizeof(float);
    if (sqlite3_reset(stmt_insert_samples.get()) != SQLITE_OK
        || sqlite3_bind_int(stmt_insert_samples.get(), 1, segment_id) != SQLITE_OK
        || sqlite3_bind_int(stmt_insert_samples.get(), 2, channel_id) != SQLITE_OK
        || sqlite3_bind_blob64(stmt_insert_samples.get(), 3, samples.data(), segment_bytes, SQLITE_STATIC) != SQLITE_OK
        || sqlite3_step(stmt_insert_samples.get()) != SQLITE_DONE
            ) {
        throw std::runtime_error("cannot insert samples into SQLite file");
    }
}

void SQLiteBookWriter::insert_segment(int segment_id, int sample_count, double segment_length, double segment_offset) {
    if (sqlite3_reset(stmt_insert_segment.get()) != SQLITE_OK
        || sqlite3_bind_int(stmt_insert_segment.get(), 1, segment_id) != SQLITE_OK
        || sqlite3_bind_int(stmt_insert_segment.get(), 2, sample_count) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_segment.get(), 3, segment_length) != SQLITE_OK
        || sqlite3_bind_double(stmt_insert_segment.get(), 4, segment_offset) != SQLITE_OK
        || sqlite3_step(stmt_insert_segment.get()) != SQLITE_DONE
            ) {
        throw std::runtime_error("cannot insert segment into SQLite file");
    }
}

/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BOOK_WRITER_H
#define EMPI_BOOK_WRITER_H

#include <cstddef>
#include <cstdint>
#include <list>
#include <memory>
#include <queue>
#include <vector>
#include "sqlite3.h"
#include "Array.h"
#include "EpochIndex.h"
#include "ExportedAtom.h"
#include "File.h"

/**
 * Base class for all writers for output (decomposition) file formats.
 */
class BookWriter {
protected:
    const double freq_sampling;
    const index_t epoch_sample_count;

    explicit BookWriter(double freq_sampling, index_t epoch_sample_count);

public:
    virtual ~BookWriter() = default;

    /**
     * Finalize writing to output file. This should always be called after all segments have been written.
     */
    virtual void finalize() {}

    /**
     * Write a single multi-channel decomposition and signal segment.
     *
     * @param data reference to multi-channel data of the analysed signal segment
     * @param atoms array of lists of exported atom data, one list for each signal channel
     */
    virtual void write(Array2D<double> data, EpochIndex epoch, const std::vector<std::list<ExportedAtom>> &atoms) = 0;
};

/**
 * Base class for all writers using a regular file handle.
 */
class FileBackedBookWriter : public BookWriter {
protected:
    File file;

    explicit FileBackedBookWriter(double freq_sampling, index_t epoch_sample_count, const std::string &path_to_book_file);

public:
    /**
     * Finalize writing to output file. This should always be called after all segments have been written.
     */
    void finalize() override;
};

/**
 * Writer for decomposition results in JSON file format.
 */
class JsonBookWriter : public FileBackedBookWriter {
    int total_segments_written = 0;

public:
    /**
     * Create a new writer to write results to the given JSON file.
     * If the file exists, it will be truncated. Otherwise, it will be created.
     *
     * @param path_to_book_file output file path
     */
    explicit JsonBookWriter(double freq_sampling, index_t epoch_sample_count, const std::string &path_to_book_file)
    : FileBackedBookWriter(freq_sampling, epoch_sample_count, path_to_book_file) {}

    /**
     * Finalize writing to output file. This should always be called after all segments have been written.
     */
    void finalize() final;

    /**
     * Write a single multi-channel decomposition and signal segment.
     *
     * @param data reference to multi-channel data of the analysed signal segment
     * @param epoch index of epoch, the first epoch being 0
     * @param atoms array of lists of exported atom data, one list for each signal channel
     */
    void write(Array2D<double> data, EpochIndex epoch, const std::vector<std::list<ExportedAtom>> &atoms) final;
};

/**
 * Writer for decomposition results in SQLite file format.
 */
class SQLiteBookWriter : public BookWriter {
    int total_segments_written = 0;

protected:
    std::shared_ptr<sqlite3> db;
    std::shared_ptr<sqlite3_stmt> stmt_insert_atom;
    std::shared_ptr<sqlite3_stmt> stmt_insert_metadata;
    std::shared_ptr<sqlite3_stmt> stmt_insert_samples;
    std::shared_ptr<sqlite3_stmt> stmt_insert_segment;

    static std::shared_ptr<sqlite3> create_database(const char *path);

    std::shared_ptr<sqlite3_stmt> prepare(const char *sql);

    void execute(const char *sql);

    void create_tables();

    void insert_atom(int segment_id, int channel_id, int iteration, double amplitude, double energy,
                     const char *envelope, double f, double phase, double scale, double t0, double t0_abs);

    template<typename B, typename V>
    void insert_metadata(const char *param, V value, B sqlite3_bind);

    void insert_samples(int segment_id, int channel_id, const std::vector<float> &samples);

    void insert_segment(int segment_id, int sample_count, double segment_length, double segment_offset);

public:
    /**
     * Create a new writer to write results to the given SQLite database file.
     * If the file exists, it will be truncated. Otherwise, it will be created.
     *
     * @param path_to_book_file output file path
     */
    explicit SQLiteBookWriter(double freq_sampling, index_t epoch_sample_count, const std::string &path_to_book_file);

    /**
     * Finalize writing to output file. This should always be called after all segments have been written.
     */
    void finalize() final;

    /**
     * Write a single multi-channel decomposition and signal segment.
     *
     * @param data reference to multi-channel data of the analysed signal segment
     * @param epoch index of epoch, the first epoch being 0
     * @param atoms array of lists of exported atom data, one list for each signal channel
     */
    void write(Array2D<double> data, EpochIndex epoch, const std::vector<std::list<ExportedAtom>> &atoms) final;
};

#endif // EMPI_BOOK_WRITER_H

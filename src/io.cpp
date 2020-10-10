/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <limits>
#include <memory>
#include <unistd.h>
#include "io.hpp"

SignalReader::SignalReader(const std::string& pathToSignalFile)
: pathToSignalFile(pathToSignalFile) {
	file = fopen(pathToSignalFile.c_str(), "rb");
	if (!file) throw Exception("couldNotOpenSignalFile");
}

SignalReader::~SignalReader(void) {
	if (file) fclose(file);
}

MultiSignal SignalReader::readEpoch(int samplesToRead) {
	MultiSignal result;
	std::vector<float> sample(channelCount);
	const int C = selectedChannels.size();
	result.channels.resize(C);
	for (int i=0; i<C; ++i) {
		result.channels[i].freqSampling = freqSampling;
	}
	int leftToStore = std::numeric_limits<int>::max();
	while (samplesToRead-- > 0 && fread(sample.data(), sizeof(float), channelCount, file) == static_cast<size_t>(channelCount)) {
		if (--leftToStore < 0) {
			throw Exception("signalFileIsTooLongForThisMachine");
		}
		for (int i=0; i<C; ++i) {
			result.channels[i].samples.push_back(sample[selectedChannels[i]-1]);
		}
	}
	result.segmentOffset = 0;
	return result;
}

void SignalReader::seek(long sampleOffset) {
	if (fseek(file, sizeof(float) * channelCount * sampleOffset, SEEK_SET)) {
		throw Exception("signalSeekFailed");
	}
}

SignalReaderForAllEpochs::SignalReaderForAllEpochs(const std::string& pathToSignalFile, int epochSize)
: SignalReader(pathToSignalFile), epochSize(epochSize), epochsRead(0) { }

MultiSignal SignalReaderForAllEpochs::read() {
	MultiSignal result = readEpoch(epochSize);
	result.segmentOffset = epochSize * epochsRead++;
	int sampleCount = result.channels[0].samples.size();
	if (sampleCount > 0 && sampleCount != epochSize) {
		throw Exception("fileIsTruncated");
	}
	return result;
}

SignalReaderForSelectedEpochs::SignalReaderForSelectedEpochs(const std::string& pathToSignalFile, int epochSize, const std::vector<int>& epochs)
: SignalReaderForAllEpochs(pathToSignalFile, epochSize) {
	for (int epoch : epochs) {
		this->epochs.push(epoch);
	}
}

MultiSignal SignalReaderForSelectedEpochs::read() {
	if (epochs.empty()) {
		return readEpoch(0);
	} else {
		int epoch = epochs.front();
		epochs.pop();
		long segmentOffset = epochSize * (epoch-1);
		seek(segmentOffset);
		MultiSignal result = readEpoch(epochSize);
		result.segmentOffset = segmentOffset;
		int sampleCount = result.channels[0].samples.size();
		if (sampleCount != epochSize) {
			throw Exception("fileIsTruncated");
		}
		return result;
	}
}

SignalReaderForWholeSignal::SignalReaderForWholeSignal(const std::string& pathToSignalFile)
: SignalReader(pathToSignalFile) { }

MultiSignal SignalReaderForWholeSignal::read() {
	return readEpoch(std::numeric_limits<int>::max());
}

//------------------------------------------------------------------------------

BookWriter::BookWriter(const std::string& pathToBookFile)
		: totalSegmentsWritten(0), pathToBookFile(pathToBookFile) { }

//------------------------------------------------------------------------------

FileBackedBookWriter::FileBackedBookWriter(const std::string& pathToBookFile) : BookWriter(pathToBookFile) {
	file = std::shared_ptr<FILE>(fopen(pathToBookFile.c_str(), "wb"), fclose);
	if (!file) throw Exception("couldNotCreateOutputFile");
}

void FileBackedBookWriter::finalize() {
	if (file) {
		fflush(file.get());
	}
}

//------------------------------------------------------------------------------

void LegacyBookWriter::write(const MultiSignal& signal, const MultiChannelResult& result) {
	BookDataHeader headerData;
	BookDataSignalHeader headerSignal;
	BookDataAtomsHeader headerAtoms;
	BookDataAtomHeader headerAtom;

	const int C = signal.channels.size();
	const int N = C ? signal.channels.front().samples.size() : 0;

	if (!totalSegmentsWritten++) {
		BookHeader headerStart;
		// TODO meaningful information
		setBE(headerStart.content.channelCount, C);
		setBE(headerStart.content.dictionarySize, 0);
		setBE(headerStart.content.energyPercent, 100.0);
		setBE(headerStart.content.iterationCount, 1);
		setBE(headerStart.content.pointsPerMicrovolt, 1.0);
		setBE(headerStart.content.freqSampling, C ? signal.channels.front().freqSampling : 0.0);
		fwrite(&headerStart, sizeof headerStart, 1, file.get());
	}

	long headerDataPosition = ftell(file.get());
	setBE(headerData.epochNumber, totalSegmentsWritten);
	setBE(headerData.sampleCount, N);
	fwrite(&headerData, sizeof headerData, 1, file.get());

	std::vector<float> sampleBuffer(N);
	std::vector<float> paramsBuffer;
	for (int c=0; c<C; ++c) {
		setBE(headerSignal.channelNumber, c+1);
		setBE(headerSignal.len, sizeof headerSignal.channelNumber + sizeof(float) * N);
		fwrite(&headerSignal, sizeof headerSignal, 1, file.get());
		for (int i=0; i<N; ++i) {
			setBE(sampleBuffer[i], signal.channels[c].samples[i]);
		}
		fwrite(sampleBuffer.data(), sizeof(float), N, file.get());

		long headerAtomsPosition = ftell(file.get());
		setBE(headerAtoms.channelNumber, c+1);
		fwrite(&headerAtoms, sizeof headerAtoms, 1, file.get());
		for (const Atom& atom : result[c]) {
			const size_t P = atom.params.size();
			setBE(headerAtom.len, sizeof(float) * P);
			setBE(headerAtom.type, atom.type);
			fwrite(&headerAtom, sizeof headerAtom, 1, file.get());
			paramsBuffer.resize(atom.params.size());
			for (size_t p=0; p<P; ++p) {
				setBE(paramsBuffer[p], atom.params[p]);
			}
			fwrite(paramsBuffer.data(), sizeof(float), P, file.get());
		}

		long currentPosition = ftell(file.get());
		fseek(file.get(), headerAtomsPosition, SEEK_SET);
		setBE(headerAtoms.len, currentPosition - headerAtomsPosition - 5);
		fwrite(&headerAtoms, sizeof headerAtoms, 1, file.get());
		fseek(file.get(), currentPosition, SEEK_SET);
	}

	long currentPosition = ftell(file.get());
	fseek(file.get(), headerDataPosition, SEEK_SET);
	setBE(headerData.len, currentPosition - headerDataPosition - 5);
	fwrite(&headerData, sizeof headerData, 1, file.get());
	fseek(file.get(), currentPosition, SEEK_SET);
}

//------------------------------------------------------------------------------

void JsonBookWriter::finalize() {
	fprintf(file.get(), "}],\n\"segment_count\": %d\n}", totalSegmentsWritten);
	FileBackedBookWriter::finalize();
}

void JsonBookWriter::write(const MultiSignal& signal, const MultiChannelResult& result) {
	const int C = signal.channels.size();
	const int N = C ? signal.channels.front().samples.size() : 0;
	const double freqSampling = signal.getFreqSampling();

	if (!totalSegmentsWritten) {
		fputs("{\n", file.get());
		fprintf(file.get(), "\"channel_count\": %d,\n", C);
		fprintf(file.get(), "\"sampling_frequency_Hz\": %.8lg,\n", freqSampling);
		fputs("\"segments\": [{\n", file.get());
	} else {
		fputs("},{\n\n", file.get());
	}
	fprintf(file.get(), "\"sample_count\": %d,\n", N);
	fprintf(file.get(), "\"segment_length_s\": %.15lg,\n", N / freqSampling);
	fprintf(file.get(), "\"segment_offset_s\": %.15lg,\n", signal.segmentOffset / freqSampling);
	fputs("\"channels\": [{\n", file.get());

	for (int c=0; c<C; ++c) {
		if (c > 0) {
			fputs("},{\n", file.get());
		}
		fputs("\"atoms\": [{\n", file.get());

		bool firstAtomWritten = false;
		for (const Atom& atom : result[c]) {
			if (firstAtomWritten) {
				fputs("},{\n", file.get());
			}

			fprintf(file.get(), "\"amplitude\": %.8lg,\n", atom.params[1]);
			fprintf(file.get(), "\"energy\": %.8lg,\n", atom.params[0] * atom.params[0] / freqSampling);
			fprintf(file.get(), "\"envelope\": \"%s\",\n", "gauss");
			fprintf(file.get(), "\"f_Hz\": %.8lg,\n", freqSampling/2 * atom.params[4]);
			fprintf(file.get(), "\"phase\": %.8lg,\n", atom.params[5]);
			fprintf(file.get(), "\"scale_s\": %.8lg,\n", atom.params[3] / freqSampling);
			fprintf(file.get(), "\"t0_s\": %.8lg,\n", atom.params[2] / freqSampling);
			fprintf(file.get(), "\"t0_abs_s\": %.8lg\n", (signal.segmentOffset + atom.params[2]) / freqSampling);

			firstAtomWritten = true;
		}

		fputs("}],\n\"samples\": [\n", file.get());

		fprintf(file.get(), "%.8lg", signal.channels[c].samples[0]);
		for (int i=1; i<N; ++i) {
			fprintf(file.get(), ",\n%.8lg", signal.channels[c].samples[i]);
		}

		fputs("\n]\n", file.get());
	}

	fputs("}]\n", file.get());

	totalSegmentsWritten++;
}

//------------------------------------------------------------------------------

SQLiteBookWriter::SQLiteBookWriter(const std::string& pathToBookFile) : BookWriter(pathToBookFile) {
	db = createDatabase(pathToBookFile.c_str());
	execute("PRAGMA journal_mode = OFF");

	createTables();

	stmtInsertAtom = prepare("INSERT INTO atoms VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
	stmtInsertMetadata = prepare("INSERT INTO metadata VALUES (?, ?)");
	stmtInsertSamples = prepare("INSERT INTO samples VALUES (?, ?, ?)");
	stmtInsertSegment = prepare("INSERT INTO segments VALUES (?, ?, ?, ?)");
}

void SQLiteBookWriter::finalize() {
	insertMetadata("segment_count", totalSegmentsWritten, sqlite3_bind_int);
}

void SQLiteBookWriter::write(const MultiSignal& signal, const MultiChannelResult& result) {
	const int C = signal.channels.size();
	const int N = C ? signal.channels.front().samples.size() : 0;
	const double freqSampling = signal.getFreqSampling();

	if (!totalSegmentsWritten) {
		insertMetadata("channel_count", C, sqlite3_bind_int);
		insertMetadata("sampling_frequency_Hz", freqSampling, sqlite3_bind_double);
	}
	insertSegment(totalSegmentsWritten, N, N / freqSampling, signal.segmentOffset / freqSampling);

	std::vector<float> data(N);
	for (int c=0; c<C; ++c) {
		const int atomCount = result[c].size();
		for (int i=0; i<atomCount; ++i) {
			const Atom& atom = result[c][i];
			insertAtom(totalSegmentsWritten, c, i,
					atom.params[1],
					atom.params[0] * atom.params[0] / freqSampling,
					"gauss",
					freqSampling/2 * atom.params[4],
					atom.params[5],
					atom.params[3] / freqSampling,
					atom.params[2] / freqSampling,
					(signal.segmentOffset + atom.params[2]) / freqSampling
			);
		}
		const std::vector<double>& samples = signal.channels[c].samples;
		for (int i=0; i<N; ++i) {
			setBE(data[i], samples[i]);
		}
		insertSamples(totalSegmentsWritten, c, data);
	}
	totalSegmentsWritten++;
}

std::shared_ptr<sqlite3> SQLiteBookWriter::createDatabase(const char* path) {
	sqlite3* db = nullptr;
	unlink(path);
	int result = sqlite3_open_v2(path, &db, SQLITE_OPEN_CREATE | SQLITE_OPEN_READWRITE, nullptr);
	if (result != SQLITE_OK || db == nullptr) {
		throw Exception("internalSQLiteErrorOpen");
	}
	return std::shared_ptr<sqlite3>(db, sqlite3_close_v2);
}

std::shared_ptr<sqlite3_stmt> SQLiteBookWriter::prepare(const char* sql) {
	sqlite3_stmt* stmt = nullptr;
	int result = sqlite3_prepare_v2(db.get(), sql, -1, &stmt, nullptr);
	if (result != SQLITE_OK || stmt == nullptr) {
		throw Exception("internalSQLiteErrorPrepare");
	}
	return std::shared_ptr<sqlite3_stmt>(stmt, sqlite3_finalize);
}

void SQLiteBookWriter::execute(const char* sql) {
	int result = sqlite3_exec(db.get(), sql, nullptr, nullptr, nullptr);
	if (result != SQLITE_OK) {
		throw Exception("internalSQLiteErrorExecute");
	}
}

void SQLiteBookWriter::createTables() {
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

void SQLiteBookWriter::insertAtom(int segmentID, int channelID, int iteration, double amplitude, double energy, const char* envelope, double f, double phase, double scale, double t0, double t0_abs) {
	if (sqlite3_reset(stmtInsertAtom.get()) != SQLITE_OK
		|| sqlite3_bind_int(stmtInsertAtom.get(), 1, segmentID) != SQLITE_OK
		|| sqlite3_bind_int(stmtInsertAtom.get(), 2, channelID) != SQLITE_OK
		|| sqlite3_bind_int(stmtInsertAtom.get(), 3, iteration) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertAtom.get(), 4, amplitude) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertAtom.get(), 5, energy) != SQLITE_OK
		|| sqlite3_bind_text(stmtInsertAtom.get(), 6, envelope, -1, SQLITE_STATIC) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertAtom.get(), 7, f) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertAtom.get(), 8, phase) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertAtom.get(), 9, scale) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertAtom.get(), 10, t0) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertAtom.get(), 11, t0_abs) != SQLITE_OK
		|| sqlite3_step(stmtInsertAtom.get()) != SQLITE_DONE
	) {
		throw Exception("internalSQLiteErrorInsertAtom");
	}
}

template<typename B, typename V>
void SQLiteBookWriter::insertMetadata(const char* param, V value, B sqlite3_bind) {
	if (sqlite3_reset(stmtInsertMetadata.get()) != SQLITE_OK
		|| sqlite3_bind_text(stmtInsertMetadata.get(), 1, param, -1, SQLITE_STATIC) != SQLITE_OK
		|| sqlite3_bind(stmtInsertMetadata.get(), 2, value) != SQLITE_OK
		|| sqlite3_step(stmtInsertMetadata.get()) != SQLITE_DONE
	) {
		throw Exception("internalSQLiteErrorInsertMetadata");
	}
}

void SQLiteBookWriter::insertSamples(int segmentID, int channelID, const std::vector<float>& samples) {
	if (sqlite3_reset(stmtInsertSamples.get()) != SQLITE_OK
			|| sqlite3_bind_int(stmtInsertSamples.get(), 1, segmentID) != SQLITE_OK
			|| sqlite3_bind_int(stmtInsertSamples.get(), 2, channelID) != SQLITE_OK
			|| sqlite3_bind_blob(stmtInsertSamples.get(), 3, samples.data(), samples.size() * sizeof(float), SQLITE_STATIC) != SQLITE_OK
			|| sqlite3_step(stmtInsertSamples.get()) != SQLITE_DONE
			) {
		throw Exception("internalSQLiteErrorInsertSamples");
	}
}

void SQLiteBookWriter::insertSegment(int segmentID, int sampleCount, double segmentLength, double segmentOffset) {
	if (sqlite3_reset(stmtInsertSegment.get()) != SQLITE_OK
		|| sqlite3_bind_int(stmtInsertSegment.get(), 1, segmentID) != SQLITE_OK
		|| sqlite3_bind_int(stmtInsertSegment.get(), 2, sampleCount) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertSegment.get(), 3, segmentLength) != SQLITE_OK
		|| sqlite3_bind_double(stmtInsertSegment.get(), 4, segmentOffset) != SQLITE_OK
		|| sqlite3_step(stmtInsertSegment.get()) != SQLITE_DONE
	) {
		throw Exception("internalSQLiteErrorInsertSegment");
	}
}

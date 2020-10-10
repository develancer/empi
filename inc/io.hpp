/**********************************************************
 * Piotr T. Różański (c) 2015–2018                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_IO_HPP
#define	EMPI_IO_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <queue>
#include <vector>

#include "base.hpp"
#include "sqlite3.h"

//------------------------------------------------------------------------------

#if defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN

template<typename X, typename T>
void setBE(X& z, T t) {
	z = static_cast<X>(t);
}

#else

template<typename X, typename T>
void setBE(X& z, T t);

template<typename T>
void setBE(uint8_t& z, T t) {
	z = static_cast<uint8_t>(t);
}

template<typename T>
void setBE(uint16_t& z, T t) {
	const uint16_t x = static_cast<uint16_t>(t);
	z = (x >> 8u) | (x << 8u);
}

template<typename T>
void setBE(uint32_t& z, T t) {
	const uint32_t x = static_cast<uint32_t>(t);
	z = (x >> 24u) | ((x >> 8u) & 0x0000ff00u) | ((x << 8u) & 0x00ff0000u) | (x << 24u);
}

template<typename T>
void setBE(float& z, T t) {
	const float x = static_cast<float>(t);
	const uint32_t* p = reinterpret_cast<const uint32_t*>(&x);
	setBE(*reinterpret_cast<uint32_t*>(&z), *p);
}

template<typename T>
T getBE(T t) {
	T z;
	setBE(z, t);
	return z;
}

#endif

//------------------------------------------------------------------------------

#pragma pack(push,1)

struct BookHeaderStart {
	char _version[6] = { 'M','P','v','5','.','0' };
};

struct BookHeaderComment {
	uint8_t _id = 1;
	uint8_t _len[4] = { 0, 0, 0, sizeof _comment };
	char _comment[4] = { 'e', 'm', 'p', 'i' };
};

struct BookHeaderContent {
	uint8_t _id = 2;
	uint8_t _len[4] = { 0, 0, 0, sizeof(BookHeaderContent)-offsetof(BookHeaderContent, _idSignalInfo) };
	uint8_t _idSignalInfo = 5;
	uint8_t _lenSignalInfo = 10;
	float freqSampling;
	float pointsPerMicrovolt;
	uint16_t channelCount;
	uint8_t _idDecompInfo = 6;
	uint8_t _lenDecompInfo = 13;
	float energyPercent;
	uint32_t iterationCount;
	uint32_t dictionarySize;
	char _dictionaryType = 'F';
};

struct BookHeader {
	BookHeaderStart _start;
	BookHeaderComment _comment;
	BookHeaderContent content;
};

struct BookDataHeader {
	uint8_t _id = 7;
	uint32_t len;
	uint16_t epochNumber;
	uint32_t sampleCount;
};

struct BookDataSignalHeader {
	uint8_t _id = 8;
	uint32_t len;
	uint16_t channelNumber;
};

struct BookDataAtomsHeader {
	uint8_t _id = 9;
	uint32_t len;
	uint16_t channelNumber;
};

struct BookDataAtomHeader {
	uint8_t type;
	uint8_t len;
};

#pragma pack(pop)

//------------------------------------------------------------------------------

class SignalReader {
	FILE* file;

protected:
	MultiSignal readEpoch(int samplesToRead);

	void seek(long sampleOffset);

public:
	const std::string pathToSignalFile;
	int channelCount;
	double freqSampling;
	std::vector<int> selectedChannels;

	SignalReader(const std::string& pathToSignalFile);

	virtual ~SignalReader(void);

	virtual MultiSignal read(void) =0;
};

class SignalReaderForAllEpochs : public SignalReader {
protected:
	const int epochSize;
	int epochsRead;

public:
	SignalReaderForAllEpochs(const std::string& pathToSignalFile, int epochSize);

	MultiSignal read(void);
};

class SignalReaderForSelectedEpochs : public SignalReaderForAllEpochs {
	std::queue<int> epochs;

public:
	SignalReaderForSelectedEpochs(const std::string& pathToSignalFile, int epochSize, const std::vector<int>& epochs);

	MultiSignal read(void);
};

class SignalReaderForWholeSignal : public SignalReader {
public:
	SignalReaderForWholeSignal(const std::string& pathToSignalFile);

	MultiSignal read(void);
};

//------------------------------------------------------------------------------

class BookWriter {
protected:
	int totalSegmentsWritten;

public:
	const std::string pathToBookFile;

	explicit BookWriter(const std::string& pathToBookFile);

	virtual ~BookWriter() = default;

	virtual void finalize() { }

	virtual void write(const MultiSignal& signal, const MultiChannelResult& result) =0;
};

//------------------------------------------------------------------------------

class FileBackedBookWriter : public BookWriter {
protected:
	std::shared_ptr<FILE> file;

public:
	explicit FileBackedBookWriter(const std::string& pathToBookFile);

	virtual void finalize();
};

//------------------------------------------------------------------------------

class LegacyBookWriter : public FileBackedBookWriter {
public:
	explicit LegacyBookWriter(const std::string& pathToBookFile) : FileBackedBookWriter(pathToBookFile) { }

	virtual void write(const MultiSignal& signal, const MultiChannelResult& result);
};

//------------------------------------------------------------------------------

class JsonBookWriter : public FileBackedBookWriter {
public:
	explicit JsonBookWriter(const std::string& pathToBookFile) : FileBackedBookWriter(pathToBookFile) { }

	void finalize();

	void write(const MultiSignal& signal, const MultiChannelResult& result);
};

//------------------------------------------------------------------------------

class SQLiteBookWriter : public BookWriter {
protected:
	std::shared_ptr<sqlite3> db;
	std::shared_ptr<sqlite3_stmt> stmtInsertAtom;
	std::shared_ptr<sqlite3_stmt> stmtInsertMetadata;
	std::shared_ptr<sqlite3_stmt> stmtInsertSamples;
	std::shared_ptr<sqlite3_stmt> stmtInsertSegment;

	static std::shared_ptr<sqlite3> createDatabase(const char* path);

	std::shared_ptr<sqlite3_stmt> prepare(const char* sql);

	void execute(const char* sql);

	void createTables();

	void insertAtom(int segmentID, int sampleCount, int iteration, double amplitude, double energy, const char* envelope, double f, double phase, double scale, double t0, double t0_abs);

	template<typename B, typename V>
	void insertMetadata(const char* param, V value, B sqlite3_bind);

	void insertSamples(int segmentID, int channelID, const std::vector<float>& samples);

	void insertSegment(int segmentID, int sampleCount, double segmentLength, double segmentOffset);

public:
	explicit SQLiteBookWriter(const std::string& pathToBookFile);

	void finalize();

	void write(const MultiSignal& signal, const MultiChannelResult& result);
};

#endif	/* EMPI_IO_HPP */

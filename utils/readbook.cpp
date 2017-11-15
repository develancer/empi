#include <cstdlib>
#include <cassert>
#include "io.hpp"

#pragma pack(push,1)

struct ReadAtom {
	BookDataAtomHeader header;
	union {
		struct {
			float modulus, amplitude, t, scale, f, phase;
		};
		float params[6];
	};
};

#pragma pack(pop)

void safe_fread(void *__restrict ptr, size_t size, size_t n, FILE *__restrict stream) {
	size_t read = fread(ptr, size, n, stream);
	assert(read == n);
}

void read(FILE* file) {
	BookHeaderStart headerStart;
	BookHeaderComment headerComment;
	BookHeaderContent headerContent;
	BookDataHeader headerData;
	BookDataSignalHeader headerSignal;
	BookDataAtomsHeader headerAtoms;
	BookDataAtomHeader headerAtom;

	safe_fread(&headerStart, sizeof headerStart, 1, file);
	int firstChar = fgetc(file);
	ungetc(firstChar, file);
	if (firstChar == 1) {
		safe_fread(&headerComment, sizeof headerComment, 1, file);
	}
	safe_fread(&headerContent, sizeof headerContent, 1, file);

	double freqSampling = getBE(headerContent.freqSampling);
	int C = getBE(headerContent.channelCount);

	int epochNumber = 0;
	while (fread(&headerData, sizeof headerData, 1, file) == 1) {
		++epochNumber;
		int epochNumberRead = getBE(headerData.epochNumber);
		int sampleCount = getBE(headerData.sampleCount);
		assert(epochNumberRead == epochNumber);

		for (int c=0; c<C; ++c) {
			safe_fread(&headerSignal, sizeof headerSignal, 1, file);
			fseek(file, getBE(headerSignal.len) - sizeof headerSignal.channelNumber, SEEK_CUR);

			safe_fread(&headerAtoms, sizeof headerAtoms, 1, file);
			int sizeOfAllAtoms = getBE(headerAtoms.len) - sizeof headerAtoms.channelNumber;
			int sizeOfSingleAtom = sizeof headerAtom + 6 * sizeof(float);
			int atomCount = sizeOfAllAtoms / sizeOfSingleAtom;
			for (int a=1; a<=atomCount; ++a) {
				ReadAtom atom;
				safe_fread(&atom, sizeof atom, 1, file);
				for (unsigned i=0; i<(sizeof atom.params)/sizeof(float); ++i) {
					atom.params[i] = getBE(atom.params[i]);
				}
				atom.t += (epochNumber-1) * sampleCount;
				atom.t /= freqSampling;
				atom.scale /= freqSampling;
				atom.f *= freqSampling / 2;
				printf("%2d %3d", c+1, a);
				for (unsigned i=0; i<(sizeof atom.params)/sizeof(float); ++i) {
					printf(" %18.15e", atom.params[i]);
				}
				putchar('\n');
			}
		}
	}
}

int main(int argc, char** argv) {
	if (argc < 2) {
		return 1;
	}
	FILE* file = fopen(argv[1], "rb");
	read(file);
	fclose(file);
}

/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <fstream>
#include <string>
#include <vector>
#include "classes.hpp"
#include "conf.hpp"
#include "gabor.hpp"
#include "io.hpp"
#include "mmp.hpp"

static std::vector<int> parseChannelSpecification(const std::string& string, int channelCount) {
	std::string scString;
	std::stringstream scStream(string);
	std::vector<int> numbers;
	while (std::getline(scStream, scString, ',')) {
		int start, end, status;
		status = sscanf(scString.c_str(), "%d-%d", &start, &end);

		if (status == 2 && start >= 1 && end <= channelCount) {
			for (int i=start; i<=end; ++i) {
				numbers.push_back(i);
			}
		} else if (status == 1 && start >= 1 && start <= channelCount) {
			numbers.push_back(start);
		} else {
			throw Exception("invalidSelectedChannels");
		}
	}
	if (!numbers.size()) {
		throw Exception("noChannelsSelected");
	}
	if (static_cast<intmax_t>(numbers.size()) > static_cast<intmax_t>(std::numeric_limits<int>::max())) {
		throw Exception("tooManySelectedChannels");
	}
	return numbers;
}

static void empi(const char* configFilePath) {
	std::unique_ptr<WorkspaceBuilder> builder;
	std::unique_ptr<Decomposition> decomposition;
	DecompositionSettings settings;
	SignalReader reader;
	BookWriter writer;

	// legacy configuration START /////////////////////////////////////
	Configuration legacyConfiguration;
	legacyConfiguration.parse(configFilePath);

	double energyError = atof(legacyConfiguration.at("energyError").c_str());
	if (energyError <= 0.0 || energyError >= 1.0) {
		throw Exception("invalidEnergyErrorValue");
	}

	builder.reset( new GaborWorkspaceBuilder(energyError) );

	settings.iterationMax = atoi(legacyConfiguration.at("maximalNumberOfIterations").c_str());
	if (settings.iterationMax <= 0) {
		throw Exception("invalidMaximalNumberOfIterations");
	}
	settings.residualEnergy = 1 - 0.01 * atof(legacyConfiguration.at("energyPercent").c_str());
	if (settings.residualEnergy < 0 || settings.residualEnergy >= 1.0) {
		throw Exception("invalidEnergyPercentValue");
	}

	reader.pathToSignalFile = legacyConfiguration.at("nameOfDataFile");
	reader.freqSampling = atof(legacyConfiguration.at("samplingFrequency").c_str());
	if (reader.freqSampling <= 0.0) {
		throw Exception("invalidSamplingFrequency");
	}
	reader.channelCount = atoi(legacyConfiguration.at("numberOfChannels").c_str());
	if (reader.channelCount <= 0) {
		throw Exception("invalidNumberOfChannels");
	}
	reader.selectedChannels = parseChannelSpecification(legacyConfiguration.at("selectedChannels"), reader.channelCount);

	std::string typeOfMP = legacyConfiguration.at("MP");
	std::transform(typeOfMP.begin(), typeOfMP.end(), typeOfMP.begin(), tolower);
	if (reader.selectedChannels.size() == 1 || typeOfMP == "smp") {
		decomposition.reset( new SmpDecomposition );
	} else if (typeOfMP == "mmp1") {
		decomposition.reset( new Mmp1Decomposition );
	} else if (typeOfMP == "mmp2") {
		decomposition.reset( new Mmp2Decomposition );
	} else if (typeOfMP == "mmp3") {
		decomposition.reset( new Mmp3Decomposition );
	} else {
		throw Exception("unsupportedDecompositionType");
	}

	std::string prefix = reader.pathToSignalFile.substr(0, reader.pathToSignalFile.find_last_of('.'));
	writer.pathToBookFile = legacyConfiguration.at("nameOfOutputDirectory")+"/"+prefix+"_"+typeOfMP+".b";

	// legacy configuration END ///////////////////////////////////////

	MultiSignal signal = reader.read();
	std::cout << "START" << '\t' << 1 << '\t' << signal.channels.size() << '\t' << settings.iterationMax << '\t' << 100*(1-settings.residualEnergy) << std::endl;
	MultiChannelResult result = decomposition->compute(settings, *builder, signal);
	writer.write(signal, result);
	puts(" END");
}

static void exception(const char* message) {
	std::cout << "ERROR\t" << message << std::endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char** argv) {
	const char* configFilePath = 0;
	for (int i=1; i<argc; ++i) {
		if (*argv[i] != '-' && *argv[i] != 0) {
			configFilePath = argv[i];
		}
	}
	if (!configFilePath) {
		fprintf(stderr, "USAGE: %s path_to_config_file\n", argv[0]);
		return 1;
	}
	try {
		empi(configFilePath);
	} catch (const Exception& e) {
		exception(e.what());
	} catch (const std::bad_alloc& e) {
		exception("insufficientMemory");
	} catch (...) {
		exception("internalError");
	}
	exit(EXIT_SUCCESS);
}

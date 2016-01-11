/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#include <algorithm>
#include <cstdio>
#include <map>
#include <fstream>
#include <string>
#include "classes.hpp"
#include "gabor.hpp"
#include "io.hpp"
#include "mmp.hpp"

typedef std::map<std::string, std::string> StringDictionary;

StringDictionary parseLegacyConfiguration(const std::string& pathToFile) {
	StringDictionary result;
	std::string line;
	std::ifstream stream(pathToFile);
	while (std::getline(stream, line)) {
		std::string::size_type iname = line.find_first_not_of(" \t");
		if (iname == std::string::npos || line[iname] == '#') {
			continue;
		}
		std::string::size_type ispace = line.find_first_of(" \t", iname+1);
		if (ispace == std::string::npos) {
			continue;
		}
		std::string::size_type ivalue = line.find_first_not_of(" \t", ispace+1);
		if (ivalue == std::string::npos) {
			continue;
		}
		std::string name = line.substr(iname, ispace-iname);
		std::string value = line.substr(ivalue);
		result[name] = value;
	}
	return result;
}

int main(int argc, char** argv) {
	std::unique_ptr<WorkspaceBuilder> builder;
	std::unique_ptr<Decomposition> decomposition;
	DecompositionSettings settings;
	SignalReader reader;
	BookWriter writer;

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

	// legacy configuration START /////////////////////////////////////
	StringDictionary legacyConfiguration = parseLegacyConfiguration(configFilePath);

	double energyError = atof(legacyConfiguration.at("energyError").c_str());
	double scaleMin = 0.1;
	double scaleMax = 5.0;
	builder.reset( new GaborWorkspaceBuilder(energyError, scaleMin, scaleMax) );

	std::string typeOfMP = legacyConfiguration.at("MP");
	std::transform(typeOfMP.begin(), typeOfMP.end(), typeOfMP.begin(), tolower);
	if (typeOfMP == "smp") {
		decomposition.reset( new SmpDecomposition );
//	} else if (typeOfMP == "mmp1") {
//		decomposition.reset( new Mmp1Decomposition );
//	} else if (typeOfMP == "mmp2") {
//		decomposition.reset( new Mmp2Decomposition );
//	} else if (typeOfMP == "mmp3") {
//		decomposition.reset( new Mmp3Decomposition );
	} else {
		throw std::runtime_error("unsupported decomposition type");
	}

	settings.iterationMax = atoi(legacyConfiguration.at("maximalNumberOfIterations").c_str());

	reader.pathToSignalFile = legacyConfiguration.at("nameOfDataFile");
	reader.freqSampling = atof(legacyConfiguration.at("samplingFrequency").c_str());
	reader.channelCount = atoi(legacyConfiguration.at("numberOfChannels").c_str());
	for (int c=1; c<=reader.channelCount; ++c) {
		reader.selectedChannels.push_back(c);
	}

	std::string prefix = reader.pathToSignalFile.substr(0, reader.pathToSignalFile.find_last_of('.'));
	writer.pathToBookFile = legacyConfiguration.at("nameOfOutputDirectory")+"/"+prefix+"_"+typeOfMP+".b";

	// legacy configuration END ///////////////////////////////////////

	MultiSignal signal = reader.read();
	printf("START\t%d\t%d\t%d\t%.2f\n", 1, 1, settings.iterationMax, 99.0);
	printf("EPOCH\t%d\n", 0);
	MultiChannelResult result = decomposition->compute(settings, *builder, signal);
	writer.write(signal, result);
	puts(" END");
}

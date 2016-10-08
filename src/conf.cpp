/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <fstream>
#include "base.hpp"
#include "conf.hpp"

const std::string& Configuration::at(const std::string& key) {
	auto it = values.find(key);
	if (it == values.end()) {
		throw Exception("missingConfigurationEntry:"+key);
	}
	return it->second;
}

bool Configuration::has(const std::string& key) const {
	return values.count(key);
}

void Configuration::parse(const std::string& pathToFile) {
	std::string line;
	std::ifstream stream(pathToFile);
	if (!stream) {
		throw Exception("couldNotOpenConfigFile");
	}
	values.clear();
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
		if (!name.empty() && !value.empty()) {
			values[name] = value;
		}
	}
}

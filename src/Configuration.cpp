/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <fstream>
#include "Configuration.h"

//////////////////////////////////////////////////////////////////////////////

const std::string &Configuration::at(const std::string &key) {
    auto it = values.find(key);
    if (it == values.end()) {
        throw std::runtime_error("missing configuration entry");
    }
    return it->second;
}

bool Configuration::has(const std::string &key) const {
    return values.count(key);
}

void Configuration::parse(const std::string &path_to_file) {
    std::string line;
    std::ifstream stream(path_to_file);
    if (!stream) {
        throw std::runtime_error("could not open config file");
    }
    values.clear();
    while (std::getline(stream, line)) {
        std::string::size_type iname = line.find_first_not_of(" \t");
        if (iname == std::string::npos || line[iname] == '#') {
            continue;
        }
        std::string::size_type ispace = line.find_first_of(" \t", iname + 1);
        if (ispace == std::string::npos) {
            continue;
        }
        std::string::size_type ivalue = line.find_first_not_of(" \t", ispace + 1);
        if (ivalue == std::string::npos) {
            continue;
        }
        std::string name = line.substr(iname, ispace - iname);
        std::string value = line.substr(ivalue);
        if (!name.empty() && !value.empty()) {
            values[name] = value;
        }
    }
}

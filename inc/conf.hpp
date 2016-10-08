/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_CONF_HPP
#define	EMPI_CONF_HPP

#include <map>
#include <string>

class Configuration {
	std::map<std::string, std::string> values;

public:
	void parse(const std::string& pathToFile);
	bool has(const std::string& key) const;
	const std::string& at(const std::string& key);
};

#endif	/* EMPI_CONF_HPP */

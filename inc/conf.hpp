/**********************************************************
 * University of Warsaw, Department of Biomedical Physics *
 * See LICENCE for details.                               *
 **********************************************************/
#ifndef EMPI_CONF_HPP
#define	EMPI_CONF_HPP

#include <map>
#include <string>

class Configuration {
	std::map<std::string, std::string> values;

public:
	void parse(const std::string& pathToFile);
	const std::string& at(const std::string& key);
};

#endif	/* EMPI_CONF_HPP */

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_CONFIGURATION_H
#define EMPI_CONFIGURATION_H

#include <map>
#include <string>

/**
 * Utility class for text-based configuration files.
 * Instances of this class are created empty,
 * and can be populated with a call to parse() method.
 */
class Configuration {
    std::map<std::string, std::string> values;

public:
    /**
     * Parse given configuration file and store results inside this object.
     *
     * @param path_to_file path to configuration file
     */
    void parse(const std::string &path_to_file);

    /**
     * Check whether given key is defined inside the configuration.
     *
     * @param key name of the configuration entry to look for
     * @return true if key is defined, false otherwise
     */
    [[nodiscard]] bool has(const std::string &key) const;

    /**
     * Return configuration value for a given key.
     * If the key is not defined, throw a runtime error.
     *
     * @param key name of the configuration entry to look for
     * @return value associated with the given key
     */
    const std::string &at(const std::string &key);
};

#endif // EMPI_CONFIGURATION_H

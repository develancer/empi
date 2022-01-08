/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_CONFIGURATION_H
#define EMPI_CONFIGURATION_H

#include <list>
#include <string>
#include "CLI11.hpp"
#include "Extractor.h"
#include "Family.h"
#include "OptimizationMode.h"

struct EnvelopeConfiguration {
    std::shared_ptr<Family> family;

    /**
     * Minimum scale (in samples) for the atom, in samples.
     * The special value of 0.0 means that the minimum scale should be automatically calculated.
     */
    double scale_min = 0.0;

    /**
     * Maximum scale (in samples) for the atom, in samples.
     * The special value of 0.0 means that the maximum scale should be automatically calculated.
     */
    double scale_max = 0.0;

    /**
     * Maximum frequency for the atom, relative to the sampling frequency represented as 1.0.
     * The default value of 0.5 represents the Nyquist frequency.
     */
    double freq_max = 0.5;
};

/**
 * Utility class for text-based configuration files.
 * Instances of this class are created empty,
 * and can be populated with a call to parse() method.
 */
struct Configuration {

    /**
     * Number of channels in the input signal.
     */
    int channel_count = 1;

    /**
     * Number of CPU threads to use.
     * The default value is taken from std::thread::hardware_concurrency().
     */
    unsigned cpu_threads;

    /**
     * Energy ε² parameter corresponding to the size of the dictionary.
     * Smaller values correspond to more fine-grained dictionaries.
     */
    double energy_error = 0.05;

    /**
     * Maximum total energy of the residual signal for the iterating procedure to stop.
     * This is represented as the fraction of the initial signal's energy.
     */
    double energy_max_residual = 0.01;

    /**
     * List of configurations for various envelope functions.
     */
    std::map<std::string, EnvelopeConfiguration> envelopes;

    /**
     * Extractor object to be used for multi-channel decomposition.
     */
    Extractor extractor = nullptr;

#ifdef HAVE_CUDA
    /**
     * List of IDs of all GPU devices to use.
     */
    std::vector<int> gpu_devices;
#endif

    /**
     * Whether to treat input data as double precision 64-bit floating point values (true)
     * or as single-precision 32-bit floating point values (false).
     */
    bool input64 = false;

    /**
     * Sampling frequency in hertz, defaults to 1 Hz.
     */
    double freq_sampling = 1.0;

    /**
     * Whether to include delta-type atoms in the dictionary.
     */
    bool include_delta_atoms = false;

    /**
     * Path for the input signal file.
     */
    std::string input_file_path;

    /**
     * Maximum number of iterations to perform.
     * The iterating procedure will stop if either
     * a) max_iterations are reached, or
     * b) energy_max_residual is reached, or
     * c) residual signal's energy reaches zero, or
     * d) no atom matches the residual.
     * The special value of 0 represents the infinite number of iterations,
     * which means that only conditions b-d will be checked.
     */
    int iterations_max = 0; // 0 means infinite

    /**
     * Parameter optimization mode.
     */
    OptimizationMode optimization_mode = OPTIMIZATION_GLOBAL;

    /**
     * Path for the output decomposition file.
     */
    std::string output_file_path;

    /**
     * Range of channels to process, e.g. "1-3,5,8-9", starting from 1.
     * Special value is an empty string representing all channels.
     */
    std::string channel_specs;

    /**
     * Number of samples in each segment.
     */
    int segment_size = 0; // 0 means not specified

    /**
     * Range of signal segments to process, e.g. "1-100,201-300", starting from 1.
     * Special value is an empty string representing all segments.
     */
    std::string segment_specs;

    /**
     * Parse given command-line configuration and store results inside this object.
     *
     * @param error_code will be filled with error code if the returned value is false
     * (in which case the error code may still be 0 in case of using the --help flag)
     * @return whether the program should continue (true) or exit with error_code (false)
     */
    bool parse(int argc, char **argv, int &error_code);

    Configuration();
};

#endif // EMPI_CONFIGURATION_H

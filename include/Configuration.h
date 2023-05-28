/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_CONFIGURATION_H
#define EMPI_CONFIGURATION_H

#include <list>
#include <string>
#include "Extractor.h"
#include "Family.h"
#include "OptimizationMode.h"

struct EnvelopeConfiguration {

    /**
     * Pointer to the shared Family implementation.
     * This field is filled by the parser's final callback.
     */
    std::shared_ptr<Family> family = nullptr;

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
     * Number of CPU threads to be used by a single worker.
     * The default value is taken from std::thread::hardware_concurrency().
     */
    unsigned cpu_threads;

    /**
     * Number of CPU workers to run in parallel.
     * Each worker will be independent and responsible for processing different data segment.
     */
    unsigned cpu_workers = 1;

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

    /**
     * List of IDs of all GPU devices to use.
     */
    std::vector<int> gpu_devices;

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
     * If set to true, part of any atom can exceed the time domain of the signal.
     * Otherwise (default), zero-valued samples are assumed on both ends.
     */
    bool full_atoms_in_signal = false;

    /**
     * Path for the input signal file.
     */
    std::string input_file_path;

    /**
     * Directory in which the residual energy log files should be created.
     */
    std::string residual_log_dir;

    /**
     * Path to create a dictionary structure XML file (optional).
     */
    std::string dictionary_output;

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
     * Maximum number of iterations for local parameter optimization.
     */
    static int optimization_max_iterations;

    /**
     * Parameter optimization mode.
     */
    OptimizationMode optimization_mode = OPTIMIZATION_GLOBAL;

    /**
     * Target simplex size (relative to the initial simplex size of 1) for local parameter optimization.
     */
    static double optimization_target;

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

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <functional>
#include <sstream>
#include <thread>
#include "CLI11.hpp"
#include "Configuration.h"
#include "GaussianFamily.h"

//////////////////////////////////////////////////////////////////////////////

struct EnvelopeTuple {
    const std::string name;
    const std::string option_name;
    const std::string description;
    const std::function<std::shared_ptr<Family>()> builder;
};

//////////////////////////////////////////////////////////////////////////////

std::string check_positive_int(const std::string &input) {
    int value = atoi(input.c_str());
    if (value <= 0) {
        return "Value " + input + " is not a positive integer";
    }
    return {};
}

std::string check_positive_number(const std::string &input) {
    double value = atof(input.c_str());
    if (value <= 0) {
        return "Value " + input + " is not a positive number";
    }
    return {};
}

//////////////////////////////////////////////////////////////////////////////

int Configuration::optimization_max_iterations = 10000;

double Configuration::optimization_target = 1.0e-5;

Configuration::Configuration() : cpu_threads(std::thread::hardware_concurrency()) {}

bool Configuration::parse(int argc, char **argv, int &error_code) {
    CLI::App app{"Enhanced Matching Pursuit Implementation (empi)", "empi"};

    double gauss_half_width = GaussianFamily::DEFAULT_MIN_MAX;
    std::list<EnvelopeTuple> envelope_tuples;
    envelope_tuples.push_back({"gauss", "gabor", "Gaussian", [gauss_half_width]() {
        return std::make_shared<GaussianFamily>(gauss_half_width);
    }});

    bool mmp1 = false, mmp3 = false;

    const std::map<std::string, OptimizationMode> optimization_mode_map = {
            {"none",   OPTIMIZATION_DISABLED},
            {"local",  OPTIMIZATION_LOCAL},
            {"global", OPTIMIZATION_GLOBAL}
    };
    std::string optimization_specs = "global";
    app.add_option("input_file", input_file_path, "Path to the input signal file or input configuration file")->required();
    app.add_option("output_file", output_file_path, "Path for the output file unless configuration file is used")->required();
    app.add_option("-c", channel_count, "Number of channels in the input signal")
            ->check(check_positive_int)->capture_default_str();
    app.add_option("-f", freq_sampling, "Sampling frequency of the input signal in hertz (default: 1 Hz)")
            ->check(check_positive_number);
    app.add_option("-i", iterations_max, "Maximum number of iterations (default: no limit)")
            ->check(check_positive_int);
    app.add_option("-o", optimization_specs, "Parameter optimization mode: none|local|global (default: global)");
    app.add_option("-r", energy_max_residual, "Energy of the residual as a fraction of the total signal energy")
            ->check(check_positive_number)->capture_default_str();
    app.add_option("--channels", channel_specs, "Range of channels to process, e.g. 1-3,5,8-9 (default: all)");
    app.add_option("--cpu-threads", cpu_threads, "Number of CPU threads for each worker")
            ->check(check_positive_int)->capture_default_str();
    app.add_option("--cpu-workers", cpu_workers, "Number of independent CPU workers to run")
            ->check(check_positive_int)->capture_default_str();
    app.add_flag("--delta", include_delta_atoms, "Include delta-type atoms");
    app.add_option("--energy-error", energy_error, "Epsilon-squared parameter corresponding to the dictionary size")
            ->check(check_positive_number)->capture_default_str();
#ifdef HAVE_CUDA
    std::string gpu_specs;
    app.add_option("--dictionary-output", dictionary_output, "Path to create a dictionary structure XML file (default: none)");
    app.add_option("--gpu-id", gpu_specs, "Comma-separated ID list of GPU device(s) to use (default: none)");
#endif
    app.add_flag("--input64", input64, "Read input data as double-precision (64-bit) floating point values (default: read as 32-bit values)");
    app.add_flag("--mmp1", mmp1, "Use multi-variate decomposition with constant phase across channels");
    app.add_flag("--mmp3", mmp3, "Use multi-variate decomposition with variable phase across channels");
    app.add_option("--opt-max-iter", optimization_max_iterations, "Maximum number of iterations for local parameter optimization")
            ->capture_default_str();
    app.add_option("--opt-target", optimization_target, "Target accuracy (relative to the initial dictionary size) for local parameter optimization")
            ->capture_default_str();
    app.add_option("--residual-log-dir", residual_log_dir, "Directory in which residual energy log files should be created (default: none)");
    app.add_option("--segment-size", segment_size, "Number of samples in each segment (default: all samples)");
    app.add_option("--segments", segment_specs, "Range of signal segments, e.g. 1-100,201-300 (default: all)");

    app.get_option("--segments")->needs("--segment-size");
    app.get_option("--mmp1")->excludes("--mmp3");

    std::string all_envelope_option_names;
    for (const auto &e : envelope_tuples) {
        const std::string option_name = "--" + e.option_name;
        if (!all_envelope_option_names.empty()) {
            all_envelope_option_names += ", ";
        }
        all_envelope_option_names += option_name;

        app.add_flag_callback(
                option_name,
                [=]() {
                    envelopes[e.name];
                },
                std::string("Include atoms with ")
                        .append(e.description)
                        .append(" envelope (not needed if any other ")
                        .append(option_name)
                        .append("-* option is given)")
        );

        app.add_option_function<double>(
                option_name + "-freq-max",
                [=](double freq_max) {
                    envelopes[e.name].freq_max = freq_max / freq_sampling;
                    if (envelopes[e.name].freq_max < 0 || envelopes[e.name].freq_max > 0.5) {
                        throw CLI::ValidationError(option_name + "-freq-max", "frequency out of range");
                    }
                },
                "Maximum frequency (in hertz) for " + e.description + " envelope (default: auto)"
        );

        app.add_option_function<double>(
                option_name + "-scale-min",
                [=](double scale_min) {
                    envelopes[e.name].scale_min = scale_min * freq_sampling;
                },
                "Minimum scale (in seconds) for " + e.description + " envelope (default: auto)"
        )->check(check_positive_number);

        app.add_option_function<double>(
                option_name + "-scale-max",
                [=](double scale_max) {
                    envelopes[e.name].scale_max = scale_max * freq_sampling;
                },
                "Maximum scale (in seconds) for " + e.description + " envelope (default: auto)"
        )->check(check_positive_number);

        if (e.name == "gauss") {
            app.add_option(
                    "--gabor-half-width",
                    gauss_half_width,
                    "Half-width of the Gaussian envelope function"
            )->check(check_positive_number)->capture_default_str();
        }
    }

    app.final_callback([&]() {
#ifdef HAVE_CUDA
        if (!gpu_specs.empty()) {
            for (const std::string &id : CLI::detail::split(gpu_specs, ',')) {
                int device_id;
                if (!CLI::detail::lexical_cast(id, device_id) || device_id < 0) {
                    throw CLI::ValidationError("Invalid GPU device ID");
                }
                gpu_devices.push_back(device_id);
            }
        }
#endif
        if (optimization_specs == "none") {
            optimization_mode = OPTIMIZATION_DISABLED;
        } else if (optimization_specs == "local") {
            optimization_mode = OPTIMIZATION_LOCAL;
        } else if (optimization_specs == "global") {
            optimization_mode = OPTIMIZATION_GLOBAL;
        } else {
            throw CLI::ConversionError("invalid value for -o switch");
        }

        if (envelopes.empty() && !include_delta_atoms) {
            throw CLI::ValidationError("You have to select at least one out of " + all_envelope_option_names + " or --delta");
        }

        for (const auto &e : envelope_tuples) {
            auto it = envelopes.find(e.name);
            if (it != envelopes.end()) {
                it->second.family = e.builder();
            }
        }

        if (mmp1) {
            extractor = extractorConstantPhase;
        } else if (mmp3) {
            extractor = extractorVariablePhase;
        } else {
            extractor = extractorSingleChannel;
        }
    });

    try {
        app.parse(argc, argv);
        error_code = 0;
        return true;
    } catch (const CLI::ParseError &e) {
        error_code = app.exit(e);
        return false;
    }
}

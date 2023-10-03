/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <memory>
#include <numeric>
#include <stdexcept>
#include "BlockDictionary.h"
#include "BookWriter.h"
#include "BufferedWriter.h"
#include "Worker.h"
#include "Configuration.h"
#include "Logger.h"
#include "PinnedArray.h"
#include "WorkerLoop.h"
#include "SignalReader.h"
#include "SpectrogramCalculatorFFTW.h"

#ifdef HAVE_CUDA
#include "SpectrogramCalculatorCUDA.h"
#endif

static std::vector<int> parseIntegerSubset(const std::string &string, int max) {
    std::string scString;
    std::stringstream scStream(string);
    std::vector<int> numbers;
    while (std::getline(scStream, scString, ',')) {
        int start, end, status;
        status = sscanf(scString.c_str(), "%d-%d", &start, &end);

        if (status == 2 && start >= 1 && end <= max) {
            for (int i = start; i <= end; ++i) {
                numbers.push_back(i);
            }
        } else if (status == 1 && start >= 1 && start <= max) {
            numbers.push_back(start);
        } else {
            throw std::runtime_error("Configuration contains invalid subset specification");
        }
    }
    if (numbers.empty()) {
        throw std::runtime_error("Configuration contains empty subset specification");
    }
    if (numbers.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Configuration contains too large subset specification");
    }
    return numbers;
}


template<typename REAL>
std::shared_ptr<SignalReader> create_signal_reader(const Configuration &configuration, const std::vector<int> &selected_channels) {
    if (configuration.segment_size > 0) {
        if (!configuration.segment_specs.empty()) {
            std::vector<int> epochs = parseIntegerSubset(configuration.segment_specs, std::numeric_limits<int>::max() / configuration.segment_size);
            return std::make_unique<SignalReaderForSelectedEpochs<REAL>>(
                    configuration.input_file_path.c_str(),
                    configuration.channel_count,
                    selected_channels,
                    configuration.segment_size,
                    std::move(epochs)
            );
        } else {
            return std::make_unique<SignalReaderForAllEpochs<REAL>>(
                    configuration.input_file_path.c_str(),
                    configuration.channel_count,
                    selected_channels,
                    configuration.segment_size
            );
        }
    } else {
        return std::make_unique<SignalReaderForWholeSignal<REAL>>(
                configuration.input_file_path.c_str(),
                configuration.channel_count,
                selected_channels
        );
    }
}

/**
 * Check whether given string ends with given suffix, if compared without case-sensitivity.
 */
static bool ci_ends_with(const std::string &string, const std::string &suffix) {
    const size_t string_size = string.size();
    const size_t suffix_size = suffix.size();
    if (string_size < suffix_size) {
        return false;
    }
    const size_t offset = string_size - suffix_size;
    for (size_t i = 0; i < suffix_size; ++i) {
        if (std::tolower(string[offset + i]) != std::tolower(suffix[i])) {
            return false;
        }
    }
    return true;
}

static int empi(const Configuration &configuration) {
#ifdef HAVE_CUDA
    if (configuration.gpu_devices.empty()) {
        cuda_host_disable();
    }
#endif

    Logger::info("Starting empi");

    double energyError = configuration.energy_error;
    if (energyError <= 0.0 || energyError >= 1.0) {
        throw std::runtime_error("Energy error parameter is invalid");
    }
    double freqSampling = configuration.freq_sampling;
    if (freqSampling <= 0.0) {
        throw std::runtime_error("Sampling frequency is invalid");
    }

    int channel_count = configuration.channel_count;
    if (channel_count <= 0) {
        throw std::runtime_error("Number of channels is invalid");
    }
    std::vector<int> selected_channels;
    if (configuration.channel_specs.empty()) {
        selected_channels.resize(channel_count);
        std::iota(selected_channels.begin(), selected_channels.end(), 1);
    } else {
        selected_channels = parseIntegerSubset(configuration.channel_specs, channel_count);
    }

    if (configuration.extractor != extractorSingleChannel && configuration.channel_count == 1) {
        throw std::runtime_error("MMP does not make sense with only one channel");
    }
    if (configuration.envelopes.empty() && !configuration.include_delta_atoms) {
        throw std::runtime_error("At least one atom type must be selected");
    }

    std::shared_ptr<SignalReader> reader;
    if (configuration.input64) {
        reader = create_signal_reader<double>(configuration, selected_channels);
    } else {
        reader = create_signal_reader<float>(configuration, selected_channels);
    }

    std::unique_ptr<BookWriter> book_writer;
    if (ci_ends_with(configuration.output_file_path, ".json")) {
        book_writer = std::make_unique<JsonBookWriter>(freqSampling, reader->get_epoch_sample_count(), configuration.output_file_path.c_str());
    } else {
        book_writer = std::make_unique<SQLiteBookWriter>(freqSampling, reader->get_epoch_sample_count(), configuration.output_file_path.c_str());
    }
    auto buffered_writer = std::make_shared<BufferedWriter>(reader->get_epoch_channel_count(), reader->get_epoch_count(), std::move(book_writer));

    if (channel_count > 1 && configuration.extractor == extractorSingleChannel) {
        reader = std::make_shared<SignalReaderSingleChannel>(std::move(reader));
    }

    std::list<BlockDictionaryStructure> structures;
    for (const auto &kv : configuration.envelopes) {
        const auto &e = kv.second;

        double scale_min = e.scale_min;
        if (scale_min == 0.0) {
            // defaults to the shortest scale for which the time step does not go below 1 sample
            const double dt_scale = e.family->inv_time_integral(1 - energyError);
            scale_min = 1 / dt_scale;
        }

        double scale_max = e.scale_max;
        if (scale_max == 0.0) {
            // defaults to the epoch length
            scale_max = static_cast<double>(reader->get_epoch_sample_count());
        }
        if (configuration.full_atoms_in_signal) {
            // if --full-atoms-in-signal is set, further restrict the max scale
            double scale_max_in_signal = static_cast<double>(reader->get_epoch_sample_count() - 1)
                                         / (e.family->max_arg() - e.family->min_arg());
            scale_max = std::min(scale_max, scale_max_in_signal);
        }

        structures.emplace_back(e.family, energyError, scale_min, scale_max, e.freq_max);
    }

    FILE* dictionary_xml_handle = nullptr;
    if (!configuration.dictionary_output.empty()) {
        dictionary_xml_handle = fopen(configuration.dictionary_output.c_str(), "w");
        if (!dictionary_xml_handle) {
            throw std::runtime_error("Could not create XML file with dictionary structure");
        }
        fputs("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
              "<dict>\n"
              "<libVersion>0.2</libVersion>\n",
              dictionary_xml_handle);
    }

    std::set<int> transform_sizes;
    for (const auto &structure: structures) {
        transform_sizes.merge(structure.get_transform_sizes());
        if (dictionary_xml_handle) {

            std::string lower_case_name = structure.family->name();
            std::string upper_case_name = lower_case_name;
            for (char &c: upper_case_name) {
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
            fprintf(dictionary_xml_handle, "<blockproperties name=\"%s-WINDOW\">\n", upper_case_name.c_str());
            fprintf(dictionary_xml_handle, "<param name=\"windowtype\" value=\"%s\"/>\n", lower_case_name.c_str());
            fprintf(dictionary_xml_handle, "</blockproperties>\n");
            if (lower_case_name == "gauss") {
                lower_case_name = "gabor";
            }
            for (const auto &bs: structure.block_structures) {
                const double scale_ratio = bs.scale / (bs.envelope_length + 1);
                const double opt = 0.5 * M_1_PI * scale_ratio * scale_ratio;
                fprintf(dictionary_xml_handle, "<block uses=\"%s-WINDOW\">\n", upper_case_name.c_str());
                fprintf(dictionary_xml_handle, "<param name=\"type\" value=\"%s\"/>\n", lower_case_name.c_str());
                fprintf(dictionary_xml_handle, "<param name=\"windowLen\" value=\"%d\"/>\n", bs.envelope_length);
                fprintf(dictionary_xml_handle, "<param name=\"windowShift\" value=\"%g\"/>\n", bs.input_shift);
                fprintf(dictionary_xml_handle, "<param name=\"windowopt\" value=\"%lf\"/>\n", opt);
                fprintf(dictionary_xml_handle, "<param name=\"fftSize\" value=\"%d\"/>\n", bs.transform_size);
                fprintf(dictionary_xml_handle, "</block>\n");
            }
        }
    }

    if (dictionary_xml_handle) {
        fputs("</dict>\n",
              dictionary_xml_handle);
        fflush(dictionary_xml_handle);
        fclose(dictionary_xml_handle);
    }

    std::unique_ptr<SpectrogramCalculatorFFTW> primary_calculator;
    if (!transform_sizes.empty()) {
        primary_calculator = std::make_unique<SpectrogramCalculatorFFTW>(reader->get_epoch_channel_count(), transform_sizes);
    }

    auto progress = std::make_shared<Progress>(reader->get_epoch_count());

    std::list<Thread> worker_threads;
    for (unsigned i=1; i<configuration.cpu_workers; ++i) {
        worker_threads.emplace_back(WorkerLoop(reader, buffered_writer, progress, std::make_unique<SpectrogramCalculatorFFTW>(*primary_calculator), structures, configuration));
    }
    worker_threads.emplace_back(WorkerLoop(reader, buffered_writer, progress, std::move(primary_calculator), structures, configuration));

    buffered_writer->finalize();
    for (auto& thread : worker_threads) {
        thread.join();
    }

    double failed_optimization_percent = BlockAtom::get_failed_optimization_percent();
    if (failed_optimization_percent >= 0.1) {
        Logger::info(
            "%.1f%% of atom optimizations could not fully converge. Unless this number is significant, it should not affect the decomposition. Tweaking --opt-max-iter and --opt-target options might help.",
            failed_optimization_percent
        );
    }

    Logger::info("Decomposition finished successfully");
    return 0;
}

int main(int argc, char **argv) {
    Configuration configuration;
    if (int error_code; !configuration.parse(argc, argv, error_code)) {
        return error_code;
    }

    try {
        return empi(configuration);
    } catch (const std::bad_alloc &e) {
        Logger::error("Could not continue due to insufficient memory");
    } catch (const std::logic_error &e) {
        Logger::internal_error(e.what());
    } catch (const std::exception &e) {
        Logger::error(e.what());
    }
    return EXIT_FAILURE;
}

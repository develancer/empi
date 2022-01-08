/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <thread>
#include "Array.h"
#include "BlockDictionary.h"
#include "BookWriter.h"
#include "Computer.h"
#include "Configuration.h"
#include "DeltaDictionary.h"
#include "Logger.h"
#include "PinnedArray.h"
#include "SignalReader.h"

#ifdef HAVE_CUDA

#include "WorkerCUDA.h"

#endif

#include "WorkerFFTW.h"

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

static double compute_total_energy(Array2D<double> data) {
    double sum2 = 0.0;
    for (int h = 0; h < data.height(); ++h) {
        for (index_t i = 0; i < data.length(); ++i) {
            const double value = data[h][i];
            sum2 += value * value;
        }
    }
    return sum2;
}

template<typename REAL>
std::unique_ptr<SignalReader> create_signal_reader(const Configuration &configuration, const std::vector<int> &selected_channels) {
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
    std::unique_ptr<SignalReader> reader;
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

    if (configuration.input64) {
        reader = create_signal_reader<double>(configuration, selected_channels);
    } else {
        reader = create_signal_reader<float>(configuration, selected_channels);
    }

    const int real_channel_count = reader->get_epoch_channel_count();

    if (channel_count > 1 && configuration.extractor == extractorSingleChannel) {
        reader = std::make_unique<SignalReaderSingleChannel>(std::move(reader));
    }

    PinnedArray2D<double> data(reader->get_epoch_channel_count(), reader->get_epoch_sample_count());
    Array2D<double> initial(real_channel_count, reader->get_epoch_sample_count());
    Computer computer(data, configuration.cpu_threads, configuration.optimization_mode);

    if (configuration.include_delta_atoms) {
        computer.add_dictionary(std::make_unique<DeltaDictionary>(data));
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
            // TODO if the entire envelope needs to fit in epoch, divide by family->max_arg()-family->min_arg()
        }

        structures.emplace_back(e.family, energyError, scale_min, scale_max, e.freq_max);
    }

    std::set<int> transform_sizes;
    for (const auto &structure : structures) {
        transform_sizes.merge(structure.get_transform_sizes());
    }

    if (!transform_sizes.empty()) {
        auto primary_calculator = std::make_unique<WorkerFFTW>(data.height(), transform_sizes);
        for (const auto &structure : structures) {
            computer.add_dictionary(std::make_unique<BlockDictionary>(structure, data, configuration.extractor, *primary_calculator));
        }

        std::list<ProtoRequest> proto_requests = computer.get_proto_requests();

#ifdef HAVE_CUDA
        for (int device : configuration.gpu_devices) {
            computer.add_calculator(std::make_unique<WorkerCUDA>(data.height(), proto_requests, device), true);
        }
#endif
        for (unsigned i = 1; i < configuration.cpu_threads; ++i) {
            computer.add_calculator(std::make_unique<WorkerFFTW>(*primary_calculator));
        }
        computer.add_calculator(std::move(primary_calculator));
    }

    std::unique_ptr<BookWriter> book_writer;
    if (ci_ends_with(configuration.output_file_path, ".json")) {
        book_writer = std::make_unique<JsonBookWriter>(configuration.output_file_path.c_str());
    } else {
        book_writer = std::make_unique<SQLiteBookWriter>(configuration.output_file_path.c_str());
    }

    std::vector<std::list<ExportedAtom>> atoms(real_channel_count);

    const size_t epochs_all = reader->get_epoch_count();
    size_t epochs_processed = 0;
    int old_progress = -1;

    int iterations_max = configuration.iterations_max ? configuration.iterations_max : std::numeric_limits<int>::max();
    while (auto epoch_index = reader->read(data)) {
        for (int c = 0; c < data.height(); ++c) {
            std::copy(data[c], data[c] + data.length(), initial[epoch_index->channel_offset + c]);
        }
        computer.reset();
        const double initial_energy = compute_total_energy(data);
        if (!epoch_index->channel_offset) {
            Logger::info("Starting segment #%d", epoch_index->epoch + 1);
        }
        for (int i = 0; i < iterations_max; ++i) {
            auto atom = computer.get_next_atom();
            if (!atom) {
                break;
            }
            std::list<ExportedAtom> *atoms_per_channel = atoms.data() + epoch_index->channel_offset;
            atom->export_atom(atoms_per_channel);

            double residual_energy = compute_total_energy(data);
            double energy_progress = std::min(1.0, std::log(residual_energy / initial_energy) / std::log(configuration.energy_max_residual));
            double this_epoch_progress = std::max(
                    energy_progress,
                    static_cast<double>(i + 1) / static_cast<double>(iterations_max)
            );
            int progress = Types::floor<int>(100.0 * (static_cast<double>(epochs_processed) + this_epoch_progress) / static_cast<double>(epochs_all));
            if (progress != old_progress) {
                std::cout << progress << "% completed... (segment #" << (epoch_index->epoch + 1);
                if (data.height() != real_channel_count) {
                    std::cout << " channel #" << selected_channels[epoch_index->channel_offset];
                }
                std::cout << " " << Types::floor<int>(100.0 * this_epoch_progress) << "% completed with " << (i + 1) << " atoms)" << std::endl;
                old_progress = progress;
            }
            if (energy_progress == 1.0) {
                break;
            }
        }

        if (epoch_index->channel_offset + data.height() == real_channel_count) {
            Logger::info("Finished segment #%d, writing to file", epoch_index->epoch + 1);
            book_writer->write(initial, epoch_index->epoch, freqSampling, atoms);
            Logger::info("Segment #%d written to file", epoch_index->epoch + 1);
            for (auto &list : atoms) {
                list.clear();
            }
        }
        ++epochs_processed;
    }

    book_writer->finalize();
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

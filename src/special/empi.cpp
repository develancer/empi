/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <thread>
#include "Array.h"
#include "BlockDictionary.h"
#include "BookWriter.h"
#include "Computer.h"
#include "Configuration.h"
#include "Family.h"
#include "PinnedArray.h"
#include "SignalReader.h"
#ifdef HAVE_CUDA
    #include "WorkerCUDA.h"
#endif
#include "WorkerFFTW.h"

struct DecompositionSettings {
    int iterationMax;
    double residualEnergy;
};

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
            throw std::runtime_error("invalid subset specification");
        }
    }
    if (!numbers.size()) {
        throw std::runtime_error("empty subset specification");
    }
    if (numbers.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("too large subset specification");
    }
    return numbers;
}

int compute_window_length(double scale, double energy_error, std::shared_ptr<Family> family) {
    const double min_window_length = scale / family->inv_freq_integral(1 - energy_error);
    index_t envelope_length = family->size_for_values(0.0, scale, nullptr);

    int window_length = 2;
    while (static_cast<double>(window_length) < min_window_length || static_cast<index_t>(window_length) < envelope_length) {
        window_length *= 2;
        if (window_length < 0) {
            // integer overflow
            throw std::runtime_error("requested scale is too large");
        }
    }

    return window_length;
}

std::unique_ptr<BlockDictionary>
create_block_dictionary(double scale_min, double scale_max, double energy_error, PinnedArray2D<double> data, std::shared_ptr<Family> family,
                        Extractor extractor, SpectrumCalculator &calculator) {
    if (scale_min <= 0 || !(scale_min <= scale_max)) { // NOLINT no, this should not be simplified due to NAN and INF
        throw std::runtime_error("requested scale range is invalid");
    }
    if (!(energy_error > 0 && energy_error < 1)) {
        throw std::logic_error("requested energy error is invalid");
    }
    auto dictionary = std::make_unique<BlockDictionary>(data, family);

    double dl = family->inv_scale_integral(1 - energy_error);
    const double l_min = std::log(scale_min);
    const double l_max = std::log(scale_max);
    int il_max = Types::ceil<int>((l_max - l_min) / dl);

    const double dt_scale = family->inv_time_integral(1 - energy_error);
    for (int il = il_max; il >= 0; --il) {
        const double scale = std::exp(l_min + (l_max - l_min) * il / il_max);
        int input_shift = Types::floor<int>(dt_scale * scale + 1.0e-12);
        if (input_shift < 1) {
            throw std::runtime_error("requested minimum atom scale is too small");
        }

        // TODO check for math error or overflow
        int window_length = compute_window_length(scale, energy_error, family);
        int output_bins = window_length / 2 + 1; // TODO only part of spectrum

        dictionary->add_block(scale, window_length, input_shift, output_bins, extractor, calculator);
    }

    return dictionary;
}

std::set<int> estimate_window_lengths(double scale_min, double scale_max, double energy_error, std::shared_ptr<Family> family) {
    std::set<int> window_lengths;

    double dl = family->inv_scale_integral(1 - energy_error);
    const double l_min = std::log(scale_min);
    const double l_max = std::log(scale_max);
    long il_max = std::lrint((l_max - l_min) / dl + 0.5);

    for (int il = il_max; il >= 0; --il) {
        const double scale = std::exp(l_min + (l_max - l_min) * il / il_max);
        int window_length = compute_window_length(scale, energy_error, family);
        window_lengths.insert(window_length);
    }

    return window_lengths;
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

static void empi(const char *configFilePath, bool json_output_mode) {
//	std::unique_ptr<Decomposition> decomposition;
    std::unique_ptr<SignalReader> reader;
    DecompositionSettings settings;

    // legacy configuration START /////////////////////////////////////
    Configuration legacyConfiguration;
    legacyConfiguration.parse(configFilePath);

    double energyError = atof(legacyConfiguration.at("energyError").c_str());
    if (energyError <= 0.0 || energyError >= 1.0) {
        throw std::runtime_error("invalidEnergyErrorValue");
    }
    double freqSampling = atof(legacyConfiguration.at("samplingFrequency").c_str());
    if (freqSampling <= 0.0) {
        throw std::runtime_error("invalidSamplingFrequency");
    }

    double scaleMin = 0.0, scaleMax = 0.0;
    if (legacyConfiguration.has("minAtomScale")) {
        scaleMin = atof(legacyConfiguration.at("minAtomScale").c_str()); // NOLINT atof is fine here
        if (!std::isnormal(scaleMin) || scaleMin < 0) {
            throw std::runtime_error("invalid minAtomScale");
        }
        scaleMin *= freqSampling;
    }
    if (legacyConfiguration.has("maxAtomScale")) {
        scaleMax = atof(legacyConfiguration.at("maxAtomScale").c_str()); // NOLINT atof is fine here
        if (!std::isnormal(scaleMax) || scaleMax < 0) {
            throw std::runtime_error("invalid maxAtomScale");
        }
        scaleMax *= freqSampling;
    }

    double freqMax = legacyConfiguration.has("maxAtomFrequency")
                     ? atof(legacyConfiguration.at("maxAtomFrequency").c_str()) / freqSampling : INFINITY;

    if (legacyConfiguration.has("maximalNumberOfIterations")) {
        settings.iterationMax = atoi(legacyConfiguration.at("maximalNumberOfIterations").c_str());
        if (settings.iterationMax <= 0) {
            throw std::runtime_error("invalidMaximalNumberOfIterations");
        }
    } else {
        settings.iterationMax = std::numeric_limits<int>::max();
    }
    if (legacyConfiguration.has("energyPercent")) {
        settings.residualEnergy = 1 - 0.01 * atof(legacyConfiguration.at("energyPercent").c_str());
        if (settings.residualEnergy < 0 || settings.residualEnergy >= 1.0) {
            throw std::runtime_error("invalidEnergyPercentValue");
        }
    } else {
        settings.residualEnergy = 0.0;
    }

    int channel_count = atoi(legacyConfiguration.at("numberOfChannels").c_str());
    if (channel_count <= 0) {
        throw std::runtime_error("invalidNumberOfChannels");
    }
    std::vector<int> selected_channels = parseIntegerSubset(legacyConfiguration.at("selectedChannels"), channel_count);

    const std::string &pathToSignalFile = legacyConfiguration.at("nameOfDataFile");
    if (legacyConfiguration.has("numberOfSamplesInEpoch")) {
        int epochSize = atoi(legacyConfiguration.at("numberOfSamplesInEpoch").c_str());
        if (epochSize <= 0) {
            throw std::runtime_error("invalidNumberOfSamplesInEpoch");
        }
        if (legacyConfiguration.has("selectedEpochs")) {
            std::vector<int> epochs = parseIntegerSubset(legacyConfiguration.at("selectedEpochs"), std::numeric_limits<int>::max() / epochSize);
            if (epochs.empty()) {
                throw std::runtime_error("noSelectedEpochs");
            }
            reader.reset(new SignalReaderForSelectedEpochs(pathToSignalFile.c_str(), channel_count, std::move(selected_channels), epochSize,
                                                           std::move(epochs)));
        } else {
            reader.reset(new SignalReaderForAllEpochs(pathToSignalFile.c_str(), channel_count, std::move(selected_channels), epochSize));
        }
    } else {
        if (legacyConfiguration.has("selectedEpochs")) {
            // entry for numberOfSamplesInEpoch is missing
            legacyConfiguration.at("numberOfSamplesInEpoch");
        }
        reader.reset(new SignalReaderForWholeSignal(pathToSignalFile.c_str(), channel_count, std::move(selected_channels)));
    }

    const int real_channel_count = reader->get_epoch_channel_count();

    Extractor extractor;
    std::string typeOfMP = legacyConfiguration.at("MP");
    std::transform(typeOfMP.begin(), typeOfMP.end(), typeOfMP.begin(), tolower);
    if (typeOfMP == "smp") {
        reader = std::make_unique<SignalReaderSingleChannel>(std::move(reader));
        extractor = extractorVariablePhase;
//  } else if (typeOfMP == "mmp1") {
//      extractor = extractorConstantPhase;
    } else if (typeOfMP == "mmp3") {
        extractor = extractorVariablePhase;
    } else {
        throw std::runtime_error("unsupportedDecompositionType");
    }

    PinnedArray2D<double> data(reader->get_epoch_channel_count(), reader->get_epoch_sample_count());
    Array2D<double> initial(real_channel_count, reader->get_epoch_sample_count());
    Computer computer(data);

    if (data.height() == 1) {
        extractor = extractorSingleChannel;
    }

    auto family = std::make_shared<GaussianFamily>();

    if (scaleMin == 0.0) {
        // defaults to the shortest scale for which the time step does not go below 1 sample
        const double dt_scale = family->inv_time_integral(1 - energyError);
        scaleMin = 1 / dt_scale;
    }
    if (scaleMax == 0.0) {
        // defaults to the epoch length
        scaleMax = static_cast<double>(reader->get_epoch_sample_count());
        // TODO if the entire envelope needs to fit in epoch, divide by family->max_arg()-family->min_arg()
    }

    const std::set<int> window_lengths = estimate_window_lengths(scaleMin, scaleMax, energyError, family);
    auto primary_calculator = std::make_unique<WorkerFFTW>(data.height(), window_lengths);

    auto dictionary = create_block_dictionary(scaleMin, scaleMax, energyError, data, family, extractor, *primary_calculator);
    size_t total_atom_count = dictionary->get_atom_count();
    computer.add_dictionary(std::move(dictionary));
    std::list<ProtoRequest> proto_requests = computer.get_proto_requests();

    #ifdef HAVE_CUDA
    computer.add_calculator(std::make_unique<WorkerCUDA>(data.height(), proto_requests), true);
    #endif
    for (unsigned i = 1; i < std::thread::hardware_concurrency(); ++i) {
        computer.add_calculator(std::make_unique<WorkerFFTW>(*primary_calculator));
    }
    computer.add_calculator(std::move(primary_calculator));

    // extracting base name of input file
    std::string prefix = pathToSignalFile;
    std::string::size_type lastSlash = prefix.find_last_of('/');
    if (lastSlash != std::string::npos) {
        prefix = prefix.substr(lastSlash + 1);
    }
    prefix = prefix.substr(0, prefix.find_last_of('.'));
    typeOfMP.resize(3);

    std::string pathToBookFile = legacyConfiguration.at("nameOfOutputDirectory") + "/" + prefix + "_" + typeOfMP;

    std::unique_ptr<BookWriter> book_writer;
    if (json_output_mode) {
        book_writer = std::make_unique<JsonBookWriter>((pathToBookFile + ".json").c_str());
    } else {
        book_writer = std::make_unique<SQLiteBookWriter>((pathToBookFile + ".db").c_str());
    }

    std::vector<std::list<ExportedAtom>> atoms(real_channel_count);

    int epochs_processed = 0;
    while (auto epoch_index = reader->read(data)) {
        for (int c = 0; c < data.height(); ++c) {
            std::copy(data[c], data[c] + data.length(), initial[epoch_index->channel_offset + c]);
        }
        computer.reset();
        const double initial_energy = compute_total_energy(data);
        if (!epoch_index->channel_offset) {
            // first channel of single-channel decomposition OR multi-channel decomposition
            std::cout << "START" << '\t' << ++epochs_processed << '\t' << real_channel_count << '\t' << settings.iterationMax << '\t'
                      << 100 * (1 - settings.residualEnergy) << std::endl;
        }
        if (data.height() == 1) {
            // single-channel decomposition
            std::cout << "CHANNEL" << '\t' << epoch_index->channel_offset + 1 << std::endl;
        }

        for (int i = 0; i < settings.iterationMax; ++i) {
            auto atom = computer.get_next_atom();
            atom->export_atom(atoms.data() + epoch_index->channel_offset);

            double residual_energy = compute_total_energy(data);
            double energy_progress = 100 * (1.0 - residual_energy / initial_energy);
            double iteration_progress = 100.0 * i / settings.iterationMax;
            std::cout << "ATOM" << '\t' << i << '\t' << total_atom_count
                      << '\t' << energy_progress << '\t' << std::max(energy_progress, iteration_progress) << std::endl;
            if (residual_energy < settings.residualEnergy * initial_energy) {
                break;
            }
        }

        if (epoch_index->channel_offset + data.height() == real_channel_count) {
            book_writer->write(initial, epoch_index->epoch, freqSampling, atoms);
            for (auto &list : atoms) {
                list.clear();
            }
        }
    }

    book_writer->finalize();
    std::cout << " END" << std::endl;
}

static void exception(const char *message) {
    std::cout << "ERROR\t" << message << std::endl;
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
    const char *config_file_path = nullptr;
    bool json_output_mode = false;
    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        if (arg[0] && arg[0] != '-') {
            config_file_path = arg;
        } else if (arg[0] == '-' && arg[1] == 'j' && !arg[2]) {
            json_output_mode = true;
        }
    }
    if (!config_file_path) {
        fprintf(stderr, "USAGE: %s [-j] path_to_config_file\n", argv[0]);
        return 1;
    }
    try {
        empi(config_file_path, json_output_mode);
    } catch (const std::bad_alloc &e) {
        exception("insufficient memory");
    } catch (const std::exception &e) {
        exception(e.what());
    }
    exit(EXIT_SUCCESS);
}
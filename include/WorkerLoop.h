/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_WORKER_LOOP_H
#define EMPI_WORKER_LOOP_H

#include <memory>
#include <list>
#include "BlockDictionary.h"
#include "BlockDictionaryStructure.h"
#include "BufferedWriter.h"
#include "Worker.h"
#include "Configuration.h"
#include "DeltaDictionary.h"
#include "Logger.h"
#include "PinnedArray.h"
#include "Progress.h"
#include "SignalReader.h"
#include "SpectrogramCalculatorFFTW.h"

#ifdef HAVE_CUDA

#include "SpectrogramCalculatorCUDA.h"

#endif

/**
 * Runnable object used with Worker instance.
 */
class WorkerLoop {
    std::shared_ptr<SignalReader> signal_reader;
    std::shared_ptr<BookWriter> book_writer;
    std::shared_ptr<Progress> progress;
    PinnedArray2D<double> data;
    Array2D<double> initial;
    Worker computer;
    const int iterations_max;
    const double energy_max_residual;
    const std::string residual_log_dir;

public:
    WorkerLoop(std::shared_ptr<SignalReader> signal_reader_, std::shared_ptr<BufferedWriter> book_writer_,
               std::shared_ptr<Progress> progress_,
               std::unique_ptr<SpectrogramCalculatorFFTW> primary_calculator,
               const std::list<BlockDictionaryStructure> &structures,
               const Configuration &configuration)
            : signal_reader(std::move(signal_reader_)),
              book_writer(std::move(book_writer_)),
              progress(std::move(progress_)),
              data(signal_reader->get_epoch_channel_count(), signal_reader->get_epoch_sample_count()),
              initial(signal_reader->get_epoch_channel_count(), signal_reader->get_epoch_sample_count()),
              computer(data, configuration.cpu_threads, configuration.optimization_mode),
              iterations_max(
                      configuration.iterations_max ? configuration.iterations_max : std::numeric_limits<int>::max()),
              energy_max_residual(configuration.energy_max_residual),
              residual_log_dir(configuration.residual_log_dir) {
        // TODO check cpu_threads > 0
        if (configuration.include_delta_atoms) {
            computer.add_dictionary(std::make_unique<DeltaDictionary>(data));
        }
        bool allow_overstep = !configuration.full_atoms_in_signal;
        for (const auto &structure: structures) {
            computer.add_dictionary(
                    std::make_unique<BlockDictionary>(structure, data, configuration.extractor, *primary_calculator, allow_overstep));
        }

        for (unsigned i = 1; i < configuration.cpu_threads; ++i) {
            computer.add_calculator(std::make_unique<SpectrogramCalculatorFFTW>(*primary_calculator), false);
        }
        computer.add_calculator(std::move(primary_calculator), false);

#ifdef HAVE_CUDA
        std::list<ProtoRequest> proto_requests = computer.get_proto_requests();
        for (int device: configuration.gpu_devices) {
            computer.add_calculator(std::make_unique<SpectrogramCalculatorCUDA>(data.height(), proto_requests, device),
                                    true);
        }
#endif
    }

    WorkerLoop(WorkerLoop&&) = default;
    WorkerLoop(const WorkerLoop&) = delete;
    void operator=(const WorkerLoop&) = delete;

    void operator()() {
        std::vector<std::list<ExportedAtom>> atoms(signal_reader->get_epoch_channel_count());

        size_t epochs_processed = 0;
        while (auto epoch_index = signal_reader->read(data)) {
            FILE* residual_log = nullptr;
            if (!residual_log_dir.empty()) {
                const std::string residual_output_path = residual_log_dir + "/"
                                                         + std::to_string(epoch_index->epoch_counter) + "-"
                                                         + std::to_string(epoch_index->channel_offset) + ".log";
                residual_log = fopen(residual_output_path.c_str(), "wt");
            }
            for (int c = 0; c < data.height(); ++c) {
                std::copy(data[c], data[c] + data.length(), initial[c]);
            }
            computer.reset();
            const double initial_energy = compute_total_energy(data);
            if (residual_log) {
                fprintf(residual_log, "%.6lf\n", initial_energy);
            }
            progress->epoch_started(epoch_index.value());
            for (int i = 0; i < iterations_max; ++i) {
                auto atom = computer.get_next_atom();
                if (!atom) {
                    break;
                }
                atom->export_atom(atoms.data());

                double residual_energy = compute_total_energy(data);
                if (residual_log) {
                    fprintf(residual_log, "%.6lf\n", residual_energy);
                }
                double energy_progress = std::min(1.0, std::log(residual_energy / initial_energy) /
                                                       std::log(energy_max_residual));
                double this_epoch_progress = std::max(
                        energy_progress,
                        static_cast<double>(i + 1) / static_cast<double>(iterations_max)
                );
                if (this_epoch_progress == 1.0) {
                    break;
                }
                progress->epoch_progress(epoch_index.value(), this_epoch_progress);
            }
            if (residual_log) {
                fclose(residual_log);
            }

            book_writer->write(initial, epoch_index.value(), atoms);
            progress->epoch_finished(epoch_index.value());

            for (auto &list: atoms) {
                list.clear();
            }
            ++epochs_processed;
        }
    }

private:
    static double compute_total_energy(const Array2D<double> &data) {
        double sum2 = 0.0;
        for (int h = 0; h < data.height(); ++h) {
            for (index_t i = 0; i < data.length(); ++i) {
                const double value = data[h][i];
                sum2 += value * value;
            }
        }
        return sum2;
    }
};

#endif //EMPI_WORKER_LOOP_H

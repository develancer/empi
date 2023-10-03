/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <thread>
#include "ExtenderLoop.h"
#include "Worker.h"

//////////////////////////////////////////////////////////////////////////////

Worker::Worker(Array2D<real> data, unsigned cpu_threads, OptimizationMode mode) : data(std::move(data)), mode(mode) {
    atom_queue = std::make_shared<TaskQueue<BasicAtomPointer>>();
    spectrogram_queue = std::make_shared<TaskQueue<SpectrogramRequest>>();
    reset();
    if (mode == OPTIMIZATION_GLOBAL) {
        for (unsigned i = 1; i < cpu_threads; ++i) {
            threads.emplace_back(ExtenderLoop{atom_queue});
        }
    }
}

Worker::~Worker() {
    if (atom_queue) {
        atom_queue->terminate();
    }
    if (spectrogram_queue) {
        spectrogram_queue->terminate();
    }
    for (auto &thread : threads) {
        thread.join();
    }
}

void Worker::add_calculator(std::unique_ptr<SpectrogramCalculator> calculator, bool prefer_long_fft) {
    if (!primary_calculator_loop) {
        primary_calculator_loop = std::make_unique<SpectrogramLoop>(spectrogram_queue, std::move(calculator), prefer_long_fft);
    } else {
        threads.emplace_back(SpectrogramLoop(spectrogram_queue, std::move(calculator), prefer_long_fft));
    }
}

void Worker::add_dictionary(std::unique_ptr<Dictionary> dictionary) {
    dictionaries.push_back(std::move(dictionary));
}

ExtendedAtomPointer Worker::get_next_atom() {
    std::list<SpectrogramRequest> requests;
    for (auto &dictionary : dictionaries) {
        dictionary->fetch_requests(updated_index_range, requests);
    }
    if (!requests.empty()) {
        if (!primary_calculator_loop) {
            throw std::logic_error("need at least one SpectrogramCalculator instance");
        }
        requests.sort([](const SpectrogramRequest &a, const SpectrogramRequest &b) {
            return a.window_length < b.window_length;
        });

        spectrogram_queue->put(std::move(requests));
        (*primary_calculator_loop)(false);
        spectrogram_queue->wait_for_tasks();
    }

    BasicAtomPointer best_match;
    for (auto &dictionary : dictionaries) {
        BasicAtomPointer match = dictionary->get_best_match();
        if (match && (!best_match || *best_match < *match)) {
            best_match = match;
        }
    }

    if (!best_match) {
        // strangely, no more atoms
        return nullptr;
    }

    ExtendedAtomPointer best_atom = best_match->extend(mode >= OPTIMIZATION_LOCAL);
    assert(best_atom);
    if (mode >= OPTIMIZATION_GLOBAL) {
        // we also have to check other candidates
        std::list<BasicAtomPointer> candidates;
        for (auto &dictionary : dictionaries) {
            candidates.splice(candidates.end(), dictionary->get_candidate_matches(best_atom->energy));
        }
        atom_queue->put(std::list<BasicAtomPointer>(candidates));
        ExtenderLoop{atom_queue}(false);
        atom_queue->wait_for_tasks();
        for (const auto &candidate : candidates) {
            ExtendedAtomPointer another_atom = candidate->extend(true);
            if (another_atom->energy > best_atom->energy) {
                best_atom = another_atom;
            }
        }
    }

    updated_index_range = best_atom->subtract_from_signal();
    return best_atom;
}

std::list<ProtoRequest> Worker::get_proto_requests() {
    std::list<ProtoRequest> proto_requests;
    for (auto &dictionary : dictionaries) {
        dictionary->fetch_proto_requests(proto_requests);
    }
    return proto_requests;
}

void Worker::reset() {
    updated_index_range = data.length();
}

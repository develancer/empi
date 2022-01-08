/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <thread>
#include "BlockAtom.h"
#include "Computer.h"
#include "ExtenderLoop.h"
#include "WorkerLoop.h"

static void set_thread_affinity(std::thread &thread, unsigned cpu_index) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(cpu_index, &cpuset);
    pthread_setaffinity_np(thread.native_handle(), sizeof(cpu_set_t), &cpuset);
}

//////////////////////////////////////////////////////////////////////////////

Computer::Computer(Array2D<real> data, unsigned cpu_threads, OptimizationMode mode) : data(std::move(data)), mode(mode) {
    atom_queue = std::make_shared<TaskQueue<BasicAtomPointer>>();
    task_queue = std::make_shared<TaskQueue<SpectrogramRequest>>();
    reset();
    if (mode == OPTIMIZATION_GLOBAL) {
        for (unsigned i = 0; i < cpu_threads; ++i) {
            threads.emplace_back(ExtenderLoop(atom_queue));
            set_thread_affinity(threads.back(), i % std::thread::hardware_concurrency());
        }
    }
}

Computer::~Computer() {
    if (atom_queue) {
        atom_queue->terminate();
    }
    if (task_queue) {
        task_queue->terminate();
    }
    for (auto &thread : threads) {
        thread.join();
    }
}

void Computer::add_calculator(std::unique_ptr<Worker> calculator, bool prefer_long_fft) {
    static int cpu_index = 0;
    threads.emplace_back(WorkerLoop(task_queue, std::move(calculator), prefer_long_fft));
    set_thread_affinity(threads.back(), cpu_index++ % std::thread::hardware_concurrency());
}

void Computer::add_dictionary(std::unique_ptr<Dictionary> dictionary) {
    dictionaries.push_back(std::move(dictionary));
}

ExtendedAtomPointer Computer::get_next_atom() {
    std::list<SpectrogramRequest> requests;
    for (auto &dictionary : dictionaries) {
        dictionary->fetch_requests(updated_index_range, requests);
    }
    requests.sort([](const SpectrogramRequest &a, const SpectrogramRequest &b) {
        return a.window_length < b.window_length;
    });

    task_queue->put(std::move(requests));

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

std::list<ProtoRequest> Computer::get_proto_requests() {
    std::list<ProtoRequest> proto_requests;
    for (auto &dictionary : dictionaries) {
        dictionary->fetch_proto_requests(proto_requests);
    }
    return proto_requests;
}

void Computer::reset() {
    updated_index_range = data.length();
}

/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <algorithm>
#include <thread>
#include "Computer.h"

static void set_thread_affinity(std::thread &thread, unsigned cpu_index) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(cpu_index, &cpuset);
    pthread_setaffinity_np(thread.native_handle(), sizeof(cpu_set_t), &cpuset);
}

//////////////////////////////////////////////////////////////////////////////

void Computer::work(Worker *calculator) {
    SpectrogramRequest request;
    while (next_cpu(request)) {
        calculator->compute(request);
        request.interface->notify();
    }
}

void Computer::work_gpu(Worker *calculator) {
    SpectrogramRequest request;
    while (next_gpu(request)) {
        calculator->compute(request);
        request.interface->notify();
    }
}

bool Computer::next_cpu(SpectrogramRequest &request) {
    std::lock_guard lock(requests_mutex);
    bool result = !requests.empty();
    if (result) {
        request = requests.front();
        requests.pop_front();
    }
    return result;
}

bool Computer::next_gpu(SpectrogramRequest &request) {
    std::lock_guard lock(requests_mutex);
    bool result = !requests.empty();
    if (result) {
        request = requests.back();
        requests.pop_back();
    }
    return result;
}

void Computer::add_calculator(std::unique_ptr<Worker> calculator, bool prefer_long_fft) {
    if (prefer_long_fft) {
        calculators_gpu.push_back(std::move(calculator));
    } else {
        calculators_cpu.push_back(std::move(calculator));
    }
}

void Computer::add_dictionary(std::unique_ptr<Dictionary> dictionary) {
    dictionaries.push_back(std::move(dictionary));
}

ExtendedAtomPointer Computer::get_next_atom() {
    assert(requests.empty());
    for (auto &dictionary : dictionaries) {
        dictionary->fetch_requests(updated_index_range, requests);
    }
    requests.sort([](const SpectrogramRequest &a, const SpectrogramRequest &b) {
        return a.window_length < b.window_length;
    });

    std::vector<std::thread> threads;
    unsigned cpu_index = 0;
    for (auto &calculator : calculators_cpu) {
        threads.emplace_back(&Computer::work, this, calculator.get());
        set_thread_affinity(threads.back(), cpu_index++ % std::thread::hardware_concurrency());
    }
    for (auto &calculator : calculators_gpu) {
        threads.emplace_back(&Computer::work_gpu, this, calculator.get());
    }
    for (auto &thread : threads) {
        thread.join();
    }

    BasicAtomPointer best_match;
    Dictionary *best_dictionary = nullptr;
    for (auto &dictionary : dictionaries) {
        BasicAtomPointer match = dictionary->get_best_match();
        if (match && (!best_match || *best_match < *match)) {
            best_match = match;
            best_dictionary = dictionary.get();
        }
    }

    if (!best_match) {
        // strangely, no more atoms
        return nullptr;
    }

    ExtendedAtomPointer best_atom = best_match->extend();
    assert(best_atom);
    updated_index_range = best_dictionary->subtract_from_signal(*best_atom);
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

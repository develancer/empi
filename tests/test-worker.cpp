/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <atomic>
#include <cstdio>
#include "SpectrogramCalculator.h"
#include "Testing.h"
#include "Worker.h"

class SpectrogramCalculatorForTest : public SpectrogramCalculator {
    void compute(const SpectrogramRequest &request) final {
        request.maxima[0].energy = 10.0 * (request.how_many * 13 % 100);
    }
};

class BlockImplForTest : public BlockInterface {
    ExtractedMaximum *maximum;

public:
    explicit BlockImplForTest(ExtractedMaximum *maximum) : maximum(maximum) {}

    void notify() final {
        maximum->bin_index = 1;
    }
};

class ExtendedAtomForTest : public ExtendedAtom {
public:
    explicit ExtendedAtomForTest(Array2D<double> data, double energy) : ExtendedAtom(std::move(data), energy) {}

    void export_atom(std::list<ExportedAtom> *) final {
        // nothing here
    }

    IndexRange subtract_from_signal() const final {
        return IndexRange();
    }
};

class BasicAtomForTest : public BasicAtom {
public:
    explicit BasicAtomForTest(Array2D<double> data, double energy) : BasicAtom(std::move(data), energy) {}

    [[nodiscard]] ExtendedAtomPointer extend(bool allow_optimization) final {
        return std::make_shared<ExtendedAtomForTest>(data, energy);
    }

    [[nodiscard]] double get_energy_upper_bound() const final {
        return energy;
    }
};

const int R = 100; // number of requests

class DictionaryImplForTest : public Dictionary, public BlockInterface {
    Array2D<double> data_;
    ExtractedMaximum maxima[R];

public:
    static std::atomic<int> notify_called_count;

    DictionaryImplForTest(Array2D<double> data) : data_(std::move(data)) {}

    size_t get_atom_count() final {
        return 1; // does not matter, really
    }

    BasicAtomPointer get_best_match() final {
        double max_energy = 0;
        for (const auto &maximum : maxima) {
            max_energy = std::max(max_energy, maximum.energy);
        }
        return std::make_shared<BasicAtomForTest>(data_, max_energy);
    }

    std::list<BasicAtomPointer> get_candidate_matches(double energy_to_exceed) final {
        return std::list<BasicAtomPointer>();
    }

    void fetch_requests(IndexRange signal_range, std::list<SpectrogramRequest> &requests) final {
        ASSERT_EQUALS(0, signal_range.first_index);
        ASSERT_EQUALS(55, signal_range.end_index);
        for (int r = 0; r < R; ++r) {
            SpectrogramRequest request;
            request.how_many = r + 1;
            request.maxima = &maxima[r];
            request.interface = this;
            requests.push_back(request);

            maxima[r].energy = 0.0;
            maxima[r].bin_index = 0;
        }
    }

    void fetch_proto_requests(std::list<ProtoRequest> &) final {
        // nothing here
    }

    void notify() final {
        notify_called_count++;
    }
};

std::atomic<int> DictionaryImplForTest::notify_called_count;

int main() {
    const int W = 7; // number of workers
    Array2D<double> data(2, 55);

    Worker worker(data, 1);
    worker.add_dictionary(std::make_unique<DictionaryImplForTest>(data));
    for (int w = 0; w < W; ++w) {
        worker.add_calculator(std::make_unique<SpectrogramCalculatorForTest>());
    }
    ExtendedAtomPointer atom = worker.get_next_atom();
    ASSERT_EQUALS(10.0 * (R - 1), atom->energy);
    ASSERT_EQUALS(R, DictionaryImplForTest::notify_called_count);
    puts("OK");
}

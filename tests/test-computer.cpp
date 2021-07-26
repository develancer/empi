#include <atomic>
#include <cassert>
#include <cstdio>

#include "Computer.h"
#include "Worker.h"

class SpectrogramCalculatorForTest : public Worker {
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
};

class BasicAtomForTest : public BasicAtom {
public:
    explicit BasicAtomForTest(Array2D<double> data, double energy) : BasicAtom(std::move(data), energy) {}

    ExtendedAtomPointer extend() override {
        return std::make_shared<ExtendedAtomForTest>(data(), get_energy());
    }
};

const int R = 100; // number of requests

class DictionaryImplForTest : public Dictionary, public BlockInterface {
    Array2D<double> data_;
    ExtractedMaximum maxima[R];

public:
    static std::atomic<int> notify_called_count;

    DictionaryImplForTest(Array2D<double> data) : data_(std::move(data))  { }

    size_t get_atom_count() final {
        return 1; // does not matter, really
    }

    BasicAtomPointer get_best_match() final {
        double max_energy = 0;
        for (const auto& maximum : maxima) {
            max_energy = std::max(max_energy, maximum.energy);
        }
        return std::make_shared<BasicAtomForTest>(data_, max_energy);
    }

    void fetch_requests(IndexRange signal_range, std::list<SpectrogramRequest>& requests) final {
        assert(signal_range.first_index == 0);
        assert(signal_range.end_index == 55);
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

    void fetch_proto_requests(std::list<ProtoRequest>&) final {
        // nothing here
    }

    IndexRange subtract_from_signal(const ExtendedAtom&) final {
        return {0, 0};
    }

    void notify() final {
        notify_called_count++;
    }
};

std::atomic<int> DictionaryImplForTest::notify_called_count;

int main() {
    const int W = 7; // number of workers
    Array2D<double> data(2, 55);

    Computer computer(data);
    computer.add_dictionary(std::make_unique<DictionaryImplForTest>(data));
    for (int w = 0; w < W; ++w) {
        computer.add_calculator(std::make_unique<SpectrogramCalculatorForTest>());
    }
    ExtendedAtomPointer atom = computer.get_next_atom();
    assert(atom->get_energy() == 10.0 * (R-1));
    assert(DictionaryImplForTest::notify_called_count == R);
    puts("OK");
}

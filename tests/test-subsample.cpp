/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include "BlockDictionary.h"
#include "BlockDictionaryStructure.h"
#include "GaussianFamily.h"
#include "SpectrogramCalculatorFFTW.h"
#include "OptimizationMode.h"
#include "Testing.h"
#include "Worker.h"

static const double ENERGY_ERROR = 0.001;
static const double SCALE = 10.0;

void test_subsample(double position, double frequency, double phase = 0.0) {
    auto family = std::make_shared<GaussianFamily>();
    PinnedArray2D<double> data(1, 1024);
    data.fill(0.0);

    index_t offset, sample_count = family->size_for_values(position, SCALE, &offset);
    family->generate_values(position, SCALE, NULL, data[0]+offset, true);

    const double omega = 2 * M_PI * frequency;
    for (int i=0; i<sample_count; ++i) {
        data[0][offset+i] *= std::cos(omega * (offset + i - position) + phase);
    }

    BlockDictionaryStructure structure(family, ENERGY_ERROR, SCALE, SCALE, 0.5);
    auto calculator = std::make_unique<SpectrogramCalculatorFFTW>(1, std::set<int>{ 512 });
    auto dictionary = std::make_unique<BlockDictionary>(structure, data, extractorSingleChannel, *calculator, true);

    Worker worker(data, 1, OPTIMIZATION_DISABLED);
    worker.add_calculator(std::move(calculator));
    worker.add_dictionary(std::move(dictionary));
    auto result = worker.get_next_atom();

    std::list<ExportedAtom> atoms;
    result->export_atom(&atoms);
    ASSERT_EQUALS(1, atoms.size());

    const ExportedAtom& atom = atoms.front();
    ASSERT_EQUALS(position, atom.position);
    ASSERT_EQUALS(frequency, atom.frequency);
    ASSERT_SAME_PHASE(phase, atom.phase, 1.0e-12);
}

int main() {
    test_subsample(511.25, 0/512.0, 0.0);
    test_subsample(512.00, 50/512.0, 1.5);
    test_subsample(512.75, 100/512.0, 1.0);
    test_subsample(513.50, 150/512.0, 0.5);

    puts("OK");
}

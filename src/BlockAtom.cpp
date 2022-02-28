/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <array>
#include <cmath>
#include <memory>
#include <vector>
#include "nelder_mead.h"
#include "Array.h"
#include "Atom.h"
#include "BlockAtom.h"
#include "BlockAtomObjective.h"
#include "Configuration.h"
#include "Corrector.h"
#include "Logger.h"

//////////////////////////////////////////////////////////////////////////////

BlockAtom::BlockAtom(Array2D<double> data, double energy, double energy_upper_bound, std::shared_ptr<Family> family_,
                     double frequency, double position, double scale,
                     Extractor extractor, std::shared_ptr<BlockAtomParamsConverter> converter_)
        : BasicAtom(std::move(data), energy), energy_upper_bound(energy_upper_bound),
          extractor(extractor), family(std::move(family_)), converter(std::move(converter_)), params(frequency, position, scale) {}

void BlockAtom::connect_cache(std::shared_ptr<BlockAtomCache> cache, size_t key) {
    this->cache_slot = BlockAtomCacheSlot(std::move(cache), key);
}

ExtendedAtomPointer BlockAtom::extend(bool allow_optimization) {
    BlockAtomObjective objective(family, data, extractor, converter);

    std::array<double, 3> array = converter->arrayFromParams(params);
    if (allow_optimization) {
        if (cache_slot) {
            auto result = cache_slot->get();
            if (result) {
                return result;
            }
        }

        std::array<double, 3> step{0.5, 0.5, 0.5};
        auto result = nelder_mead<double, 3>(
                objective,
                array,
                Configuration::optimization_target,
                step,
                1,
                Configuration::optimization_max_iterations
        );
        if (result.ifault) {
            Logger::info("Parameter optimization could not converge. Unless it happens very often, it should not affect the decomposition. Tweaking --opt-max-iter and --opt-target options might help.");
        } else {
            array = result.xmin;
        }
    }

    double norm;
    const int channel_count = data.height();
    Array1D<ExtraData> extra_data(channel_count);
    double fit_energy = objective.calculate_energy(array, &norm, extra_data.get());
    BlockAtomParams fit = converter->paramsFromArray(array).first;

    const double amplitude_factor = norm * family->value(0.0);
    for (int i = 0; i < channel_count; ++i) {
        extra_data[i].amplitude *= amplitude_factor;
    }

    auto result = std::make_shared<BlockExtendedAtom>(
            data, fit_energy, family,
            fit.frequency, fit.position, fit.scale,
            std::move(extra_data)
    );
    if (allow_optimization && cache_slot) {
        IndexRange range = family->compute_range(result->params.position, result->params.scale);
        cache_slot->set(range, result);
    }

    return std::static_pointer_cast<ExtendedAtom>(result);
}

/**
 * @return estimate for maximum possible energy of the atom that can be obtained by locally optimizing the coefficients
 */
[[nodiscard]] double BlockAtom::get_energy_upper_bound() const {
    return energy_upper_bound;
}

//////////////////////////////////////////////////////////////////////////////

BlockExtendedAtom::BlockExtendedAtom(Array2D<double> data, double energy, std::shared_ptr<Family> family, double frequency, double position,
                                     double scale,
                                     Array1D<ExtraData> extra)
        : ExtendedAtom(std::move(data), energy), family(std::move(family)), params(frequency, position, scale), extra(std::move(extra)) {}

void BlockExtendedAtom::export_atom(std::list<ExportedAtom> *atoms) {
    for (index_t c = 0; c < extra.length(); ++c) {
        const ExtraData &e = extra[c];
        ExportedAtom atom(e.energy);
        atom.amplitude = e.amplitude;
        atom.envelope = family->name();
        atom.frequency = params.frequency;
        atom.scale = params.scale;
        atom.phase = e.phase;
        atom.position = params.position;
        atoms[c].push_back(std::move(atom));
    }
}

IndexRange BlockExtendedAtom::subtract_from_signal() const {
    index_t first_sample_offset;
    const index_t sample_count = family->size_for_values(params.position, params.scale, &first_sample_offset);
    index_t end_sample_offset = first_sample_offset + sample_count;
    if (end_sample_offset <= 0) {
        return {0, 0};
    }

    const index_t first_valid_sample_offset = std::max<index_t>(0, first_sample_offset);
    const index_t end_valid_sample_offset = std::min<index_t>(data.length(), end_sample_offset);

    Array1D<double> samples(sample_count);
    family->generate_values(params.position, params.scale, nullptr, samples.get(), false);
    double common_amplitude_factor = 1.0 / family->value(0.0);

    const int channel_count = data.height();
    const double omega = 2 * M_PI * params.frequency;
    for (int c = 0; c < channel_count; ++c) {
        const double channel_amplitude_factor = extra[c].amplitude * common_amplitude_factor;
        for (index_t i = first_valid_sample_offset; i < end_valid_sample_offset; ++i) {
            data[c][i] -= channel_amplitude_factor * samples[i - first_sample_offset] *
                           std::cos(omega * (static_cast<double>(i) - params.position) + extra[c].phase);
        }
    }

    return {first_valid_sample_offset, end_valid_sample_offset};
}

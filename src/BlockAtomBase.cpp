/**********************************************************
 * Piotr T. Różański (c) 2015-2021                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include "BlockAtomBase.h"

//////////////////////////////////////////////////////////////////////////////

BlockAtomParamsConverter::BlockAtomParamsConverter(double frequency_step, double position_step, double log_scale_step, double frequency_max)
        : frequency_step(frequency_step), position_step(position_step), log_scale_step(log_scale_step),
          max_scaled_frequency(frequency_max / frequency_step) {}

std::array<double, 3> BlockAtomParamsConverter::arrayFromParams(const BlockAtomParams &params) const {
    return {
            params.frequency / frequency_step,
            params.position / position_step,
            std::log(params.scale) / log_scale_step
    };
}

std::pair<BlockAtomParams, double> BlockAtomParamsConverter::paramsFromArray(const std::array<double, 3> &array) const {
    double scaled_frequency = array[0];
    double scaled_position = array[1];
    double scaled_log_scale = array[2];

    double penalty = fix_frequency(scaled_frequency) + fix_log_scale(scaled_log_scale);

    return std::make_pair(
            BlockAtomParams(
                    scaled_frequency * frequency_step,
                    scaled_position * position_step,
                    std::exp(scaled_log_scale * log_scale_step)
            ),
            penalty
    );
}

double BlockAtomParamsConverter::fix_frequency(double &scaled_frequency) const {
    return fix_scaled_argument(scaled_frequency, 0.0, max_scaled_frequency);
}

double BlockAtomParamsConverter::fix_log_scale(double &scaled_log_scale) const {
    return 0.0;
}

double BlockAtomParamsConverter::fix_scaled_argument(double &value, double min, double max) {
    double penalty = 0.0;
    if (value < min) {
        penalty += min - value;
        value = min;
    } else if (value > max) {
        penalty += value - max;
        value = max;
    }
    return penalty;
}

//////////////////////////////////////////////////////////////////////////////

BlockAtomParamsConverterBounded::BlockAtomParamsConverterBounded(double frequency_step, double position_step, double log_scale_step,
                                                                 double frequency_max, double min_scale, double max_scale)
        : BlockAtomParamsConverter(frequency_step, position_step, log_scale_step, frequency_max),
          min_scaled_log_scale(std::log(min_scale) / log_scale_step),
          max_scaled_log_scale(std::log(max_scale) / log_scale_step) {}

double BlockAtomParamsConverterBounded::fix_log_scale(double &scaled_log_scale) const {
    return fix_scaled_argument(scaled_log_scale, min_scaled_log_scale, max_scaled_log_scale);
}

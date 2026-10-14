/*!
 * \file secondary_code_synchronizer.h
 * \brief Statistical (multi-cycle) phase search for a known secondary/overlay code.
 * \author GNSS-SDR contributors, 2026
 *
 * -----------------------------------------------------------------------------
 *
 * GNSS-SDR is a Global Navigation Satellite System software-defined receiver.
 * This file is part of GNSS-SDR.
 *
 * Copyright (C) 2010-2026  (see AUTHORS file for a list of contributors)
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * -----------------------------------------------------------------------------
 */

#ifndef GNSS_SDR_SECONDARY_CODE_SYNCHRONIZER_H
#define GNSS_SDR_SECONDARY_CODE_SYNCHRONIZER_H

#include <complex>
#include <cstdint>
#include <string>
#include <vector>

/** \addtogroup Tracking
 * \{ */
/** \addtogroup Tracking_libs
 * \{ */


/*!
 * \brief Finds the phase of a known secondary/overlay code (e.g. GPS L5Q's
 * NH20) via differential demodulation, the same way a single-shot
 * exact-match correlator does, but accumulates a running correlation score
 * per candidate phase hypothesis across many code cycles instead of
 * requiring one perfect window.
 *
 * A wrong hypothesis's score behaves like a random walk around zero; the
 * true hypothesis's score grows roughly linearly with the number of
 * epochs processed (apart from an occasional symbol error). Lock is
 * declared once one hypothesis's score is both dominant over the
 * runner-up and stable across several consecutive evaluations -- the same
 * statistical-confidence shape as HistogramBitSynchronizer, applied to a
 * known pattern instead of generic transition events.
 *
 * A tracking-quality gap (bad prompt magnitude/quality) resets all
 * accumulated state: the running differential-bit chain is a single
 * global reference, so splicing scores accumulated before and after a gap
 * would silently mix two segments with an unknown relative polarity.
 */
class SecondaryCodeSynchronizer
{
public:
    struct Config
    {
        double dominance_ratio;
        double runner_up_margin;
        int stable_best_required;
        int min_score_for_lock;
        float min_prompt_mag;

        Config()
            : dominance_ratio(0.7),
              runner_up_margin(0.3),
              stable_best_required(3),
              min_score_for_lock(10),
              min_prompt_mag(0.0F)
        {
        }
    };

    SecondaryCodeSynchronizer(const std::string& secondary_code_string, const Config& cfg);

    void reset();

    /*!
     * \brief Update with one new coherent Prompt correlator output (one
     * secondary-code chip period, e.g. 1 ms for GPS L5Q).
     * \return true only on the epoch when lock is first declared.
     */
    bool update(const std::complex<float>& prompt, bool tracking_quality_ok);

    bool locked() const { return locked_; }

    /*!
     * \brief Secondary-code string index the *next* symbol (the one after
     * the one just passed to update()) should be compared/sign-flipped
     * against -- same contract dll_pll_veml_tracking::d_current_symbol
     * already has after a successful acquire_secondary().
     */
    int next_symbol_phase() const;

    /*!
     * \brief 180-degree carrier-phase ambiguity resolved as a side effect
     * of the correlation, same meaning/sign convention as
     * dll_pll_veml_tracking::d_Flag_PLL_180_deg_phase_locked.
     */
    bool phase_180_deg_locked() const { return phase_180_deg_locked_; }

private:
    void evaluate_lock();

    std::string secondary_code_;
    int N_;
    Config cfg_;

    std::vector<std::int32_t> score_;

    bool has_prev_prompt_;
    std::complex<float> prev_prompt_;
    bool running_bit_;

    std::int64_t epoch_count_;
    std::int64_t n_scored_;

    bool locked_;
    int locked_hypothesis_;
    bool phase_180_deg_locked_;

    bool has_last_best_hypothesis_;
    int last_best_hypothesis_;
    int stable_best_count_;
};


/** \} */
/** \} */
#endif  // GNSS_SDR_SECONDARY_CODE_SYNCHRONIZER_H

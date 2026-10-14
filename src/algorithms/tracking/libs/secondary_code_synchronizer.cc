/*!
 * \file secondary_code_synchronizer.cc
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

#include "secondary_code_synchronizer.h"
#include <algorithm>
#include <cmath>

SecondaryCodeSynchronizer::SecondaryCodeSynchronizer(const std::string& secondary_code_string, const Config& cfg)
    : secondary_code_(secondary_code_string),
      N_(static_cast<int>(secondary_code_string.size())),
      cfg_(cfg),
      has_prev_prompt_(false),
      running_bit_(false),
      epoch_count_(0),
      n_scored_(0),
      locked_(false),
      locked_hypothesis_(-1),
      phase_180_deg_locked_(false),
      has_last_best_hypothesis_(false),
      last_best_hypothesis_(0),
      stable_best_count_(0)
{
    score_.assign(N_ > 0 ? N_ : 0, 0);
}


void SecondaryCodeSynchronizer::reset()
{
    std::fill(score_.begin(), score_.end(), 0);
    has_prev_prompt_ = false;
    running_bit_ = false;
    epoch_count_ = 0;
    n_scored_ = 0;
    locked_ = false;
    locked_hypothesis_ = -1;
    phase_180_deg_locked_ = false;
    has_last_best_hypothesis_ = false;
    last_best_hypothesis_ = 0;
    stable_best_count_ = 0;
}


bool SecondaryCodeSynchronizer::update(const std::complex<float>& prompt, bool tracking_quality_ok)
{
    if (N_ <= 0 || locked_)
        {
            return false;
        }

    const std::int64_t t = epoch_count_;
    ++epoch_count_;

    if (!tracking_quality_ok || (std::abs(prompt) < cfg_.min_prompt_mag))
        {
            // A gap breaks the single global differential-bit chain that
            // every hypothesis's score is built on: splicing scores from
            // before and after an untrustworthy sample would silently mix
            // two segments with an unknown relative polarity. Restart
            // clean rather than risk corrupting the accumulated evidence.
            reset();
            return false;
        }

    if (!has_prev_prompt_)
        {
            prev_prompt_ = prompt;
            has_prev_prompt_ = true;
            // First symbol of a (re)started chain is the arbitrary polarity
            // reference; nothing to score against yet.
            return false;
        }

    const bool bit_flipped = (prompt * std::conj(prev_prompt_)).real() < 0.0F;
    running_bit_ = (running_bit_ != bit_flipped);
    prev_prompt_ = prompt;
    ++n_scored_;

    for (int h = 0; h < N_; h++)
        {
            const bool expected_bit = (secondary_code_[static_cast<std::size_t>((t + h) % N_)] == '1');
            score_[h] += (running_bit_ == expected_bit) ? 1 : -1;
        }

    evaluate_lock();
    return locked_;
}


void SecondaryCodeSynchronizer::evaluate_lock()
{
    if (n_scored_ < cfg_.min_score_for_lock)
        {
            return;
        }

    int best_h = 0;
    std::int32_t best_abs = std::abs(score_[0]);
    std::int32_t second_abs = 0;
    for (int h = 1; h < N_; h++)
        {
            const std::int32_t a = std::abs(score_[h]);
            if (a > best_abs)
                {
                    second_abs = best_abs;
                    best_abs = a;
                    best_h = h;
                }
            else if (a > second_abs)
                {
                    second_abs = a;
                }
        }

    if (best_abs < cfg_.min_score_for_lock)
        {
            return;
        }

    const auto total = static_cast<double>(n_scored_);
    const double dominance = static_cast<double>(best_abs) / total;
    const double margin = static_cast<double>(best_abs - second_abs) / total;
    if (dominance < cfg_.dominance_ratio || margin < cfg_.runner_up_margin)
        {
            has_last_best_hypothesis_ = false;
            stable_best_count_ = 0;
            return;
        }

    if (!has_last_best_hypothesis_ || best_h != last_best_hypothesis_)
        {
            last_best_hypothesis_ = best_h;
            has_last_best_hypothesis_ = true;
            stable_best_count_ = 1;
        }
    else
        {
            ++stable_best_count_;
        }

    if (stable_best_count_ >= cfg_.stable_best_required)
        {
            locked_ = true;
            locked_hypothesis_ = best_h;
            const bool last_symbol_negative = prev_prompt_.real() < 0.0F;
            const bool reconstructed_last_negative = !running_bit_;
            const bool invert_reconstruction = last_symbol_negative != reconstructed_last_negative;
            phase_180_deg_locked_ = (score_[best_h] < 0) != invert_reconstruction;
        }
}


int SecondaryCodeSynchronizer::next_symbol_phase() const
{
    if (!locked_ || N_ <= 0)
        {
            return 0;
        }
    return static_cast<int>((epoch_count_ + locked_hypothesis_) % N_);
}

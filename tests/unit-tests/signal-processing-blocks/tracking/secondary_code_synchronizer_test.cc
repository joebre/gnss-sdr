/*!
 * \file secondary_code_synchronizer_test.cc
 * \brief Deterministic tests for the statistical secondary-code phase synchronizer.
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
#include <gtest/gtest.h>
#include <complex>
#include <string>

namespace
{
// Real GPS L5Q NH20 pilot secondary code (IS-GPS-705), same string the
// production code uses (GPS_L5Q_NH_CODE_STR).
const std::string kGpsL5qNh20 = "00000100110101001110";

// Generates one symbol's sign for code-bit convention '0' -> +polarity,
// '1' -> -polarity, so a sign flip between consecutive symbols occurs
// exactly when the underlying code bit changes -- matching what a real
// differential (delta-phase) demodulator would see.
float symbol_sign(const std::string& code, int phase_offset, std::int64_t t, float base_polarity)
{
    const int n = static_cast<int>(code.size());
    const char bit = code[static_cast<std::size_t>(((t + phase_offset) % n + n) % n)];
    return (bit == '0') ? base_polarity : -base_polarity;
}

SecondaryCodeSynchronizer::Config fast_test_config()
{
    SecondaryCodeSynchronizer::Config cfg;
    cfg.min_score_for_lock = 8;
    cfg.dominance_ratio = 0.7;
    cfg.runner_up_margin = 0.3;
    cfg.stable_best_required = 2;
    cfg.min_prompt_mag = 0.0F;
    return cfg;
}
}  // namespace


TEST(SecondaryCodeSynchronizerTest, LocksAtKnownPhaseCleanSignal)
{
    SecondaryCodeSynchronizer sync(kGpsL5qNh20, fast_test_config());
    const int true_phase = 7;
    bool lock_event = false;
    for (std::int64_t t = 0; t < 5 * static_cast<std::int64_t>(kGpsL5qNh20.size()); ++t)
        {
            const float s = symbol_sign(kGpsL5qNh20, true_phase, t, 1.0F);
            lock_event = sync.update(std::complex<float>(s, 0.0F), true) || lock_event;
            if (sync.locked())
                {
                    break;
                }
        }
    EXPECT_TRUE(lock_event);
    ASSERT_TRUE(sync.locked());
}


TEST(SecondaryCodeSynchronizerTest, ReportsCorrectNextSymbolPhase)
{
    SecondaryCodeSynchronizer sync(kGpsL5qNh20, fast_test_config());
    const int true_phase = 3;
    std::int64_t last_t = -1;
    for (std::int64_t t = 0; t < 5 * static_cast<std::int64_t>(kGpsL5qNh20.size()); ++t)
        {
            const float s = symbol_sign(kGpsL5qNh20, true_phase, t, 1.0F);
            sync.update(std::complex<float>(s, 0.0F), true);
            last_t = t;
            if (sync.locked())
                {
                    break;
                }
        }
    ASSERT_TRUE(sync.locked());
    // The symbol right after the one just processed (index last_t+1) should
    // be compared against code position (last_t+1+true_phase) mod N.
    const int n = static_cast<int>(kGpsL5qNh20.size());
    const int expected = static_cast<int>(((last_t + 1 + true_phase) % n + n) % n);
    EXPECT_EQ(sync.next_symbol_phase(), expected);
}


TEST(SecondaryCodeSynchronizerTest, TwoPolaritiesGiveOppositePhaseLockFlag)
{
    SecondaryCodeSynchronizer sync_pos(kGpsL5qNh20, fast_test_config());
    SecondaryCodeSynchronizer sync_neg(kGpsL5qNh20, fast_test_config());
    const int true_phase = 5;
    for (std::int64_t t = 0; t < 5 * static_cast<std::int64_t>(kGpsL5qNh20.size()); ++t)
        {
            sync_pos.update(std::complex<float>(symbol_sign(kGpsL5qNh20, true_phase, t, 1.0F), 0.0F), true);
            sync_neg.update(std::complex<float>(symbol_sign(kGpsL5qNh20, true_phase, t, -1.0F), 0.0F), true);
        }
    ASSERT_TRUE(sync_pos.locked());
    ASSERT_TRUE(sync_neg.locked());
    EXPECT_NE(sync_pos.phase_180_deg_locked(), sync_neg.phase_180_deg_locked());
}


TEST(SecondaryCodeSynchronizerTest, TOLERATES_SCATTERED_SYMBOL_ERRORS)
{
    SecondaryCodeSynchronizer::Config cfg = fast_test_config();
    cfg.min_score_for_lock = 30;
    SecondaryCodeSynchronizer sync(kGpsL5qNh20, cfg);
    const int true_phase = 11;
    const int n = static_cast<int>(kGpsL5qNh20.size());

    // ~10% of symbols wrong, spread across many code cycles -- would never
    // achieve a single-window exact match, but should still lock given
    // enough accumulated evidence.
    for (std::int64_t t = 0; t < 15 * static_cast<std::int64_t>(n); ++t)
        {
            float s = symbol_sign(kGpsL5qNh20, true_phase, t, 1.0F);
            if ((t % 10) == 3)
                {
                    s = -s;  // inject an error
                }
            sync.update(std::complex<float>(s, 0.0F), true);
            if (sync.locked())
                {
                    break;
                }
        }
    ASSERT_TRUE(sync.locked());
}


TEST(SecondaryCodeSynchronizerTest, DoesNotFalseLockOnConstantInput)
{
    SecondaryCodeSynchronizer::Config cfg = fast_test_config();
    cfg.min_score_for_lock = 8;
    SecondaryCodeSynchronizer sync(kGpsL5qNh20, cfg);
    const int n = static_cast<int>(kGpsL5qNh20.size());

    // A constant-sign input never produces a differential transition, so
    // every phase hypothesis grows at the same rate (the code's own 0/1
    // imbalance) -- no hypothesis should ever become dominant.
    for (std::int64_t t = 0; t < 20 * static_cast<std::int64_t>(n); ++t)
        {
            sync.update(std::complex<float>(1.0F, 0.0F), true);
        }
    EXPECT_FALSE(sync.locked());
}


TEST(SecondaryCodeSynchronizerTest, TrackingQualityGapResetsAccumulatedState)
{
    SecondaryCodeSynchronizer sync(kGpsL5qNh20, fast_test_config());
    const int true_phase = 2;
    const int n = static_cast<int>(kGpsL5qNh20.size());

    // Feed less than min_score_for_lock worth of clean, correct evidence...
    for (std::int64_t t = 0; t < 3; ++t)
        {
            sync.update(std::complex<float>(symbol_sign(kGpsL5qNh20, true_phase, t, 1.0F), 0.0F), true);
        }
    EXPECT_FALSE(sync.locked());

    // ...then a quality gap. State must reset; immediately resuming the
    // *same* phase alignment that was interrupted should not let old and
    // new evidence silently combine into an early, wrong lock.
    sync.update(std::complex<float>(0.0F, 0.0F), false);

    for (std::int64_t t = 0; t < 3 * static_cast<std::int64_t>(n); ++t)
        {
            sync.update(std::complex<float>(symbol_sign(kGpsL5qNh20, true_phase, t, 1.0F), 0.0F), true);
            if (sync.locked())
                {
                    break;
                }
        }
    ASSERT_TRUE(sync.locked());
}


TEST(SecondaryCodeSynchronizerTest, ResetClearsStateAndAllowsReacquisition)
{
    SecondaryCodeSynchronizer sync(kGpsL5qNh20, fast_test_config());
    const int n = static_cast<int>(kGpsL5qNh20.size());

    for (std::int64_t t = 0; t < 5 * static_cast<std::int64_t>(n); ++t)
        {
            sync.update(std::complex<float>(symbol_sign(kGpsL5qNh20, 4, t, 1.0F), 0.0F), true);
            if (sync.locked())
                {
                    break;
                }
        }
    ASSERT_TRUE(sync.locked());

    sync.reset();
    EXPECT_FALSE(sync.locked());

    for (std::int64_t t = 0; t < 5 * static_cast<std::int64_t>(n); ++t)
        {
            sync.update(std::complex<float>(symbol_sign(kGpsL5qNh20, 9, t, 1.0F), 0.0F), true);
            if (sync.locked())
                {
                    break;
                }
        }
    EXPECT_TRUE(sync.locked());
}

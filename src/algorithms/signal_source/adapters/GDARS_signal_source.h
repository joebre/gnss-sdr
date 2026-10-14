/*!
 * \file GDARS_signal_source.h
 * \brief GDARS 4-channel record/playback signal source
 *
 * Reads a shared, interleaved GDARS capture file (see GDARS_source.h) with
 * up to 4 independently enabled raw channels (CH0=E5, CH1=L2, CH2=L1/E1,
 * CH3=E6), each fanning out to one or more RF-output ports (direct copy,
 * or CPU frequency-translated/decimated -- e.g. splitting a wideband
 * channel into E5a and E5b). One GDARS_source block instance per GLOBAL
 * RF-out port (each genuinely single-port, so no adapter-level tap is
 * needed -- see GDARS_source.h's top doc comment), exposed to the rest of
 * GNSS-SDR as a flat set of RF_channels via getRfChannels()/
 * get_right_block(int), same mechanism already used for any other
 * multi-RF-channel signal source (see gnss_flowgraph.cc's
 * connect_signal_sources_to_signal_conditioners()).
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

#ifndef GNSS_SDR_GDARS_SIGNAL_SOURCE_H
#define GNSS_SDR_GDARS_SIGNAL_SOURCE_H

#include "GDARS_source.h"
#include "concurrent_queue.h"
#include "signal_source_base.h"
#include <cstdint>
#include <string>
#include <vector>

class ConfigurationInterface;

/** \addtogroup Signal_Source
 * \{ */
/** \addtogroup Signal_Source_adapters
 * \{ */


/*!
 * \brief This class reads a GDARS 4-channel capture file.
 */
class GDARSSignalSource : public SignalSourceBase
{
public:
    GDARSSignalSource(
        const ConfigurationInterface* configuration,
        const std::string& role,
        unsigned int in_stream,
        unsigned int out_stream,
        Concurrent_Queue<pmt::pmt_t>* queue);

    ~GDARSSignalSource() = default;

    inline size_t item_size() override
    {
        return item_size_;
    }

    size_t getRfChannels() const override;

    void connect(gr::top_block_sptr top_block) override;
    void disconnect(gr::top_block_sptr top_block) override;
    gr::basic_block_sptr get_right_block() override;
    gr::basic_block_sptr get_right_block(int RF_channel) override;

private:
    // The one shared decode engine (single background thread) feeding
    // every port below -- see GDARS_source.h.
    GDARSDecodeEngine_sptr engine_;

    // One entry per GLOBAL rf_out port (flat 0..N-1 across all enabled
    // channels). Each is already a true single-port block, directly
    // addressable via get_right_block(int) -- no adapter-level tap needed.
    std::vector<GDARS_source_sptr> ports_;

    size_t item_size_;
};


/** \} */
/** \} */
#endif  // GNSS_SDR_GDARS_SIGNAL_SOURCE_H

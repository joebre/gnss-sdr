/*!
 * \file GDARS_signal_source.cc
 * \brief GDARS 4-channel record/playback signal source; see GDARS_signal_source.h.
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

#include "GDARS_signal_source.h"
#include "configuration_interface.h"
#include <gnuradio/blocks/file_sink.h>
#include <gnuradio/filter/firdes.h>
#include <array>
#include <cmath>
#include <utility>

#if USE_GLOG_AND_GFLAGS
#include <glog/logging.h>
#else
#include <absl/log/log.h>
#endif

using namespace std::string_literals;

namespace
{
constexpr unsigned int kNumRawChannels = 4;

gr::fft::window::win_type parse_window_type(const std::string& role, const std::string& name)
{
    if (name == "hamming")
        {
            return gr::fft::window::win_type::WIN_HAMMING;
        }
    if (name == "hann")
        {
            return gr::fft::window::win_type::WIN_HANN;
        }
    if (name == "blackman")
        {
            return gr::fft::window::win_type::WIN_BLACKMAN;
        }
    if (name == "blackman_harris")
        {
            return gr::fft::window::win_type::WIN_BLACKMAN_HARRIS;
        }
    if (name == "rectangular")
        {
            return gr::fft::window::win_type::WIN_RECTANGULAR;
        }
    if (name == "kaiser")
        {
            return gr::fft::window::win_type::WIN_KAISER;
        }
    LOG(WARNING) << role << ": unknown window type '" << name << "', falling back to hamming";
    return gr::fft::window::win_type::WIN_HAMMING;
}

// Same as evk1029_signal_source.cc's helper -- Kaiser's own empirical
// formula (1974) relating a target stopband attenuation to the window
// beta that actually achieves it (firdes::low_pass_2() sizes the tap count
// from attenuation_dB but doesn't derive beta from it itself).
double kaiser_beta_from_attenuation(double attenuation_dB)
{
    if (attenuation_dB > 50.0)
        {
            return 0.1102 * (attenuation_dB - 8.7);
        }
    if (attenuation_dB >= 21.0)
        {
            return 0.5842 * std::pow(attenuation_dB - 21.0, 0.4) + 0.07886 * (attenuation_dB - 21.0);
        }
    return 0.0;
}

struct ChannelConfig
{
    unsigned int rf_out = 0;
    int64_t fs = 0;
};
}  // namespace


GDARSSignalSource::GDARSSignalSource(
    const ConfigurationInterface* configuration,
    const std::string& role,
    unsigned int in_stream,
    unsigned int out_stream,
    Concurrent_Queue<pmt::pmt_t>* queue)
    : SignalSourceBase(configuration, role, "GDARS_Signal_Source"s),
      item_size_(sizeof(gr_complex))
{
    const std::string filename = configuration->property(role + ".filename"s, "data.bin"s);
    const int64_t default_fs = configuration->property("GNSS-SDR.internal_fs_sps"s, int64_t(0));

    // First pass: validate each channel (same checks as the second pass
    // below) and build a FLAT list of global ports, one entry per port in
    // order, each tagged with its own raw channel id. The one shared
    // decode engine -- opened exactly once regardless of how many ports
    // read from the same channel -- needs this mapping up front to know
    // which channel ids to decode and which global port(s) to fan each
    // one out to.
    std::array<ChannelConfig, kNumRawChannels> channel_config{};
    std::array<bool, kNumRawChannels> channel_enabled{};
    std::vector<uint32_t> port_channel_ids;
    std::vector<uint32_t> channel_ids;
    std::vector<double> channel_sampling_frequencies;

    for (unsigned int ch = 0; ch < kNumRawChannels; ch++)
        {
            const std::string ch_prefix = role + ".ch"s + std::to_string(ch) + "_"s;
            if (!configuration->property(ch_prefix + "en"s, false))
                {
                    continue;
                }
            const unsigned int rf_out = configuration->property(ch_prefix + "rf_out"s, 0U);
            if (rf_out == 0)
                {
                    LOG(ERROR) << role << ": " << ch_prefix << "en is true but " << ch_prefix << "rf_out is 0 (or unset) -- skipping channel " << ch;
                    continue;
                }
            const int64_t fs = configuration->property(ch_prefix + "sampling_frequency"s, default_fs);
            if (fs <= 0)
                {
                    LOG(ERROR) << role << ": " << ch_prefix << "sampling_frequency is not set (and GNSS-SDR.internal_fs_sps is unavailable) -- skipping channel " << ch;
                    continue;
                }
            channel_enabled[ch] = true;
            channel_config[ch] = ChannelConfig{rf_out, fs};
            for (unsigned int p = 0; p < rf_out; p++)
                {
                    port_channel_ids.push_back(ch);
                }
            channel_ids.push_back(ch);
            channel_sampling_frequencies.push_back(static_cast<double>(fs));
        }

    engine_ = std::make_shared<GDARSDecodeEngine>(filename, port_channel_ids, channel_ids, channel_sampling_frequencies);

    unsigned int global_rf_out = 0;
    for (unsigned int ch = 0; ch < kNumRawChannels; ch++)
        {
            if (!channel_enabled[ch])
                {
                    continue;
                }
            const unsigned int rf_out = channel_config[ch].rf_out;
            const int64_t fs = channel_config[ch].fs;

            for (unsigned int p = 0; p < rf_out; p++)
                {
                    const std::string port_prefix = role + ".rf_out"s + std::to_string(global_rf_out) + "_"s;
                    GDARSFreqXlatingBand band;
                    band.enable_freq_xlating = configuration->property(port_prefix + "enable_freq_xlating"s, false);
                    if (band.enable_freq_xlating)
                        {
                            const double intermediate_freq = configuration->property(port_prefix + "IF"s, 0.0);
                            const int decimation_factor = configuration->property(port_prefix + "decimation_factor"s, 1);
                            const double default_bw = (static_cast<double>(fs) / decimation_factor) / 2;
                            const double bw = configuration->property(port_prefix + "bw"s, default_bw);
                            const double default_tw = bw / 10.0;
                            const double tw = configuration->property(port_prefix + "tw"s, default_tw);
                            const std::string window_name = configuration->property(port_prefix + "window"s, "hamming"s);
                            const gr::fft::window::win_type window = parse_window_type(role, window_name);
                            band.decimation_factor = decimation_factor;
                            band.center_freq = intermediate_freq;
                            band.use_cuda = configuration->property(port_prefix + "cuda"s, false);
                            if (configuration->is_present(port_prefix + "attenuation_dB"s))
                                {
                                    const double attenuation_dB = configuration->property(port_prefix + "attenuation_dB"s, 53.0);
                                    const double default_kaiser_beta = kaiser_beta_from_attenuation(attenuation_dB);
                                    const double kaiser_beta = configuration->property(port_prefix + "kaiser_beta"s, default_kaiser_beta);
                                    band.taps = gr::filter::firdes::low_pass_2(1.0, static_cast<double>(fs), bw, tw, attenuation_dB, window, kaiser_beta);
                                }
                            else
                                {
                                    const double kaiser_beta = configuration->property(port_prefix + "kaiser_beta"s, 6.76);
                                    band.taps = gr::filter::firdes::low_pass(1.0, static_cast<double>(fs), bw, tw, window, kaiser_beta);
                                }
                            LOG(INFO) << role << ": " << port_prefix << "freq-xlating, IF=" << intermediate_freq
                                      << " decimation=" << decimation_factor << " window=" << window_name << " taps=" << band.taps.size();
                        }
                    else
                        {
                            LOG(INFO) << role << ": " << port_prefix << "direct (no frequency translation)";
                        }

                    ports_.push_back(GDARS_make_source(engine_, global_rf_out, static_cast<double>(fs), queue, band));
                    global_rf_out++;
                }
        }

    if (ports_.empty())
        {
            LOG(ERROR) << role << ": no raw channel enabled (ch0_en/ch1_en/ch2_en/ch3_en all false, or all skipped due to a configuration error above)";
        }

    if (in_stream > 0)
        {
            LOG(ERROR) << "A signal source does not have an input stream";
        }
    if (out_stream > 1)
        {
            LOG(ERROR) << "This implementation only supports one output stream";
        }
}


size_t GDARSSignalSource::getRfChannels() const
{
    return ports_.size();
}


void GDARSSignalSource::connect(gr::top_block_sptr top_block __attribute__((unused)))
{
    // Nothing to wire here: every port is already its own true single-port
    // GDARS_source block (see GDARS_source.h), connected directly by
    // gnss_flowgraph.cc via get_right_block(int). No adapter-internal tap
    // is needed any more.
}


void GDARSSignalSource::disconnect(gr::top_block_sptr top_block __attribute__((unused)))
{
}


gr::basic_block_sptr GDARSSignalSource::get_right_block()
{
    return get_right_block(0);
}


gr::basic_block_sptr GDARSSignalSource::get_right_block(int RF_channel)
{
    if (RF_channel < 0 || static_cast<std::size_t>(RF_channel) >= ports_.size())
        {
            if (RF_channel == 0)
                {
                    LOG(ERROR) << "GDARSSignalSource::get_right_block: no RF output port is available";
                }
            return gr::block_sptr();
        }
    return ports_.at(static_cast<std::size_t>(RF_channel));
}

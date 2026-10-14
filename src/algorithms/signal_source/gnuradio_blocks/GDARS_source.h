/*!
 * \file GDARS_source.h
 * \brief Reads a GDARS 4-channel record/playback capture file (samples
 * interleaved across channels in 2048-byte blocks, each tagged with a
 * 3-bit channel id and a 3-bit resolution code in an 8-byte header;
 * resolution-dependent LUT decode into native quadrature I/Q) and fans
 * each enabled raw channel out to one or more frequency-translated/
 * decimated RF-output ports (e.g. a wideband channel carrying both E5a
 * and E5b).
 *
 * Threading model: ONE persistent background thread (owned by
 * GDARSDecodeEngine) does all file I/O and decode for every enabled raw
 * channel, fanning each channel's decoded samples out into the
 * GDARSRawPortBuffer of every RF-out port that reads it. ONE additional
 * persistent thread per RF-out port (owned by that port's own GDARS_source
 * instance) drains its raw buffer and either copies it straight through or
 * performs frequency translation (FIR + mixer), pushing the result into
 * that port's own output ring buffer; GNU Radio's work() call for that
 * port then only has to drain an already-filled buffer, so it is cheap
 * regardless of how small noutput_items is and never needs a raised
 * output_multiple.
 *
 * This replaces an earlier design (one gr::sync_decimator instance per
 * *raw channel*, each with its own ping-pong decode-ahead thread and
 * per-port BandWorker threads, plus a separate gr::blocks::multiply_const_cc
 * tap block per RF-out port) that needed the tap for two reasons: letting
 * gnss_flowgraph.cc's get_right_block(int) address each port as a true
 * single-port block despite several ports living on one multi-port
 * instance, and isolating that instance's raised output_multiple from the
 * tens of acq+trk readers sharing a single port (removing the tap once,
 * without also removing the raised output_multiple, reproduced the exact
 * "many readers + raised output_multiple" GNU Radio scheduler deadlock
 * Evk1029Source's own code documents avoiding -- confirmed live: a real
 * run hung permanently mid-startup, never reaching acquisition). The
 * current design sidesteps both reasons at the root: each GDARS_source
 * instance is already genuinely single-port (so get_right_block(int) can
 * return it directly, no tap needed), and its work() never needs a raised
 * output_multiple (so there is nothing to isolate).
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

#ifndef GNSS_SDR_GDARS_SOURCE_H
#define GNSS_SDR_GDARS_SOURCE_H

#include "concurrent_queue.h"
#include "gnss_block_interface.h"
#include <gnuradio/sync_block.h>
#include <pmt/pmt.h>
#include <volk/volk_alloc.hh>
#include <array>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#if CUDA_GPU_ACCEL
#include "cuda_ddc_engine.h"
#endif

/** \addtogroup Signal_Source
 * \{ */
/** \addtogroup Signal_Source_gnuradio_blocks
 * \{ */


class GDARS_source;
using GDARS_source_sptr = gnss_shared_ptr<GDARS_source>;

//! Per-port (rf_out) configuration, parsed by the adapter from config.
struct GDARSFreqXlatingBand
{
    bool enable_freq_xlating = false;  // false => direct: no FIR/mixer, just a copy
    int decimation_factor = 1;
    double center_freq = 0.0;  // IF, relative to this port's own raw channel baseband; unused when !enable_freq_xlating
    std::vector<float> taps;   // real low-pass prototype; empty when !enable_freq_xlating
    bool use_cuda = false;     // only takes effect when built with -DENABLE_CUDA=ON; see GDARS_source.cc
};

// Bounded, blocking single-producer/single-consumer queue of raw decoded
// gr_complex samples for ONE global RF-out port. Filled by
// GDARSDecodeEngine's single decode thread (push, blocking until there is
// room -- bounding this buffer's size is what bounds how far the decode
// thread can run ahead of this port's own processor thread); drained by
// that port's own GDARS_source processor thread (pop, blocking until data
// is available or the engine has signaled no more will ever arrive).
class GDARSRawPortBuffer
{
public:
    explicit GDARSRawPortBuffer(std::size_t capacity);

    void push(const gr_complex* data, std::size_t count);
    // Pops up to max_count samples into out; returns the count actually
    // popped. Blocks until at least min_count samples are available (so a
    // caller can batch many small pushes into one bigger pop, e.g.
    // GDARS_source::processor_loop() -- see its own quantum comment),
    // UNLESS input is done (or this buffer was abandoned) first, in which
    // case it returns whatever is available immediately (even less than
    // min_count, down to 0 with eof_out set) rather than blocking forever
    // waiting for a batch that will never complete.
    std::size_t pop(gr_complex* out, std::size_t max_count, std::size_t min_count, bool& eof_out);
    void mark_input_done();
    // Forcibly unblocks any in-progress or future push()/pop() call
    // without waiting for real EOF -- used by GDARS_source::stop() (the
    // GNU Radio shutdown hook) and GDARSDecodeEngine's destructor so a
    // blocked producer or consumer thread can exit promptly instead of
    // deadlocking flowgraph teardown.
    void abandon();

private:
    std::mutex mutex_;
    std::condition_variable cv_;
    std::deque<gr_complex> buffer_;
    std::size_t capacity_;
    bool input_done_ = false;
    bool abandoned_ = false;
};

// Owns the one shared, interleaved GDARS capture file and the single
// background thread that reads and decodes it for every enabled raw
// channel, fanning each channel's decoded samples out into the
// GDARSRawPortBuffer of every GLOBAL rf_out port that reads it (several
// ports may share one raw channel -- e.g. ch0_rf_out=2 -- each gets its
// own copy of that channel's decoded stream, since each port processes/
// decimates it independently downstream). A block whose header id is
// disabled or has no port mapped to it is discarded after reading just
// its header, same as before.
class GDARSDecodeEngine
{
public:
    // port_channel_ids[g] = the raw channel id (0-3) global port g reads
    // from. channel_ids/channel_sampling_frequencies: parallel arrays, one
    // entry per enabled raw channel, giving that channel's own native
    // sample rate -- used only to size this channel's own decode-
    // accumulation batch (see run()'s doc comment), same quantum
    // GDARS_source::processor_loop() uses downstream of it.
    GDARSDecodeEngine(const std::string& filename, const std::vector<uint32_t>& port_channel_ids,
        const std::vector<uint32_t>& channel_ids, const std::vector<double>& channel_sampling_frequencies);
    ~GDARSDecodeEngine();

    GDARSRawPortBuffer& port_buffer(std::size_t global_port);

private:
    void run();

    std::array<std::size_t, 8> channel_quantum_{};  // per-channel decode-accumulation threshold, see run()
    std::array<std::vector<std::size_t>, 8> ports_by_channel_;  // channel id -> global port indices reading it
    std::vector<std::unique_ptr<GDARSRawPortBuffer>> port_buffers_;

    std::ifstream file_;
    std::atomic<bool> stop_requested_{false};
    std::thread thread_;
};
using GDARSDecodeEngine_sptr = std::shared_ptr<GDARSDecodeEngine>;


GDARS_source_sptr GDARS_make_source(GDARSDecodeEngine_sptr engine, std::size_t global_port, double sampling_frequency,
    Concurrent_Queue<pmt::pmt_t>* queue, const GDARSFreqXlatingBand& band);

// One GNU Radio block PER GLOBAL RF-OUT PORT -- see this file's top doc
// comment for the threading model and why it replaced the earlier
// one-instance-per-raw-channel design.
class GDARS_source : public gr::sync_block
{
public:
    ~GDARS_source() override;

private:
    friend GDARS_source_sptr GDARS_make_source(GDARSDecodeEngine_sptr engine, std::size_t global_port, double sampling_frequency,
        Concurrent_Queue<pmt::pmt_t>* queue, const GDARSFreqXlatingBand& band);
    GDARS_source(GDARSDecodeEngine_sptr engine, std::size_t global_port, double sampling_frequency,
        Concurrent_Queue<pmt::pmt_t>* queue, const GDARSFreqXlatingBand& band);

    int work(int noutput_items,
        gr_vector_const_void_star& input_items,
        gr_vector_void_star& output_items) override;

    // GNU Radio's shutdown hook, called BEFORE the scheduler tries to join
    // this block's own work()-calling thread. work() can block indefinitely
    // inside out_buffer_.pop() waiting for processor_loop() to produce more
    // -- without this override unblocking that wait, flowgraph teardown
    // would deadlock (the TPB thread can't be joined until work() returns,
    // but work() would otherwise only return on genuine EOF).
    bool stop() override;

    // Persistent thread: drains raw_buffer_, does the copy/FIR+mixer work,
    // pushes the result into out_buffer_. Runs continuously, decoupled
    // from GNU Radio's own work() calling cadence -- so work() itself is
    // always just a cheap drain of already-ready output.
    void processor_loop();

    GDARSDecodeEngine_sptr engine_;  // keeps the shared decode thread alive as long as any port using it exists
    GDARSRawPortBuffer& raw_buffer_;
    std::size_t global_port_;
    Concurrent_Queue<pmt::pmt_t>* queue_;

    // Minimum raw (pre-decimation) samples processor_loop() accumulates
    // from raw_buffer_ before doing any copy/FIR work -- same ~0.7 ms
    // target and nearest-power-of-two rounding as evk1029_source.cc's own
    // output_multiple_from_sampling_frequency(), computed from this port's
    // raw channel sampling_frequency. Without this, processor_loop() was
    // woken and did real work on every single decoded raw channel-block
    // (at most ~8160 samples, often far fewer), nowhere near EVK1029's
    // quantum.
    std::size_t batch_quantum_;

    // Direct mode (no FIR/mixer): a pure copy at the raw channel's own
    // native rate (decimation_ == 1). Freq-xlating mode: composite FIR
    // (frequency-shifted taps, reversed to match fir_filter's convention)
    // + post-FIR mixer rotation, same math as the earlier Band struct.
    bool direct_ = true;
    volk::vector<gr_complex> composite_fir_;
    std::size_t history_len_ = 0;  // composite_fir_.size() - 1, 0 when direct_
    gr_complex phase_{1.0F, 0.0F};
    gr_complex phase_incr_{1.0F, 0.0F};
    int samples_since_phase_renorm_ = 0;
    int decimation_ = 1;
#if CUDA_GPU_ACCEL
    std::unique_ptr<CudaDdcEngine> cuda_engine_;  // null => CPU path (construction, set_taps(), or a prior compute_complex() failed)
#endif

    // Trailing history_len_ raw samples carried from the previous
    // processor_loop() iteration, for FIR continuity across pop() calls
    // (overlap-save), followed by freshly-popped new samples.
    std::vector<gr_complex> raw_window_;

    // Final, ready-to-serve output, filled by processor_loop(), drained by
    // work() -- reusing GDARSRawPortBuffer here too (not just for the raw
    // engine->port feed) gives it the same bounded-blocking-push
    // backpressure and the same abandon() escape hatch stop() needs.
    GDARSRawPortBuffer out_buffer_;
    std::atomic<bool> stop_{false};
    std::thread processor_thread_;
};


/** \} */
/** \} */
#endif  // GNSS_SDR_GDARS_SOURCE_H

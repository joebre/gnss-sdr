/*!
 * \file GDARS_source.cc
 * \brief Unpacks one RF-out port of a GDARS capture file; see GDARS_source.h.
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

#include "GDARS_source.h"
#include "command_event.h"
#include <gnuradio/io_signature.h>
#include <gnuradio/math.h>
#include <volk/volk.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>

#if defined(__linux__) || defined(__APPLE__)
#include <pthread.h>
#endif

namespace
{
// One raw "channel block" as recorded by the GDARS hardware: an 8-byte
// header (channel id, resolution, attenuation, ...) followed by
// SAMPLE_BLK_PER_CH_BLK 8-byte sub-blocks, each packing two DWORDs (Q then
// I) of variable-bit-depth samples. Decode logic ported verbatim from the
// pre-"next" GDARS_source -- it was already correct.
constexpr std::size_t BYTES_PER_CH_BLK = 2048;
constexpr std::size_t SAMPLE_BLK_PER_CH_BLK = 255;

// Read a whole number of channel-blocks at a time, amortizing I/O overhead
// (one syscall per ~512 blocks instead of one per block, including the
// blocks this decode engine's own filtering will just discard after
// reading their header). Must stay an exact multiple of BYTES_PER_CH_BLK
// so a block can never straddle a refill.
constexpr std::size_t kReadChunkBytes = 512 * BYTES_PER_CH_BLK;  // 1 MiB
static_assert(kReadChunkBytes % BYTES_PER_CH_BLK == 0, "must be a whole number of channel-blocks");

// The true maximum samples one channel-block can decode to (resolution 0,
// i.e. 1-bit samples: num_of_samples[0] == 32 per block). NOTE: the old
// GDARS_source's own MAX_SAMPLE_PER_BLK (16 * 255) was inconsistent with
// its own num_of_samples table (max 32, not 16) -- recomputed correctly
// here since it's actually relied on below to size buffers.
constexpr std::size_t kMaxSamplesPerBlock = 32U * SAMPLE_BLK_PER_CH_BLK;

// Capacity (in gr_complex samples) of every GDARSRawPortBuffer, both the
// engine->port raw feed and each port's own processed output buffer.
// Bounds how far the decode thread (or a port's processor thread) can run
// ahead of a slow consumer, and -- since read from a file is effectively
// unthrottled -- is what actually provides this adapter's pacing/
// backpressure. Matches the historical set_min_output_buffer() sizing
// used throughout this codebase (e.g. Evk1029Source).
constexpr std::size_t kPortBufferCapacity = 4'000'000;

// Same tuning/rationale as evk1029_source.cc's own
// output_multiple_from_sampling_frequency(): target ~0.7 ms of raw
// samples per processor_loop() batch, rounded to the nearest power of
// two, so GDARS_source::processor_loop() isn't woken and redoing FIR
// setup on every single tiny raw channel-block the decode engine happens
// to push. Unlike EVK1029's use of this quantum (gr::sync_block's own
// output_multiple, which raises a real GNU Radio scheduler deadlock risk
// when many readers share one port -- see this file's and the adapter's
// header doc comments), this is purely an internal threshold on our own
// background thread's pop() from raw_buffer_: it has no GNU Radio
// scheduler involvement, so none of that risk applies here.
constexpr double kBatchTargetSeconds = 0.0007;
constexpr std::size_t kBatchFallback = 131072;

std::size_t batch_quantum_from_sampling_frequency(double sampling_frequency)
{
    if (!(sampling_frequency > 0.0))
        {
            return kBatchFallback;
        }
    const double target_samples = sampling_frequency * kBatchTargetSeconds;
    const int lower_exp = static_cast<int>(std::floor(std::log2(std::max(target_samples, 1.0))));
    const double lower = std::exp2(static_cast<double>(lower_exp));
    const double upper = std::exp2(static_cast<double>(lower_exp + 1));
    return static_cast<std::size_t>((target_samples - lower <= upper - target_samples) ? lower : upper);
}

const float dec_LUT1[2] = {1, -1};
const float dec_LUT2[4] = {1, 3, -1, -3};
const float dec_LUT3[8] = {1, 3, 5, 7, -1, -3, -5, -7};
const float dec_LUT4[16] = {1, 3, 5, 7, 9, 11, 13, 15, -1, -3, -5, -7, -9, -11, -13, -15};
const float dec_LUT5[32] = {
    1, 3, 5, 7, 9, 11, 13, 15,
    17, 19, 21, 23, 25, 27,
    29, 31, -1, -3, -5, -7,
    -9, -11, -13, -15, -17,
    -19, -21, -23, -25, -27,
    -29, -31};
const float dec_LUT6[64] = {
    1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25,
    27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47,
    49, 51, 53, 55, 57, 59, 61, 63, -1, -3, -5,
    -7, -9, -11, -13, -15, -17, -19, -21, -23,
    -25, -27, -29, -31, -33, -35, -37, -39,
    -41, -43, -45, -47, -49, -51, -53, -55,
    -57, -59, -61, -63};
const float dec_LUT7[128] = {
    1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25,
    27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47,
    49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69,
    71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91,
    93, 95, 97, 99, 101, 103, 105, 107, 109,
    111, 113, 115, 117, 119, 121, 123, 125,
    127, -1, -3, -5, -7, -9, -11, -13, -15, -17,
    -19, -21, -23, -25, -27, -29, -31, -33,
    -35, -37, -39, -41, -43, -45, -47, -49,
    -51, -53, -55, -57, -59, -61, -63, -65,
    -67, -69, -71, -73, -75, -77, -79, -81,
    -83, -85, -87, -89, -91, -93, -95, -97,
    -99, -101, -103, -105, -107, -109, -111,
    -113, -115, -117, -119, -121, -123, -125, -127};
const float dec_LUT8[256] = {
    1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27,
    29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51,
    53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75,
    77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99,
    101, 103, 105, 107, 109, 111, 113, 115, 117,
    119, 121, 123, 125, 127, 129, 131, 133, 135,
    137, 139, 141, 143, 145, 147, 149, 151, 153,
    155, 157, 159, 161, 163, 165, 167, 169, 171,
    173, 175, 177, 179, 181, 183, 185, 187, 189,
    191, 193, 195, 197, 199, 201, 203, 205, 207,
    209, 211, 213, 215, 217, 219, 221, 223, 225,
    227, 229, 231, 233, 235, 237, 239, 241, 243,
    245, 247, 249, 251, 253, 255, -1, -3, -5, -7,
    -9, -11, -13, -15, -17, -19, -21, -23, -25,
    -27, -29, -31, -33, -35, -37, -39, -41, -43,
    -45, -47, -49, -51, -53, -55, -57, -59, -61,
    -63, -65, -67, -69, -71, -73, -75, -77, -79,
    -81, -83, -85, -87, -89, -91, -93, -95, -97,
    -99, -101, -103, -105, -107, -109, -111, -113,
    -115, -117, -119, -121, -123, -125, -127, -129,
    -131, -133, -135, -137, -139, -141, -143, -145,
    -147, -149, -151, -153, -155, -157, -159, -161,
    -163, -165, -167, -169, -171, -173, -175, -177,
    -179, -181, -183, -185, -187, -189, -191, -193,
    -195, -197, -199, -201, -203, -205, -207, -209,
    -211, -213, -215, -217, -219, -221, -223, -225,
    -227, -229, -231, -233, -235, -237, -239, -241,
    -243, -245, -247, -249, -251, -253, -255};

const float* const dec_LUTMatrix[8] = {
    dec_LUT1, dec_LUT2, dec_LUT3, dec_LUT4,
    dec_LUT5, dec_LUT6, dec_LUT7, dec_LUT8};

const unsigned int num_of_samples[8] = {32, 16, 10, 8, 6, 5, 4, 4};

struct block_header
{
    uint32_t id : 3;
    uint32_t fk : 1;
    uint32_t res : 3;
    uint32_t dec : 6;
    uint32_t out_att : 11;
    uint32_t reserved : 7;
    uint32_t pps_data : 1;
};

void name_current_thread(const std::string& name)
{
#if defined(__linux__)
    pthread_setname_np(pthread_self(), name.c_str());
#elif defined(__APPLE__)
    pthread_setname_np(name.c_str());
#endif
}

// Decodes one already-read 2048-byte block (whose header has already been
// checked to belong to a channel at least one port is reading) into out.
// Returns the sample count produced (num_of_samples[res] * SAMPLE_BLK_PER_CH_BLK).
std::size_t decode_channel_block(gr_complex* out, const uint8_t* block)
{
    const auto* hdr = reinterpret_cast<const block_header*>(&block[4]);
    const unsigned int res = hdr->res + 1;
    const unsigned int samples = num_of_samples[hdr->res];
    const float* const lut = dec_LUTMatrix[hdr->res];
    const uint32_t mask = (0x000000ffU >> (8U - res));
    const auto* samp_blk = reinterpret_cast<const uint32_t*>(&block[8]);

    std::size_t item_index = 0;
    for (std::size_t i = 0; i < 2 * SAMPLE_BLK_PER_CH_BLK; i += 2)
        {
            uint32_t sblk_q = samp_blk[i + 0];
            uint32_t sblk_i = samp_blk[i + 1];
            for (unsigned int j = 0; j < samples; j++)
                {
                    out[item_index] = gr_complex(lut[sblk_i & mask], lut[sblk_q & mask]);
                    sblk_q >>= res;
                    sblk_i >>= res;
                    item_index++;
                }
        }
    return item_index;
}
}  // namespace


// ---------------------------------------------------------------------------
// GDARSRawPortBuffer
// ---------------------------------------------------------------------------

GDARSRawPortBuffer::GDARSRawPortBuffer(std::size_t capacity)
    : capacity_(capacity)
{
}


void GDARSRawPortBuffer::push(const gr_complex* data, std::size_t count)
{
    std::unique_lock<std::mutex> lock(mutex_);
    std::size_t offset = 0;
    while (offset < count)
        {
            cv_.wait(lock, [this] { return buffer_.size() < capacity_ || abandoned_; });
            if (abandoned_)
                {
                    return;
                }
            const std::size_t room = capacity_ - buffer_.size();
            const std::size_t take = std::min(room, count - offset);
            buffer_.insert(buffer_.end(), data + offset, data + offset + take);
            offset += take;
            cv_.notify_all();
        }
}


std::size_t GDARSRawPortBuffer::pop(gr_complex* out, std::size_t max_count, std::size_t min_count, bool& eof_out)
{
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this, min_count] { return buffer_.size() >= min_count || input_done_ || abandoned_; });
    if (buffer_.empty())
        {
            eof_out = true;
            return 0;
        }
    eof_out = false;
    const std::size_t take = std::min(max_count, buffer_.size());
    std::copy(buffer_.begin(), buffer_.begin() + static_cast<std::ptrdiff_t>(take), out);
    buffer_.erase(buffer_.begin(), buffer_.begin() + static_cast<std::ptrdiff_t>(take));
    lock.unlock();
    cv_.notify_all();
    return take;
}


void GDARSRawPortBuffer::mark_input_done()
{
    {
        std::lock_guard<std::mutex> lock(mutex_);
        input_done_ = true;
    }
    cv_.notify_all();
}


void GDARSRawPortBuffer::abandon()
{
    {
        std::lock_guard<std::mutex> lock(mutex_);
        abandoned_ = true;
    }
    cv_.notify_all();
}


// ---------------------------------------------------------------------------
// GDARSDecodeEngine
// ---------------------------------------------------------------------------

GDARSDecodeEngine::GDARSDecodeEngine(const std::string& filename, const std::vector<uint32_t>& port_channel_ids,
    const std::vector<uint32_t>& channel_ids, const std::vector<double>& channel_sampling_frequencies)
{
    for (std::size_t g = 0; g < port_channel_ids.size(); g++)
        {
            const uint32_t id = port_channel_ids[g];
            if (id < ports_by_channel_.size())
                {
                    ports_by_channel_[id].push_back(g);
                }
            port_buffers_.push_back(std::make_unique<GDARSRawPortBuffer>(kPortBufferCapacity));
        }
    for (std::size_t i = 0; i < channel_ids.size(); i++)
        {
            const uint32_t id = channel_ids[i];
            if (id < channel_quantum_.size())
                {
                    channel_quantum_[id] = batch_quantum_from_sampling_frequency(channel_sampling_frequencies[i]);
                }
        }

    file_.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (file_.is_open())
        {
            std::cout << "GDARS_source: file '" << filename << "' successfully opened.\n";
        }
    else
        {
            std::cerr << "GDARS_source: failed to open file: " << filename << '\n';
            exit(1);
        }

    thread_ = std::thread(&GDARSDecodeEngine::run, this);
}


GDARSDecodeEngine::~GDARSDecodeEngine()
{
    stop_requested_ = true;
    for (auto& pb : port_buffers_)
        {
            pb->abandon();
        }
    if (thread_.joinable())
        {
            thread_.join();
        }
}


GDARSRawPortBuffer& GDARSDecodeEngine::port_buffer(std::size_t global_port)
{
    return *port_buffers_.at(global_port);
}


void GDARSDecodeEngine::run()
{
    name_current_thread("gdars_dec");
    std::vector<uint8_t> buffer(kReadChunkBytes);
    std::size_t buffer_valid = 0;
    std::size_t buffer_pos = 0;
    std::vector<gr_complex> decoded(kMaxSamplesPerBlock);

    // Per-channel-id decode accumulator: a block's decoded samples are
    // appended here, not pushed to that channel's ports immediately --
    // same rationale as GDARS_source::processor_loop()'s own batching,
    // applied one step earlier so push() (mutex lock + deque insert +
    // notify) isn't paid once per ~1-8k-sample raw block, which at
    // realistic capture rates can be tens of thousands of times a second.
    // Flushed once a channel's own accumulation reaches channel_quantum_
    // (computed from that channel's own sampling_frequency, same formula
    // GDARS_source uses), and unconditionally at EOF below so no decoded
    // tail is ever lost to a partially-filled accumulator.
    std::array<std::vector<gr_complex>, 8> pending;

    while (!stop_requested_.load(std::memory_order_relaxed))
        {
            if (buffer_pos + BYTES_PER_CH_BLK > buffer_valid)
                {
                    file_.read(reinterpret_cast<char*>(buffer.data()), static_cast<std::streamsize>(buffer.size()));
                    buffer_valid = static_cast<std::size_t>(file_.gcount());
                    buffer_pos = 0;
                    if (buffer_valid < BYTES_PER_CH_BLK)
                        {
                            break;  // EOF (or a short trailing read)
                        }
                }

            const uint8_t* block = &buffer[buffer_pos];
            const auto* hdr = reinterpret_cast<const block_header*>(&block[4]);
            buffer_pos += BYTES_PER_CH_BLK;

            if (hdr->id >= ports_by_channel_.size() || ports_by_channel_[hdr->id].empty())
                {
                    continue;  // disabled/unmapped channel -- discarded after reading just the header
                }

            const std::size_t n = decode_channel_block(decoded.data(), block);
            auto& acc = pending[hdr->id];
            acc.insert(acc.end(), decoded.begin(), decoded.begin() + static_cast<std::ptrdiff_t>(n));
            if (acc.size() >= channel_quantum_[hdr->id])
                {
                    for (std::size_t g : ports_by_channel_[hdr->id])
                        {
                            port_buffers_[g]->push(acc.data(), acc.size());
                        }
                    acc.clear();
                }
        }

    std::cout << "GDARS_source: end of file reached.\n";
    for (std::size_t id = 0; id < pending.size(); id++)
        {
            if (!pending[id].empty())
                {
                    for (std::size_t g : ports_by_channel_[id])
                        {
                            port_buffers_[g]->push(pending[id].data(), pending[id].size());
                        }
                }
        }
    for (auto& pb : port_buffers_)
        {
            pb->mark_input_done();
        }
}


// ---------------------------------------------------------------------------
// GDARS_source
// ---------------------------------------------------------------------------

GDARS_source_sptr GDARS_make_source(GDARSDecodeEngine_sptr engine, std::size_t global_port, double sampling_frequency,
    Concurrent_Queue<pmt::pmt_t>* queue, const GDARSFreqXlatingBand& band)
{
    return GDARS_source_sptr(new GDARS_source(std::move(engine), global_port, sampling_frequency, queue, band));
}


GDARS_source::GDARS_source(GDARSDecodeEngine_sptr engine, std::size_t global_port, double sampling_frequency,
    Concurrent_Queue<pmt::pmt_t>* queue, const GDARSFreqXlatingBand& band)
    : gr::sync_block("GDARS_source",
          gr::io_signature::make(0, 0, 0),
          gr::io_signature::make(1, 1, sizeof(gr_complex))),
      engine_(std::move(engine)),
      raw_buffer_(engine_->port_buffer(global_port)),
      global_port_(global_port),
      queue_(queue),
      batch_quantum_(batch_quantum_from_sampling_frequency(sampling_frequency)),
      out_buffer_(kPortBufferCapacity)
{
    // Same rationale as Evk1029Source's own set_min_output_buffer() call:
    // every acq+trk reader of every Channel mapped to this port (tens, at
    // realistic channel counts) connects directly to THIS block, so GNU
    // Radio's own output buffer for this port -- not out_buffer_ above,
    // which is this block's internal, pre-GNU-Radio staging buffer -- is
    // what those readers actually read from. Left at GNU Radio's small
    // default, that buffer forces work() to be called far more often with
    // much smaller batches to avoid overflow, which is a real throughput
    // cost (not a deadlock risk -- min_output_buffer is the lever
    // Evk1029Source itself documents as the safe one). This was missed in
    // the port-per-instance rewrite and measured directly: with it unset,
    // Rec/Wall stayed capped around 0.74x on a config EVK1029 runs at
    // 2.1x despite neither the GDARS threads nor the CPU as a whole being
    // saturated -- i.e. exactly this class of symptom.
    set_min_output_buffer(0, 4000000);

    if (band.enable_freq_xlating && !band.taps.empty())
        {
            direct_ = false;
            decimation_ = band.decimation_factor;
            composite_fir_.resize(band.taps.size());
            const float fwT0 = static_cast<float>(2.0 * GR_M_PI * band.center_freq / sampling_frequency);
            for (std::size_t i = 0; i < band.taps.size(); i++)
                {
                    composite_fir_[i] = band.taps[i] * std::exp(gr_complex(0, static_cast<float>(i) * fwT0));
                }
            std::reverse(composite_fir_.begin(), composite_fir_.end());  // match gr::filter::kernel::fir_filter::set_taps()'s internal reversal
            history_len_ = composite_fir_.size() - 1;
            phase_incr_ = std::exp(gr_complex(0, -fwT0 * static_cast<float>(decimation_)));
            const float group_delay_phase = -fwT0 * static_cast<float>(band.taps.size() - 1) / 2.0F;
            phase_ = std::exp(gr_complex(0, group_delay_phase));
#if CUDA_GPU_ACCEL
            if (band.use_cuda)
                {
                    auto cuda_engine = std::make_unique<CudaDdcEngine>();
                    if (cuda_engine->is_valid() && cuda_engine->set_taps(reinterpret_cast<const std::complex<float>*>(composite_fir_.data()), static_cast<int>(composite_fir_.size())))
                        {
                            std::cout << "GDARS_source: RF output " << global_port_ << " using CUDA device " << cuda_engine->device_name() << '\n';
                            cuda_engine_ = std::move(cuda_engine);
                        }
                    else
                        {
                            std::cerr << "GDARS_source: RF output " << global_port_ << ": CUDA DDC engine unavailable (" << cuda_engine->last_error() << "), falling back to CPU\n";
                        }
                }
#endif
        }
    else
        {
            direct_ = true;
            decimation_ = 1;
            history_len_ = 0;
        }

    processor_thread_ = std::thread(&GDARS_source::processor_loop, this);
}


GDARS_source::~GDARS_source()
{
    stop_ = true;
    out_buffer_.abandon();
    raw_buffer_.abandon();
    if (processor_thread_.joinable())
        {
            processor_thread_.join();
        }
}


bool GDARS_source::stop()
{
    stop_ = true;
    out_buffer_.abandon();
    raw_buffer_.abandon();
    return true;
}


void GDARS_source::processor_loop()
{
    name_current_thread("gdars_out" + std::to_string(global_port_));
    // Pop at least batch_quantum_ raw samples per iteration (see its own
    // doc comment) -- max_count is the same size so one pop() drains a
    // full accumulated batch in one call rather than needing several.
    const std::size_t kChunk = batch_quantum_;
    std::vector<gr_complex> new_samples(kChunk);

    while (!stop_.load(std::memory_order_relaxed))
        {
            bool eof = false;
            const std::size_t n = raw_buffer_.pop(new_samples.data(), kChunk, batch_quantum_, eof);
            if (n == 0)
                {
                    if (eof)
                        {
                            out_buffer_.mark_input_done();
                            return;
                        }
                    continue;
                }

            if (direct_)
                {
                    out_buffer_.push(new_samples.data(), n);
                    continue;
                }

            raw_window_.resize(history_len_ + n);
            std::copy(new_samples.begin(), new_samples.begin() + static_cast<std::ptrdiff_t>(n), raw_window_.begin() + static_cast<std::ptrdiff_t>(history_len_));

            const std::size_t ntaps = composite_fir_.size();
            const std::size_t usable = (n - 1) / static_cast<std::size_t>(decimation_) + 1;
            std::vector<gr_complex> produced(usable);

            bool computed_on_gpu = false;
#if CUDA_GPU_ACCEL
            if (cuda_engine_)
                {
                    const int n_data = static_cast<int>((usable - 1) * static_cast<std::size_t>(decimation_) + ntaps);
                    computed_on_gpu = cuda_engine_->compute_complex(reinterpret_cast<const std::complex<float>*>(raw_window_.data()), n_data, decimation_,
                        reinterpret_cast<std::complex<float>*>(produced.data()), static_cast<int>(usable));
                    if (!computed_on_gpu)
                        {
                            std::cerr << "GDARS_source: RF output " << global_port_ << ": CUDA DDC compute failed (" << cuda_engine_->last_error() << "), falling back to CPU\n";
                            cuda_engine_.reset();  // stay on CPU for every subsequent iteration too
                        }
                }
#endif
            if (!computed_on_gpu)
                {
                    for (std::size_t i = 0; i < usable; i++)
                        {
                            volk_32fc_x2_dot_prod_32fc(&produced[i], composite_fir_.data(), raw_window_.data() + i * static_cast<std::size_t>(decimation_), ntaps);
                        }
                }

            volk_32fc_s32fc_x2_rotator2_32fc(produced.data(), produced.data(), &phase_incr_, &phase_, static_cast<int>(usable));
            samples_since_phase_renorm_ += static_cast<int>(usable);
            if (samples_since_phase_renorm_ > 4096)
                {
                    phase_ /= std::abs(phase_);
                    samples_since_phase_renorm_ = 0;
                }

            out_buffer_.push(produced.data(), produced.size());

            // carry the trailing history_len_ raw samples forward for the
            // next iteration's FIR continuity (overlap-save)
            raw_window_.erase(raw_window_.begin(), raw_window_.end() - static_cast<std::ptrdiff_t>(history_len_));
        }
}


int GDARS_source::work(int noutput_items,
    gr_vector_const_void_star& input_items __attribute__((unused)),
    gr_vector_void_star& output_items)
{
    // min_count=1: GNU Radio's own noutput_items already encodes how much
    // downstream wants and has room for -- an artificial minimum here
    // would only add latency to data that's already ready, unlike
    // processor_loop()'s batching of raw_buffer_, which exists purely to
    // amortize FIR/copy work and has no such downstream-readiness signal
    // to respect.
    bool eof = false;
    const std::size_t n = out_buffer_.pop(reinterpret_cast<gr_complex*>(output_items[0]), static_cast<std::size_t>(noutput_items), 1, eof);
    if (n == 0 && eof)
        {
            std::cout << "GDARS_source: RF output " << global_port_ << ": EOF\n";
            queue_->push(pmt::make_any(command_event_make(200, 0)));
            return this->WORK_DONE;
        }
    return static_cast<int>(n);
}

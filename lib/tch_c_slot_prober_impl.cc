/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tch_c_format.h"
#include "tch_c_slot_prober_impl.h"
#include <gnuradio/block_detail.h>
#include <gnuradio/buffer_reader.h>
#include <gnuradio/io_signature.h>
#include <gnuradio/tpb_detail.h>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

namespace gr {
namespace iridium {

// The same limit as burst_ack_counter::MAX_BACKLOG: the frame sorter needs the frames
// within 2 s
static const uint64_t PROBE_MAX_BACKLOG = 1000000000; // ns

// Acknowledged probes are forgotten this long (of the time of the bursts) after the
// newest one: tch_c_demod can report the frame after the acknowledgement
static const uint64_t PROBE_KEEP = 3000000000; // ns

// Input appended to the ring buffer per step, and the ring buffer holds this much more
// than the history: one step, the longest probe window and some slack
static const double CHUNK_S = 0.005;
static const double RING_MARGIN_S = 0.03;

// Search range of the survey around the predicted carrier of a stale track (no
// measurement for a while); otherwise +-1 bin
static const double STALE_SEARCH_HZ = 3000;

// Two tracks probe the same slot (one satellite) if their predictions are this close
static const double SAME_SLOT_NS = 20000;
static const double SAME_SLOT_HZ = 2000;

/*
 * Default thresholds of rho2 per type. On sat-147 (full file, 95,000 probes) below them
 * hardly a probe passes the CRC: TCH_C1 rho2 0.35-0.6 7 of 36,000, TCH_C2 0.35-0.55 0
 * of 910, TCH_C8 0.35-0.55 1 of 2,760; just above them 1-2% do.
 */
static const float DEFAULT_THRESHOLD[3] = { 0.6f, 0.55f, 0.55f };

// tch_c_demod tries all its carrier recovery variants on probes with at least this
// rho2, only the fit to the known symbols on the others. Every variant is one more
// chance for a failed decode to pass the CRC-10, and below 0.7 most of the decodes fail
// (sat-147: 0.6-0.7 the other variants add 36 frames to 391 for about 11,000 more
// failed decodes, above 0.7 760 to 38,000 for 14,000).
static const float VARIANTS_RHO2 = 0.7f;

// Strong hits (rho2 of the survey) update the timing bias of their type
static const float BIAS_MIN_RHO2 = 0.6f;

/*
 * Uplink. The UL slots are probed this long after the DL slots, so the DL frames of the
 * frames before have come back from tch_c_demod. A channel is searched in a UL slot if
 * its DL slot carried a frame within UL_DL_ACTIVE_NS, widely (UL_WIDE_*: around the
 * nominal slot time from W0 and the carrier nominal * (1 - Doppler / F_RING) + 0..5.5
 * kHz seen on sat-147) until a UL frame was found within UL_TRACK_NS, then around that
 * one. UL_THRESHOLD: on sat-147 the UL bursts have rho2 0.8-1 (all CRC ok); the wide
 * search over noise reaches 0.4 in 1 of 10,000 channels, 0.45 never.
 */
static const double UL_DELAY_S = 0.25;
static const int64_t UL_DL_ACTIVE_NS = 3000000000;
static const int64_t UL_TRACK_NS = 5000000000;
static const double UL_WIDE_S = 350e-6;
static const double UL_WIDE_F_HZ = 2750;
static const double UL_WIDE_SEARCH_HZ = 3750;
static const double UL_TRACK_S = 50e-6;
static const double UL_TRACK_SEARCH_HZ = 1500;
static const float UL_THRESHOLD = 0.45f;
// UL frames of the energy path: tolerance of their carrier for the channel match
static const double UL_MATCH_TOL_HZ = 8000;
// Acquisition sweep while W0 is unknown: half width of a window, and the slow down
static const double UL_SWEEP_S = 1.9e-3;
static const int UL_SWEEP_WINDOWS = 2;
static const uint64_t UL_SWEEP_FRAMES = 600;
static const uint64_t UL_SWEEP_SLOW = 8;
// W0 from the median of this many UL frames of the last UL_W_WINDOW_NS
static const size_t UL_W_MEASUREMENTS = 32;
static const int64_t UL_W_WINDOW_NS = 30000000000;
// No frame yet (far from any time, no overflow in differences)
static const int64_t UL_NEVER = INT64_MIN / 4;

tch_c_slot_prober::sptr tch_c_slot_prober::make(double center_frequency,
                                                double sample_rate,
                                                const std::vector<int>& types,
                                                float history_s,
                                                float threshold,
                                                int outstanding_limit,
                                                bool drop_overflow,
                                                bool uplink)
{
    return gnuradio::make_block_sptr<tch_c_slot_prober_impl>(center_frequency,
                                                             sample_rate,
                                                             types,
                                                             history_s,
                                                             threshold,
                                                             outstanding_limit,
                                                             drop_overflow,
                                                             uplink);
}

// Thread CPU time in s, for the statistics of IRIDIUM_PROBE_DEBUG
static double cpu_s()
{
    timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static bool env_flag(const char* name)
{
    const char* v = std::getenv(name);
    return v && *v && strcmp(v, "0") != 0;
}

static double env_double(const char* name, double fallback)
{
    const char* v = std::getenv(name);
    return v && *v ? atof(v) : fallback;
}

// Raster position of a type nearest to a nominal carrier
static int nearest_position(int type, double nominal_hz)
{
    const double k = type == 8 ? 8 : type == 2 ? 2 : 1;
    const double o = type == 8 ? 4 : type == 2 ? 1 : 0.5;
    const double c = (nominal_hz - tch_c_slot_survey::BAND_START) /
                     tch_c_slot_survey::CHANNEL_WIDTH;
    const int pos = (int)std::lround((c - o) / k);
    return std::max(0, std::min(tch_c_slot_survey::n_positions(type) - 1, pos));
}

tch_c_slot_prober_impl::tch_c_slot_prober_impl(double center_frequency,
                                               double sample_rate,
                                               const std::vector<int>& types,
                                               float history_s,
                                               float threshold,
                                               int outstanding_limit,
                                               bool drop_overflow,
                                               bool uplink)
    : gr::sync_block("tch_c_slot_prober",
                     gr::io_signature::make(1, 1, sizeof(gr_complex)),
                     gr::io_signature::make(0, 0, 0)),
      d_center_frequency(center_frequency),
      d_sample_rate(sample_rate),
      d_types(types),
      d_history_s(history_s),
      d_threshold{ threshold > 0 ? threshold : DEFAULT_THRESHOLD[0],
                   threshold > 0 ? threshold : DEFAULT_THRESHOLD[1],
                   threshold > 0 ? threshold : DEFAULT_THRESHOLD[2] },
      d_outstanding_limit(outstanding_limit),
      d_drop_overflow(drop_overflow),
      d_uplink(uplink),
      d_debug(env_flag("IRIDIUM_PROBE_DEBUG")),
      d_debug_rho2(d_debug ? env_double("IRIDIUM_PROBE_DEBUG_RHO2", 2) : 2),
      d_offset_ns(std::llround(env_double("IRIDIUM_PROBE_OFFSET_US", 0) * 1e3)),
      d_variants_rho2(env_double("IRIDIUM_PROBE_VARIANTS_RHO2", VARIANTS_RHO2)),
      d_ring(std::max<size_t>(
          std::ceil((history_s + RING_MARGIN_S + (uplink ? UL_DELAY_S : 0)) *
                    sample_rate),
          1)),
      d_end(0),
      d_chunk(std::max<uint64_t>(1, CHUNK_S * sample_rate)),
      d_lag(std::max<double>(0, history_s) * sample_rate),
      d_flushed(false),
      // Not in the false alarm control (shifted probe times): the acquisition sweep
      // searches all of the frame and would find the real UL bursts
      d_probe_uplink(uplink && d_offset_ns == 0 &&
                     std::find(types.begin(), types.end(), 1) != types.end()),
      d_ul_w_valid(false),
      d_ul_w0(0),
      d_n_ul_sweep_frames(0),
      d_n_ul_slots(0),
      d_n_ul_searches(0),
      d_n_ul_probes(0),
      d_n_ul_crc_ok(0),
      d_survey(sample_rate, center_frequency),
      d_symbols(TCH_C_MAX_SYMBOLS + 1),
      d_next_id(TCH_C_PROBE_ID_BASE),
      d_outstanding(0),
      d_newest(0),
      d_n_probes(0),
      d_n_dropped(0),
      d_n_iras(0),
      d_n_hits(0),
      d_n_probe_hits(0),
      d_n_probe_crc_ok(0),
      d_n_energy_associated(0),
      d_n_missed_slots(0),
      d_n_shared_slots(0),
      d_ira_delay_sum(0),
      d_ira_delay_max(0)
{
    for (int i = 0; i < 3; i++) {
        d_n_surveys[i] = d_n_type_probes[i] = d_n_type_crc_ok[i] = 0;
    }
    for (int i = 0; i < 4; i++) {
        d_cpu[i] = 0;
    }

    message_port_register_in(pmt::mp("frames"));
    set_msg_handler(pmt::mp("frames"),
                    [this](const pmt::pmt_t& msg) { this->handle_frame(msg); });
    message_port_register_in(pmt::mp("hits"));
    set_msg_handler(pmt::mp("hits"),
                    [this](const pmt::pmt_t& msg) { this->handle_hit(msg); });
    message_port_register_in(pmt::mp("frame_handled"));
    set_msg_handler(pmt::mp("frame_handled"),
                    [this](const pmt::pmt_t& msg) { this->handle_frame_handled(msg); });
    message_port_register_out(pmt::mp("cpdus"));

    /*
     * The message neighbours (demodulators) report "done" when they finish. GNU Radio
     * then ends this block at the next work() call which consumes nothing, which is
     * what work() does while it waits for acknowledgements, and the rest of the stream
     * would be lost. This block ends with its stream input.
     */
    set_msg_handler(pmt::mp("system"), [this](const pmt::pmt_t& msg) {
        if (pmt::is_pair(msg) && pmt::eqv(pmt::car(msg), pmt::mp("done"))) {
            return;
        }
        this->system_handler(msg);
    });
}

tch_c_slot_prober_impl::~tch_c_slot_prober_impl() {}

bool tch_c_slot_prober_impl::stop()
{
    if (d_debug) {
        fprintf(stderr,
                "tch_c_slot_prober: %" PRIu64 " ring alerts (delay mean %.3f s, max %.3f "
                "s), %" PRIu64 " probes (%" PRIu64 " dropped), %" PRIu64
                " probe frames (%" PRIu64 " CRC ok), %" PRIu64
                " frames from the energy path (%" PRIu64 " associated), %d outstanding\n",
                d_n_iras,
                d_n_iras ? d_ira_delay_sum / d_n_iras : 0.,
                d_ira_delay_max,
                d_n_probes,
                d_n_dropped,
                d_n_probe_hits,
                d_n_probe_crc_ok,
                d_n_hits - d_n_probe_hits,
                d_n_energy_associated,
                d_outstanding);
        fprintf(stderr,
                "tch_c_slot_prober: tracks %" PRIu64 " created, %" PRIu64
                " merged, %" PRIu64 " expired, %" PRIu64
                " evicted; slots %" PRIu64 " missed, %" PRIu64 " shared\n",
                d_tracks.n_created(),
                d_tracks.n_merged(),
                d_tracks.n_expired(),
                d_tracks.n_evicted(),
                d_n_missed_slots,
                d_n_shared_slots);
        fprintf(stderr,
                "tch_c_slot_prober: CPU %.2f s in work(): survey %.2f s, fine timing "
                "and extraction %.2f s, uplink %.2f s\n",
                d_cpu[2],
                d_cpu[0],
                d_cpu[1],
                d_cpu[3]);
        for (int i = 0; i < 3; i++) {
            fprintf(stderr,
                    "tch_c_slot_prober: TCH_C%d %" PRIu64 " slot surveys, %" PRIu64
                    " probes, %" PRIu64 " CRC ok\n",
                    i == 2 ? 8 : i + 1,
                    d_n_surveys[i],
                    d_n_type_probes[i],
                    d_n_type_crc_ok[i]);
        }
        if (d_probe_uplink) {
            fprintf(stderr,
                    "tch_c_slot_prober: TCH_C1 uplink W0 %.3f ms (%s), %" PRIu64
                    " slots, %" PRIu64 " sweep windows, %" PRIu64
                    " channel searches, %" PRIu64 " probes, %" PRIu64 " CRC ok\n",
                    d_ul_w0,
                    d_ul_w_valid ? "valid" : "not found",
                    d_n_ul_slots,
                    d_n_ul_sweep_frames,
                    d_n_ul_searches,
                    d_n_ul_probes,
                    d_n_ul_crc_ok);
        }
    }
    return true;
}

/*
 * Time base
 */
double tch_c_slot_prober_impl::sample_index(uint64_t t_ns) const
{
    if (d_time_bases.empty()) {
        return (double)t_ns * d_sample_rate / 1e9;
    }
    // The newest time base whose first sample is not after t
    const time_base* b = &d_time_bases.front();
    for (auto it = d_time_bases.rbegin(); it != d_time_bases.rend(); ++it) {
        const double t0 = it->timestamp + ((double)it->start - (double)it->offset) * 1e9 /
                                              d_sample_rate;
        if (t0 <= (double)t_ns) {
            b = &*it;
            break;
        }
    }
    return b->offset + (double)(int64_t)(t_ns - b->timestamp) * d_sample_rate / 1e9;
}

uint64_t tch_c_slot_prober_impl::time_ns(double sample_index) const
{
    if (d_time_bases.empty()) {
        return (uint64_t)std::floor(sample_index * 1e9 / d_sample_rate);
    }
    const time_base* b = &d_time_bases.front();
    for (auto it = d_time_bases.rbegin(); it != d_time_bases.rend(); ++it) {
        if ((double)it->start <= sample_index) {
            b = &*it;
            break;
        }
    }
    return b->timestamp +
           (int64_t)std::floor((sample_index - b->offset) * 1e9 / d_sample_rate);
}

void tch_c_slot_prober_impl::update_time_base(uint64_t start, uint64_t end)
{
    std::vector<tag_t> tags;
    get_tags_in_range(tags, 0, start, end, pmt::mp("time_base"));
    std::sort(tags.begin(), tags.end(), tag_t::offset_compare);
    for (const auto& tag : tags) {
        d_time_bases.push_back({ tag.offset,
                                 pmt::to_uint64(pmt::tuple_ref(tag.value, 0)),
                                 pmt::to_uint64(pmt::tuple_ref(tag.value, 1)) });
    }
    // Keep those still needed for the samples in the buffer
    const uint64_t begin = ring_begin();
    size_t drop = 0;
    while (drop + 1 < d_time_bases.size() && d_time_bases[drop + 1].start <= begin) {
        drop++;
    }
    d_time_bases.erase(d_time_bases.begin(), d_time_bases.begin() + drop);
}

/*
 * Ring buffer
 */
uint64_t tch_c_slot_prober_impl::ring_begin() const
{
    return d_end > d_ring.size() ? d_end - d_ring.size() : 0;
}

const gr_complex* tch_c_slot_prober_impl::samples(uint64_t start, size_t len)
{
    assert(start >= ring_begin() && start + len <= d_end);
    const size_t size = d_ring.size();
    const size_t i = start % size;
    if (i + len <= size) {
        return &d_ring[i];
    }
    d_scratch.resize(len);
    const size_t first = size - i;
    memcpy(&d_scratch[0], &d_ring[i], first * sizeof(gr_complex));
    memcpy(&d_scratch[first], &d_ring[0], (len - first) * sizeof(gr_complex));
    return &d_scratch[0];
}

void tch_c_slot_prober_impl::append(const gr_complex* in, size_t n)
{
    const size_t size = d_ring.size();
    if (n > size) {
        in += n - size;
        d_end += n - size;
        n = size;
    }
    const size_t i = d_end % size;
    const size_t first = std::min(n, size - i);
    memcpy(&d_ring[i], in, first * sizeof(gr_complex));
    memcpy(&d_ring[0], in + first, (n - first) * sizeof(gr_complex));
    d_end += n;
}

/*
 * Probes
 */
bool tch_c_slot_prober_impl::blocked() const
{
    if (d_outstanding_limit > 0 && d_outstanding >= d_outstanding_limit) {
        return true;
    }
    return !d_unacked_times.empty() && d_newest > *d_unacked_times.begin() &&
           d_newest - *d_unacked_times.begin() > PROBE_MAX_BACKLOG;
}

uint64_t tch_c_slot_prober_impl::publish_probe(int type,
                                               uint64_t t_first_uw_ns,
                                               double freq_hz,
                                               const gr_complex* symbols,
                                               size_t n,
                                               float magnitude,
                                               float noise,
                                               float rho2,
                                               int track,
                                               int position,
                                               int slot,
                                               bool uplink)
{
    // Reading from a file work() waits instead
    if (d_drop_overflow && blocked()) {
        d_n_dropped++;
        return 0;
    }

    const tch_c_format& f = tch_c_format_by_type(type);
    const uint64_t id = d_next_id;
    d_next_id += 10;

    pmt::pmt_t meta = pmt::make_dict();
    meta = pmt::dict_add(meta, pmt::mp("timestamp"), pmt::from_uint64(t_first_uw_ns));
    meta = pmt::dict_add(meta, pmt::mp("uw_start"), pmt::from_float(0));
    meta = pmt::dict_add(meta, pmt::mp("center_frequency"), pmt::from_double(freq_hz));
    meta = pmt::dict_add(
        meta, pmt::mp("sample_rate"), pmt::from_float(f.symbols_per_second));
    meta = pmt::dict_add(meta, pmt::mp("tch_type"), pmt::from_long(type));
    meta = pmt::dict_add(meta, pmt::mp("id"), pmt::from_uint64(id));
    meta = pmt::dict_add(meta, pmt::mp("burst_id"), pmt::from_uint64(id));
    meta = pmt::dict_add(meta, pmt::mp("noise"), pmt::from_float(noise));
    meta = pmt::dict_add(meta, pmt::mp("magnitude"), pmt::from_float(magnitude));
    meta = pmt::dict_add(meta, pmt::mp("rho2"), pmt::from_float(rho2));
    meta = pmt::dict_add(meta, pmt::mp("probe"), pmt::PMT_T);
    meta = pmt::dict_add(
        meta, pmt::mp("variants"), pmt::from_long(rho2 >= d_variants_rho2 ? 3 : 1));

    d_probes[id] = { type,     t_first_uw_ns, freq_hz, rho2, track,
                     position, slot,          uplink,  false };
    d_unacked_times.insert(t_first_uw_ns);
    d_outstanding++;
    d_newest = std::max(d_newest, t_first_uw_ns);
    d_n_probes++;
    if (uplink) {
        d_n_ul_probes++;
    } else {
        d_n_type_probes[tch_c_timing::type_index(type)]++;
    }

    message_port_pub(pmt::mp("cpdus"), pmt::cons(meta, pmt::init_c32vector(n, symbols)));
    return id;
}

void tch_c_slot_prober_impl::handle_frame_handled(const pmt::pmt_t& msg)
{
    auto it = d_probes.find(pmt::to_uint64(msg));
    if (it == d_probes.end() || it->second.acked) {
        return;
    }
    it->second.acked = true;
    d_unacked_times.erase(d_unacked_times.find(it->second.t_ns));
    d_outstanding--;
}

void tch_c_slot_prober_impl::forget_old_probes()
{
    for (auto it = d_probes.begin(); it != d_probes.end();) {
        if (it->second.acked && it->second.t_ns + PROBE_KEEP < d_newest) {
            it = d_probes.erase(it);
        } else {
            ++it;
        }
    }
}

/*
 * Ring alerts: downlink frames on the simplex ring channel, long and clean enough to
 * be ring alerts. They need not be decoded.
 */
void tch_c_slot_prober_impl::handle_frame(const pmt::pmt_t& msg)
{
    const pmt::pmt_t meta = pmt::car(msg);
    if (pmt::to_long(pmt::dict_ref(meta, pmt::mp("direction"), pmt::from_long(-1))) !=
        0) {
        return;
    }
    const double f =
        pmt::to_double(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));
    if (std::abs(f - TCH_C_F_RING) > 45000) {
        return;
    }
    if (pmt::to_long(pmt::dict_ref(meta, pmt::mp("n_symbols"), pmt::from_long(0))) < 80 ||
        pmt::to_long(pmt::dict_ref(meta, pmt::mp("confidence"), pmt::from_long(0))) <
            80) {
        return;
    }
    const uint64_t t =
        pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("timestamp"), pmt::PMT_NIL));
    d_iras.emplace_back(t, f);
    d_n_iras++;
}

/*
 * TCH_C frames from tch_c_demod: our probes (ids from TCH_C_PROBE_ID_BASE) and those of
 * the energy path. tch_c_demod only reports frames with a UW, CRC ok or not.
 */
void tch_c_slot_prober_impl::handle_hit(const pmt::pmt_t& msg)
{
    const pmt::pmt_t meta = pmt::car(msg);
    tch_c_hit hit;
    hit.id = pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("id"), pmt::PMT_NIL));
    hit.t_ns = pmt::to_uint64(pmt::dict_ref(meta, pmt::mp("timestamp"), pmt::PMT_NIL));
    hit.f_hz =
        pmt::to_double(pmt::dict_ref(meta, pmt::mp("center_frequency"), pmt::PMT_NIL));
    hit.type = pmt::to_long(pmt::dict_ref(meta, pmt::mp("tch_type"), pmt::PMT_NIL));
    hit.bits_per_symbol =
        pmt::to_long(pmt::dict_ref(meta, pmt::mp("bits_per_symbol"), pmt::from_long(2)));
    hit.uplink = pmt::to_bool(pmt::dict_ref(meta, pmt::mp("uplink"), pmt::PMT_F));
    hit.crc_ok = pmt::to_bool(pmt::dict_ref(meta, pmt::mp("crc_ok"), pmt::PMT_F));
    hit.confidence =
        pmt::to_long(pmt::dict_ref(meta, pmt::mp("confidence"), pmt::from_long(0)));
    hit.level = pmt::to_float(pmt::dict_ref(meta, pmt::mp("level"), pmt::from_float(0)));
    hit.magnitude =
        pmt::to_float(pmt::dict_ref(meta, pmt::mp("magnitude"), pmt::from_float(0)));
    hit.noise = pmt::to_float(pmt::dict_ref(meta, pmt::mp("noise"), pmt::from_float(0)));
    hit.variant =
        pmt::to_long(pmt::dict_ref(meta, pmt::mp("variant"), pmt::from_long(0)));
    hit.probe = nullptr;

    d_n_hits++;
    if (hit.id >= TCH_C_PROBE_ID_BASE) {
        auto it = d_probes.find(hit.id);
        if (it != d_probes.end()) {
            hit.probe = &it->second;
        }
        d_n_probe_hits++;
        if (hit.crc_ok) {
            d_n_probe_crc_ok++;
            if (hit.uplink) {
                d_n_ul_crc_ok++;
            } else {
                d_n_type_crc_ok[tch_c_timing::type_index(hit.type)]++;
            }
        }
        if (d_debug) {
            fprintf(stderr,
                    "PH %" PRIu64 " crc=%d uplink=%d variant=%d bps=%d df=%.1f\n",
                    hit.id,
                    hit.crc_ok,
                    hit.uplink,
                    hit.variant,
                    hit.bits_per_symbol,
                    hit.probe ? hit.f_hz - hit.probe->f_hz : 0.);
        }
    }

    on_hit(hit);
}

/*
 * A ring alert at t_ns (time of its first UW symbol) on f_hz: a measurement of the
 * frame phase and Doppler of its satellite, or a new satellite.
 */
void tch_c_slot_prober_impl::on_ira(uint64_t t_ns, double f_hz)
{
    const double delay = ((double)time_ns(d_end) - (double)t_ns) * 1e-9;
    d_ira_delay_sum += delay;
    d_ira_delay_max = std::max(d_ira_delay_max, delay);
    const int id = d_tracks.add_ira(t_ns, f_hz);
    if (d_debug) {
        fprintf(stderr,
                "IRA t_ns=%" PRIu64 " f=%.1f delay=%.3f track=%d\n",
                t_ns,
                f_hz,
                delay,
                id);
    }
}

/*
 * A TCH_C frame. CRC ok downlink frames update the track table: our own probes the
 * track and slot they were made for (strong ones also the timing bias), those of the
 * energy path the track whose slot they fit.
 */
void tch_c_slot_prober_impl::on_hit(const tch_c_hit& hit)
{
    if (!hit.crc_ok) {
        return;
    }
    if (hit.uplink) {
        on_ul_hit(hit);
        return;
    }
    if (hit.id >= TCH_C_PROBE_ID_BASE) {
        // In the control run (shifted probe times) the probes are not where the
        // bursts are
        if (!hit.probe || hit.probe->track < 0 || d_offset_ns != 0) {
            return;
        }
        const tch_c_probe_info& p = *hit.probe;
        if (p.type == 1) {
            ul_link_dl(p.track, p.slot, p.position, hit.t_ns);
        }
        d_tracks.add_hit(p.track,
                         hit.t_ns,
                         hit.f_hz,
                         p.type,
                         p.slot,
                         tch_c_slot_survey::nominal_frequency(p.type, p.position),
                         p.rho2 >= BIAS_MIN_RHO2);
        return;
    }

    // Energy path: the raster channel depends on the Doppler of the track
    int best = -1;
    tch_c_track_table::slot_match best_match;
    double best_nominal = 0;
    for (const auto& tr : d_tracks) {
        if (tr.empty()) {
            continue;
        }
        const double nominal = tch_c_slot_survey::nominal_frequency(
            hit.type,
            nearest_position(hit.type,
                             hit.f_hz / (1 + tr.doppler(hit.t_ns) / TCH_C_F_RING)));
        tch_c_track_table::slot_match m;
        if (d_tracks.match_hit(hit.t_ns, hit.f_hz, hit.type, nominal, m) &&
            m.track == tr.id() &&
            (best < 0 || std::fabs(m.dt_ns) < std::fabs(best_match.dt_ns))) {
            best = tr.id();
            best_match = m;
            best_nominal = nominal;
        }
    }
    if (best >= 0) {
        if (hit.type == 1) {
            ul_link_dl(
                best, best_match.slot, nearest_position(1, best_nominal), hit.t_ns);
        }
        d_tracks.add_hit(
            best, hit.t_ns, hit.f_hz, hit.type, best_match.slot, best_nominal, false);
        d_n_energy_associated++;
    }
}

void tch_c_slot_prober_impl::ul_link_dl(int track, int slot, int position, int64_t t_ns)
{
    if (!d_probe_uplink) {
        return;
    }
    auto it = d_ul_links.find(std::make_tuple(track, slot, position));
    if (it == d_ul_links.end()) {
        d_ul_links[std::make_tuple(track, slot, position)] = { t_ns, UL_NEVER, 0, 0 };
    } else {
        it->second.last_dl_ns = std::max(it->second.last_dl_ns, t_ns);
    }
}

/*
 * Probing
 */
void tch_c_slot_prober_impl::slot_window(int type,
                                         double p_uw,
                                         int64_t* first,
                                         int64_t* end) const
{
    // Survey window with the margin of the fine timing, and the extraction window
    // for any offset the fine timing can find
    const int m = d_survey.timing_margin(type);
    const int64_t i0 = std::lround(p_uw) - d_survey.uw_offset(type);
    long ef, ee, ef2, ee2;
    d_survey.extract_span(type, p_uw - m, &ef, &ee);
    d_survey.extract_span(type, p_uw + m, &ef2, &ee2);
    *first = std::min<int64_t>(i0 - m, ef);
    *end = std::max<int64_t>(i0 + d_survey.window_length(type) + m, ee2);
}

bool tch_c_slot_prober_impl::already_probed(int type,
                                            int64_t t_ns,
                                            double doppler_hz) const
{
    for (const auto& p : d_probed) {
        if (p.type == type && std::fabs((double)(p.t_ns - t_ns)) < SAME_SLOT_NS &&
            std::fabs(p.doppler_hz - doppler_hz) < SAME_SLOT_HZ) {
            return true;
        }
    }
    return false;
}

/*
 * Every track has a cursor per type: the next DL slot to probe. A slot is probed once
 * its window lies in the buffer and starts before `limit`. work() keeps limit
 * history_s behind the input, as far as the buffer allows: by then the ring alerts of
 * the following history_s have arrived too (they come late, behind the legacy
 * demodulators), so the slot time is interpolated rather than extrapolated and a
 * satellite found by its first ring alert is probed back to the start of the buffer
 * (a new track starts its cursors there). As all cursors move together, merging two
 * tracks loses no slot; a slot whose window has left the buffer is skipped (missed).
 */
void tch_c_slot_prober_impl::probe_slots(uint64_t limit, uint64_t ul_limit)
{
    const uint64_t begin = ring_begin();
    const int64_t begin_ns = time_ns(begin);
    d_tracks.expire(time_ns(std::min(limit, d_end)));

    // Forget the cursors of tracks gone and the slots which left the buffer
    for (auto it = d_cursors.begin(); it != d_cursors.end();) {
        if (!d_tracks.find(it->first)) {
            it = d_cursors.erase(it);
        } else {
            ++it;
        }
    }
    for (auto it = d_ul_cursors.begin(); it != d_ul_cursors.end();) {
        if (!d_tracks.find(it->first)) {
            it = d_ul_cursors.erase(it);
        } else {
            ++it;
        }
    }
    for (auto it = d_ul_sweeps.begin(); it != d_ul_sweeps.end();) {
        if (!d_tracks.find(it->first)) {
            it = d_ul_sweeps.erase(it);
        } else {
            ++it;
        }
    }
    // UL links of tracks gone or without frames for a while
    for (auto it = d_ul_links.begin(); it != d_ul_links.end();) {
        const int64_t last = std::max(it->second.last_dl_ns, it->second.last_ul_ns);
        if (!d_tracks.find(std::get<0>(it->first)) ||
            last < begin_ns - std::max(UL_DL_ACTIVE_NS, UL_TRACK_NS)) {
            it = d_ul_links.erase(it);
        } else {
            ++it;
        }
    }
    d_probed.erase(std::remove_if(d_probed.begin(),
                                  d_probed.end(),
                                  [&](const probed_slot& p) {
                                      return p.t_ns < begin_ns - 100000000;
                                  }),
                   d_probed.end());

    // Confirmed tracks first, then the older ones: if two tracks of one satellite
    // exist, the better one probes
    std::vector<const tch_c_sat_track*> order;
    for (const auto& tr : d_tracks) {
        if (!tr.empty()) {
            order.push_back(&tr);
        }
    }
    std::sort(order.begin(),
              order.end(),
              [](const tch_c_sat_track* a, const tch_c_sat_track* b) {
                  if (a->confirmed() != b->confirmed()) {
                      return a->confirmed();
                  }
                  return a->id() < b->id();
              });

    for (const tch_c_sat_track* tr : order) {
        std::vector<cursor>& cursors = d_cursors[tr->id()];
        cursors.resize(d_types.size(), { 0, 0, false });
        for (size_t ti = 0; ti < d_types.size(); ti++) {
            const int type = d_types[ti];
            cursor& c = cursors[ti];
            if (!c.valid) {
                // The first slot whose window starts in the buffer
                tr->next_slot(begin_ns + 1000000 - d_offset_ns, type, c.frame, c.slot);
                c.valid = true;
            }
            for (;;) {
                const int64_t t = tr->slot_time_ns(c.frame, c.slot, type) + d_offset_ns;
                int64_t first, end;
                slot_window(type, sample_index(t), &first, &end);
                if (end > (int64_t)d_end || first > (int64_t)limit) {
                    break;
                }
                if (first < (int64_t)begin) {
                    d_n_missed_slots++;
                } else {
                    const double doppler = tr->doppler(t);
                    if (already_probed(type, t, doppler)) {
                        d_n_shared_slots++;
                    } else {
                        d_probed.push_back({ type, t, doppler });
                        probe_slot(*tr, type, c.frame, c.slot, t);
                    }
                }
                if (++c.slot == tch_c_timing::N_SLOTS) {
                    c.slot = 0;
                    c.frame++;
                }
            }
        }

        if (d_probe_uplink && d_ul_w_valid) {
            probe_ul_slots(*tr, ul_limit, begin, begin_ns);
        } else if (d_probe_uplink) {
            probe_ul_sweep(*tr, ul_limit, begin, begin_ns);
        }
    }
}

/*
 * UL slots of the tracks, once W0 is known: per track a cursor over (frame, slot) with
 * the nominal times v = W0 + UL_SLOT_MS * slot. If W0 moves by more than a ms (a new
 * estimate), the cursor restarts after the last slot probed.
 */
void tch_c_slot_prober_impl::probe_ul_slots(const tch_c_sat_track& tr,
                                            uint64_t ul_limit,
                                            uint64_t begin,
                                            int64_t begin_ns)
{
    auto ins = d_ul_cursors.insert({ tr.id(), { { 0, 0, false }, begin_ns, d_ul_w0 } });
    ul_cursor& uc = ins.first->second;
    cursor& c = uc.c;
    if (c.valid && std::fabs(tch_c_timing::wrap_ms(uc.w0 - d_ul_w0)) > 1) {
        c.valid = false;
    }
    if (!c.valid) {
        const int64_t t0 = std::max(begin_ns + 2000000, uc.last_ns + 2000000);
        c.frame = tr.ul_frame_at(t0, d_ul_w0) - 1;
        c.slot = 0;
        while (tr.ul_time_ns(c.frame, d_ul_w0 + tch_c_timing::UL_SLOT_MS * c.slot) < t0) {
            if (++c.slot == tch_c_timing::N_SLOTS) {
                c.slot = 0;
                c.frame++;
            }
        }
        c.valid = true;
    }
    uc.w0 = d_ul_w0;
    for (;;) {
        const int64_t t =
            tr.ul_time_ns(c.frame, d_ul_w0 + tch_c_timing::UL_SLOT_MS * c.slot);
        int64_t first, end;
        ul_slot_window(sample_index(t), UL_WIDE_S, &first, &end);
        if (end > (int64_t)d_end || first > (int64_t)ul_limit) {
            break;
        }
        // Two tracks of one satellite: its UL slots are probed once (type -1)
        const double doppler = tr.doppler(t);
        if (first >= (int64_t)begin && !already_probed(-1, t, doppler)) {
            d_probed.push_back({ -1, t, doppler });
            probe_ul_slot(tr, c.frame, c.slot, t);
        }
        uc.last_ns = t;
        if (++c.slot == tch_c_timing::N_SLOTS) {
            c.slot = 0;
            c.frame++;
        }
    }
}

/*
 * The samples a UL probe around p (search +-s seconds) reads: the search window (and
 * the margin of fine_timing()) and the extraction of a burst at its end.
 */
void tch_c_slot_prober_impl::ul_slot_window(double p,
                                            double s,
                                            int64_t* first,
                                            int64_t* end) const
{
    const double S = s * d_sample_rate;
    *first = (int64_t)std::floor(p - S) - d_survey.uw_offset(1) -
             d_survey.timing_margin(1) - 1;
    long ef, ee;
    d_survey.extract_span(1, p + S + d_survey.timing_margin(1), &ef, &ee);
    *end = std::max<int64_t>(*first + d_survey.window_capacity(s > UL_WIDE_S), ee);
}

// Is the (track, slot, position) link active at t: a UL frame within UL_TRACK_NS
// (tracked) or a DL frame within UL_DL_ACTIVE_NS
static bool ul_tracked(int64_t t_ns, int64_t last_ul_ns)
{
    return std::llabs(t_ns - last_ul_ns) < UL_TRACK_NS;
}

static bool dl_active(int64_t t_ns, int64_t last_dl_ns)
{
    return std::llabs(t_ns - last_dl_ns) < UL_DL_ACTIVE_NS;
}

/*
 * Search one channel for a UL burst near p_c (+-s seconds, sample index of the stream,
 * x holds the window from `first`), carrier f_c +- f_s, and publish it.
 */
uint64_t tch_c_slot_prober_impl::ul_search(const gr_complex* x,
                                           int64_t first,
                                           const tch_c_sat_track& track,
                                           int slot,
                                           int position,
                                           double p_c,
                                           double s,
                                           double f_c,
                                           double f_s,
                                           tch_c_slot_survey::channel_candidate* cc,
                                           double* p_found)
{
    d_n_ul_searches++;
    const float threshold = std::min(UL_THRESHOLD, d_debug_rho2);
    const double S = s * d_sample_rate;
    if (!d_survey.search_channel(
            1, true, f_c, f_s, p_c - S - first, p_c + S - first, threshold, cc)) {
        return UINT64_MAX;
    }
    double p_uw = first + cc->p_uw;
    uint64_t id = 0;
    if (cc->rho2 >= UL_THRESHOLD) {
        const int D = d_survey.uw_offset(1);
        tch_c_slot_survey::candidate c;
        c.position = position;
        c.nominal_frequency = tch_c_slot_survey::nominal_frequency(1, position);
        c.frequency = f_c;
        c.estimated_frequency = cc->frequency;
        c.rho2 = cc->rho2;
        c.energy = 0;
        c.magnitude = cc->magnitude;
        c.bin = 0;
        const int64_t i0 = std::lround(p_uw) - D;
        p_uw = i0 + D + d_survey.fine_timing(x + (i0 - first), 1, true, c, nullptr);
        const int n =
            d_survey.extract(x, p_uw - first, 1, cc->frequency, true, &d_symbols[0]);
        id = publish_probe(1,
                           time_ns(p_uw),
                           cc->frequency,
                           &d_symbols[0],
                           n,
                           cc->magnitude,
                           d_survey.window_noise(),
                           cc->rho2,
                           track.id(),
                           position,
                           slot,
                           true);
    }
    *p_found = p_uw;
    return id;
}

/*
 * One UL slot of a track: search the channels of its links (see ul_link) and publish
 * what reaches UL_THRESHOLD after fine timing and extraction.
 */
void tch_c_slot_prober_impl::probe_ul_slot(const tch_c_sat_track& track,
                                           int64_t frame,
                                           int slot,
                                           int64_t t_ns)
{
    d_n_ul_slots++;
    auto first_link = d_ul_links.lower_bound(std::make_tuple(track.id(), slot, 0));
    auto in_slot = [&](decltype(first_link) it) {
        return it != d_ul_links.end() && std::get<0>(it->first) == track.id() &&
               std::get<1>(it->first) == slot;
    };
    bool any = false;
    for (auto it = first_link; in_slot(it); ++it) {
        if (ul_tracked(t_ns, it->second.last_ul_ns) ||
            dl_active(t_ns, it->second.last_dl_ns)) {
            any = true;
            break;
        }
    }
    if (!any) {
        return;
    }

    const double t0 = d_debug ? cpu_s() : 0;
    const double p_nom = sample_index(t_ns);
    int64_t first, end;
    ul_slot_window(p_nom, UL_WIDE_S, &first, &end);
    const gr_complex* x = samples(first, end - first);
    d_survey.set_window(x, d_survey.window_capacity());
    const double doppler = track.doppler(t_ns);
    const double v_nom = d_ul_w0 + tch_c_timing::UL_SLOT_MS * slot;

    for (auto it = first_link; in_slot(it); ++it) {
        const ul_link& l = it->second;
        const int position = std::get<2>(it->first);
        const bool tracked = ul_tracked(t_ns, l.last_ul_ns);
        if (!tracked && !dl_active(t_ns, l.last_dl_ns)) {
            continue;
        }
        const double f0 = tch_c_slot_survey::nominal_frequency(1, position) *
                          (1 - doppler / TCH_C_F_RING);
        double p_c, s, f_c, f_s;
        if (tracked) {
            // Where the last UL frame of the link was, within the wide search
            const double dv = std::max(-(UL_WIDE_S - UL_TRACK_S) * 1e3,
                                       std::min((UL_WIDE_S - UL_TRACK_S) * 1e3,
                                                tch_c_timing::wrap_ms(l.v_ms - v_nom)));
            p_c = p_nom + dv * 1e-3 * d_sample_rate;
            s = UL_TRACK_S;
            f_c = f0 + l.f_offset_hz;
            f_s = UL_TRACK_SEARCH_HZ;
        } else {
            p_c = p_nom;
            s = UL_WIDE_S;
            f_c = f0 + UL_WIDE_F_HZ;
            f_s = UL_WIDE_SEARCH_HZ;
        }
        tch_c_slot_survey::channel_candidate cc;
        double p_uw;
        const uint64_t id =
            ul_search(x, first, track, slot, position, p_c, s, f_c, f_s, &cc, &p_uw);
        if (id != UINT64_MAX && cc.rho2 >= d_debug_rho2) {
            fprintf(stderr,
                    "PU t_ns=%" PRId64 " track=%d frame=%" PRId64
                    " slot=%d pos=%d tracked=%d f=%.0f df=%.0f rho2=%.3f mag=%.1f "
                    "dt_us=%.2f id=%" PRIu64 "\n",
                    t_ns,
                    track.id(),
                    frame,
                    slot,
                    position,
                    tracked,
                    cc.frequency,
                    cc.frequency - f0,
                    cc.rho2,
                    cc.magnitude,
                    (p_uw - p_nom) / d_sample_rate * 1e6,
                    id);
        }
    }
    if (d_debug) {
        d_cpu[3] += cpu_s() - t0;
    }
}

/*
 * Acquisition while W0 is unknown: per track and frame a window of v (2 UL_SWEEP_S
 * wide) searched on all channels with a DL frame recently, the next frame the next
 * window (UL_SWEEP_WINDOWS per frame): all of the 90 ms frame in ~12 frames. After
 * UL_SWEEP_FRAMES windows without a
 * UL frame only every UL_SWEEP_SLOW-th frame is searched.
 */
void tch_c_slot_prober_impl::probe_ul_sweep(const tch_c_sat_track& tr,
                                            uint64_t ul_limit,
                                            uint64_t begin,
                                            int64_t begin_ns)
{
    // Half width of the search: what the long window holds
    const double sps = d_survey.samples_per_symbol(1);
    const double S =
        std::min(UL_SWEEP_S * d_sample_rate,
                 0.5 * (d_survey.window_capacity(true) - d_survey.uw_offset(1) -
                        d_survey.timing_margin(1) - (TCH_C_UW_LENGTH + 6) * sps) -
                     2);
    const double span_ms = 2 * S / d_sample_rate * 1e3;

    auto ins = d_ul_sweeps.insert({ tr.id(), { 0, 0, 0, 0, false } });
    ul_sweep& sw = ins.first->second;
    if (!sw.valid) {
        sw.frame = tr.frame_at(begin_ns + 10000000);
        sw.valid = true;
    }
    for (;;) {
        const double v_c = std::fmod(sw.v + span_ms / 2, tch_c_timing::FRAME_MS);
        const int64_t t = tr.ul_time_ns(sw.frame, v_c);
        const double p_c = sample_index(t);
        int64_t first, end;
        ul_slot_window(p_c, S / d_sample_rate, &first, &end);
        if (end > (int64_t)d_end || first > (int64_t)ul_limit) {
            break;
        }
        // UL_SWEEP_WINDOWS windows per frame
        if (++sw.window == UL_SWEEP_WINDOWS) {
            sw.window = 0;
            sw.frame++;
        }
        const double doppler = tr.doppler(t);
        if (first < (int64_t)begin || already_probed(-2, t, doppler)) {
            continue;
        }
        d_probed.push_back({ -2, t, doppler });
        d_n_ul_sweep_frames++;
        sw.n++;
        if (sw.n > UL_SWEEP_FRAMES && sw.n % UL_SWEEP_SLOW != 0) {
            continue;
        }
        sw.v = std::fmod(sw.v + span_ms, tch_c_timing::FRAME_MS);

        // The channels with DL frames, each once (the slot of its newest DL frame)
        std::map<int, std::pair<int, int64_t>> channels;
        for (auto it = d_ul_links.lower_bound(std::make_tuple(tr.id(), 0, 0));
             it != d_ul_links.end() && std::get<0>(it->first) == tr.id();
             ++it) {
            if (!dl_active(t, it->second.last_dl_ns)) {
                continue;
            }
            const int pos = std::get<2>(it->first);
            auto c = channels.find(pos);
            if (c == channels.end() || it->second.last_dl_ns > c->second.second) {
                channels[pos] = { std::get<1>(it->first), it->second.last_dl_ns };
            }
        }
        if (channels.empty()) {
            continue;
        }
        const double t0 = d_debug ? cpu_s() : 0;
        const gr_complex* x = samples(first, end - first);
        d_survey.set_window(x, d_survey.window_capacity(true));
        for (const auto& ch : channels) {
            const int position = ch.first;
            const double f0 = tch_c_slot_survey::nominal_frequency(1, position) *
                              (1 - doppler / TCH_C_F_RING);
            tch_c_slot_survey::channel_candidate cc;
            double p_uw;
            const uint64_t id = ul_search(x,
                                          first,
                                          tr,
                                          ch.second.first,
                                          position,
                                          p_c,
                                          S / d_sample_rate,
                                          f0 + UL_WIDE_F_HZ,
                                          UL_WIDE_SEARCH_HZ,
                                          &cc,
                                          &p_uw);
            if (id != UINT64_MAX && cc.rho2 >= d_debug_rho2) {
                fprintf(stderr,
                        "PW t_ns=%" PRId64 " track=%d slot=%d pos=%d f=%.0f df=%.0f "
                        "rho2=%.3f mag=%.1f v=%.3f id=%" PRIu64 "\n",
                        t,
                        tr.id(),
                        ch.second.first,
                        position,
                        cc.frequency,
                        cc.frequency - f0,
                        cc.rho2,
                        cc.magnitude,
                        tr.ul_phase_ms(time_ns(p_uw)),
                        id);
            }
        }
        if (d_debug) {
            d_cpu[3] += cpu_s() - t0;
        }
    }
}

/*
 * The track and TCH_C1 position of a UL burst at t on f: a track with a DL link on the
 * channel whose UL carrier (nominal * (1 - Doppler / F_RING) + 0..5.5 kHz) is within
 * tol_hz. *slot: the slot of the link whose v fits W0, if known, else the DL slot with
 * the newest frame on the channel.
 */
bool tch_c_slot_prober_impl::ul_match(
    int64_t t_ns, double f_hz, double tol_hz, int* track, int* position, int* slot) const
{
    for (const auto& tr : d_tracks) {
        if (tr.empty()) {
            continue;
        }
        const double d = tr.doppler(t_ns);
        const int pos =
            nearest_position(1, (f_hz - UL_WIDE_F_HZ) / (1 - d / TCH_C_F_RING));
        const double f0 =
            tch_c_slot_survey::nominal_frequency(1, pos) * (1 - d / TCH_C_F_RING);
        if (std::fabs(f_hz - f0 - UL_WIDE_F_HZ) > tol_hz) {
            continue;
        }
        int best = -1;
        int64_t best_t = UL_NEVER;
        for (int k = 0; k < tch_c_timing::N_SLOTS; k++) {
            auto it = d_ul_links.find(std::make_tuple(tr.id(), k, pos));
            if (it == d_ul_links.end()) {
                continue;
            }
            const int64_t last = std::max(it->second.last_dl_ns, it->second.last_ul_ns);
            if ((dl_active(t_ns, it->second.last_dl_ns) ||
                 ul_tracked(t_ns, it->second.last_ul_ns)) &&
                (best < 0 || last > best_t)) {
                best = k;
                best_t = last;
            }
        }
        if (best < 0) {
            continue;
        }
        if (d_ul_w_valid) {
            const double x = tch_c_timing::wrap_ms(tr.ul_phase_ms(t_ns) - d_ul_w0);
            const int k = (int)std::lround(x / tch_c_timing::UL_SLOT_MS);
            if (k >= 0 && k < tch_c_timing::N_SLOTS &&
                std::fabs(x - k * tch_c_timing::UL_SLOT_MS) < 1) {
                best = k;
            }
        }
        *track = tr.id();
        *position = pos;
        *slot = best;
        return true;
    }
    return false;
}

/*
 * A CRC ok TCH_C1 uplink frame: where the UL burst of its link is, and a measurement of
 * W0. Our own probes know their track, slot and channel; frames of the energy path are
 * matched to a channel with DL frames (ul_match()).
 */
void tch_c_slot_prober_impl::on_ul_hit(const tch_c_hit& hit)
{
    if (!d_probe_uplink || hit.type != 1 || d_offset_ns != 0) {
        return;
    }
    int track, slot, position;
    if (hit.id >= TCH_C_PROBE_ID_BASE) {
        if (!hit.probe || hit.probe->track < 0 || !hit.probe->uplink) {
            return;
        }
        track = hit.probe->track;
        slot = hit.probe->slot;
        position = hit.probe->position;
    } else if (!ul_match(hit.t_ns, hit.f_hz, UL_MATCH_TOL_HZ, &track, &position, &slot)) {
        return;
    }
    const tch_c_sat_track* tr = d_tracks.find(track);
    if (!tr) {
        return;
    }
    const double v = tr->ul_phase_ms(hit.t_ns);
    const double nominal = tch_c_slot_survey::nominal_frequency(1, position);
    auto it = d_ul_links
                  .insert({ std::make_tuple(track, slot, position),
                            { UL_NEVER, UL_NEVER, 0, 0 } })
                  .first;
    ul_link& l = it->second;
    if ((int64_t)hit.t_ns >= l.last_ul_ns) {
        l.last_ul_ns = hit.t_ns;
        l.v_ms = v;
        l.f_offset_hz = hit.f_hz - nominal * (1 - tr->doppler(hit.t_ns) / TCH_C_F_RING);
    }

    // W0: median of the recent measurements, unwrapped around the newest one
    const double w =
        std::fmod(v - tch_c_timing::UL_SLOT_MS * slot + 2 * tch_c_timing::FRAME_MS,
                  tch_c_timing::FRAME_MS);
    d_ul_w.push_back({ (int64_t)hit.t_ns, w });
    while (d_ul_w.size() > UL_W_MEASUREMENTS) {
        d_ul_w.pop_front();
    }
    std::vector<std::pair<double, double>> m;
    for (const auto& p : d_ul_w) {
        if (p.first > (int64_t)hit.t_ns - UL_W_WINDOW_NS) {
            m.push_back({ w + tch_c_timing::wrap_ms(p.second - w), 1.0 });
        }
    }
    const double w0 = tch_c_timing::weighted_median(m);
    const double w0_old = d_ul_w0;
    d_ul_w0 = std::fmod(w0 + tch_c_timing::FRAME_MS, tch_c_timing::FRAME_MS);
    if (d_debug &&
        (!d_ul_w_valid || std::fabs(tch_c_timing::wrap_ms(d_ul_w0 - w0_old)) > 0.05)) {
        fprintf(stderr,
                "UW0 t_ns=%" PRIu64 " w0=%.3f (was %.3f, valid %d) from track=%d slot=%d "
                "pos=%d v=%.3f\n",
                hit.t_ns,
                d_ul_w0,
                w0_old,
                d_ul_w_valid,
                track,
                slot,
                position,
                v);
    }
    d_ul_w_valid = true;
}

/*
 * One DL slot of a track and type: survey all raster positions at the predicted time
 * and Doppler, and publish those above the threshold after fine timing and extraction.
 */
void tch_c_slot_prober_impl::probe_slot(
    const tch_c_sat_track& track, int type, int64_t frame, int slot, int64_t t_ns)
{
    const double p = sample_index(t_ns);
    int64_t first, end;
    slot_window(type, p, &first, &end);
    const gr_complex* x = samples(first, end - first);
    const int D = d_survey.uw_offset(type);
    const int64_t i0 = std::lround(p) - D;
    const gr_complex* seg = x + (i0 - first);

    const double doppler = track.doppler(t_ns);
    const bool stale = track.stale(t_ns);
    const float threshold = d_threshold[tch_c_timing::type_index(type)];
    d_candidates.clear();
    double t0 = d_debug ? cpu_s() : 0;
    d_survey.survey(seg,
                    type,
                    false,
                    doppler,
                    stale ? STALE_SEARCH_HZ : 0,
                    std::min(threshold, d_debug_rho2),
                    d_candidates);
    d_n_surveys[tch_c_timing::type_index(type)]++;
    if (d_debug) {
        const double t1 = cpu_s();
        d_cpu[0] += t1 - t0;
        t0 = t1;
    }

    for (const auto& c : d_candidates) {
        uint64_t id = 0;
        double tau = 0;
        if (c.rho2 >= threshold) {
            float rho2_fine;
            tau = d_survey.fine_timing(seg, type, false, c, &rho2_fine);
            const double p_uw = i0 + D + tau;
            const int n = d_survey.extract(
                x, p_uw - first, type, c.estimated_frequency, false, &d_symbols[0]);
            id = publish_probe(type,
                               time_ns(p_uw),
                               c.estimated_frequency,
                               &d_symbols[0],
                               n,
                               c.magnitude,
                               d_survey.noise(),
                               c.rho2,
                               track.id(),
                               c.position,
                               slot);
            if (d_debug) {
                const double t1 = cpu_s();
                d_cpu[1] += t1 - t0;
                t0 = t1;
            }
        }
        if (c.rho2 >= d_debug_rho2) {
            fprintf(stderr,
                    "PS t_ns=%" PRId64 " type=%d track=%d n_ira=%" PRIu64
                    " frame=%" PRId64 " slot=%d pos=%d f=%.0f rho2=%.3f mag=%.1f "
                    "stale=%d tau_us=%.2f id=%" PRIu64 "\n",
                    t_ns,
                    type,
                    track.id(),
                    track.n_ira(),
                    frame,
                    slot,
                    c.position,
                    c.estimated_frequency,
                    c.rho2,
                    c.magnitude,
                    stale,
                    tau / d_sample_rate * 1e6,
                    id);
        }
    }
}

/*
 * Wait until the input or a message arrives, at most ms milliseconds: what the
 * scheduler does for a block blocked on its input. work() returning 0 does not make it
 * wait (a block without outputs never produces anything, so the scheduler calls it
 * again at once).
 */
void tch_c_slot_prober_impl::wait_for_change(int ms)
{
    gr::tpb_detail& tpb = detail()->d_tpb;
    gr::thread::scoped_lock guard(tpb.mutex);
    if (!tpb.input_changed) {
        tpb.input_cond.timed_wait(guard, boost::posix_time::milliseconds(ms));
    }
}

int tch_c_slot_prober_impl::work(int noutput_items,
                                 gr_vector_const_void_star& input_items,
                                 gr_vector_void_star& output_items)
{
    // Reading from a file: wait for tch_c_demod. The acknowledgements arrive as
    // messages, handled between the calls of work().
    if (!d_drop_overflow && blocked()) {
        wait_for_change(10);
        return 0;
    }

    /*
     * The last input sample is never consumed: GNU Radio does not call work() for a
     * finished input without samples, but it does for this one. That call probes the
     * rest of the buffer and waits for tch_c_demod to handle the probes, then ends the
     * block.
     */
    const bool input_done = detail()->input(0)->done();
    const int n = noutput_items - 1;
    const gr_complex* in = (const gr_complex*)input_items[0];
    assert(nitems_read(0) == d_end);

    const double t_work = d_debug ? cpu_s() : 0;
    int consumed = 0;
    while (consumed < n) {
        const int chunk = std::min<int>(n - consumed, d_chunk);
        update_time_base(d_end, d_end + chunk);
        append(in + consumed, chunk);
        consumed += chunk;

        while (!d_iras.empty()) {
            on_ira(d_iras.front().first, d_iras.front().second);
            d_iras.pop_front();
        }
        if (d_end > d_lag) {
            const uint64_t ul_lag = d_lag + UL_DELAY_S * d_sample_rate;
            probe_slots(d_end - d_lag, d_end > ul_lag ? d_end - ul_lag : 0);
        }
        if (!d_drop_overflow && blocked()) {
            break;
        }
    }
    forget_old_probes();
    if (d_debug) {
        d_cpu[2] += cpu_s() - t_work;
    }

    if (input_done && consumed == n) {
        if (!d_flushed) {
            while (!d_iras.empty()) {
                on_ira(d_iras.front().first, d_iras.front().second);
                d_iras.pop_front();
            }
            probe_slots(d_end, d_end);
            d_flushed = true;
        }
        if (d_outstanding == 0 || d_drop_overflow) {
            return WORK_DONE;
        }
    }
    if (consumed == 0) {
        wait_for_change(input_done ? 10 : 100);
    }
    return consumed;
}

} /* namespace iridium */
} /* namespace gr */

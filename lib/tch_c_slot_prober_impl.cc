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

// tch_c_demod tries all its carrier recovery variants on probes with at least this
// rho2 (a UW for sure), only the fit to the known symbols on the others
static const float VARIANTS_RHO2 = 0.6f;

// Strong hits (rho2 of the survey) update the timing bias of their type
static const float BIAS_MIN_RHO2 = 0.6f;

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
      d_threshold(threshold),
      d_outstanding_limit(outstanding_limit),
      d_drop_overflow(drop_overflow),
      d_uplink(uplink),
      d_debug(env_flag("IRIDIUM_PROBE_DEBUG")),
      d_debug_rho2(env_double("IRIDIUM_PROBE_DEBUG_RHO2", 2)),
      d_offset_ns(std::llround(env_double("IRIDIUM_PROBE_OFFSET_US", 0) * 1e3)),
      d_variants_rho2(env_double("IRIDIUM_PROBE_VARIANTS_RHO2", VARIANTS_RHO2)),
      d_ring(std::max<size_t>(std::ceil((history_s + RING_MARGIN_S) * sample_rate), 1)),
      d_end(0),
      d_chunk(std::max<uint64_t>(1, CHUNK_S * sample_rate)),
      d_lag(std::max<double>(0, history_s) * sample_rate),
      d_flushed(false),
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
                "and extraction %.2f s\n",
                d_cpu[2],
                d_cpu[0],
                d_cpu[1]);
        for (int i = 0; i < 3; i++) {
            fprintf(stderr,
                    "tch_c_slot_prober: TCH_C%d %" PRIu64 " slot surveys, %" PRIu64
                    " probes, %" PRIu64 " CRC ok\n",
                    i == 2 ? 8 : i + 1,
                    d_n_surveys[i],
                    d_n_type_probes[i],
                    d_n_type_crc_ok[i]);
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
                                               int slot)
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

    d_probes[id] = { type, t_first_uw_ns, freq_hz, rho2, track, position, slot, false };
    d_unacked_times.insert(t_first_uw_ns);
    d_outstanding++;
    d_newest = std::max(d_newest, t_first_uw_ns);
    d_n_probes++;
    d_n_type_probes[tch_c_timing::type_index(type)]++;

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
            d_n_type_crc_ok[tch_c_timing::type_index(hit.type)]++;
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
    if (hit.uplink || !hit.crc_ok) {
        return;
    }
    if (hit.id >= TCH_C_PROBE_ID_BASE) {
        // In the control run (shifted probe times) the probes are not where the
        // bursts are
        if (!hit.probe || hit.probe->track < 0 || d_offset_ns != 0) {
            return;
        }
        const tch_c_probe_info& p = *hit.probe;
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
        d_tracks.add_hit(
            best, hit.t_ns, hit.f_hz, hit.type, best_match.slot, best_nominal, false);
        d_n_energy_associated++;
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
void tch_c_slot_prober_impl::probe_slots(uint64_t limit)
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
    }
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
    d_candidates.clear();
    double t0 = d_debug ? cpu_s() : 0;
    d_survey.survey(seg,
                    type,
                    false,
                    doppler,
                    stale ? STALE_SEARCH_HZ : 0,
                    std::min(d_threshold, d_debug_rho2),
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
        if (c.rho2 >= d_threshold) {
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
            probe_slots(d_end - d_lag);
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
            probe_slots(d_end);
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

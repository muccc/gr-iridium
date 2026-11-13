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
#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <thread>

namespace gr {
namespace iridium {

// The same limit as burst_ack_counter::MAX_BACKLOG: the frame sorter needs the frames
// within 2 s
static const uint64_t PROBE_MAX_BACKLOG = 1000000000; // ns

// Acknowledged probes are forgotten this long (of the time of the bursts) after the
// newest one: tch_c_demod can report the frame after the acknowledgement
static const uint64_t PROBE_KEEP = 3000000000; // ns

// The ring buffer holds this much more than the history: one extraction window
static const double RING_MARGIN_S = 0.02;

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

static bool env_flag(const char* name)
{
    const char* v = std::getenv(name);
    return v && *v && strcmp(v, "0") != 0;
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
      d_ring(std::max<size_t>(std::ceil((history_s + RING_MARGIN_S) * sample_rate), 1)),
      d_end(0),
      d_time_base_offset(0),
      d_time_base_timestamp(0),
      d_next_id(TCH_C_PROBE_ID_BASE),
      d_outstanding(0),
      d_newest(0),
      d_n_probes(0),
      d_n_dropped(0),
      d_n_iras(0),
      d_n_hits(0),
      d_n_probe_hits(0),
      d_n_probe_crc_ok(0),
      d_selftest_next(0),
      d_selftest_missed(0)
{
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

    const char* selftest = std::getenv("IRIDIUM_PROBE_SELFTEST");
    if (selftest && *selftest) {
        selftest_load(selftest);
    }
}

tch_c_slot_prober_impl::~tch_c_slot_prober_impl() {}

bool tch_c_slot_prober_impl::stop()
{
    if (d_debug || !d_selftest.empty()) {
        fprintf(stderr,
                "tch_c_slot_prober: %" PRIu64 " ring alerts, %" PRIu64 " probes (%" PRIu64
                " dropped), %" PRIu64 " probe frames (%" PRIu64 " CRC ok), %" PRIu64
                " frames from the energy path, %d outstanding\n",
                d_n_iras,
                d_n_probes,
                d_n_dropped,
                d_n_probe_hits,
                d_n_probe_crc_ok,
                d_n_hits - d_n_probe_hits,
                d_outstanding);
        if (!d_selftest.empty()) {
            fprintf(
                stderr,
                "tch_c_slot_prober: self-test: %zu known frames, %zu reached, %" PRIu64
                " no longer in the buffer\n",
                d_selftest.size(),
                d_selftest_next,
                d_selftest_missed);
        }
    }
    return true;
}

/*
 * Time base
 */
double tch_c_slot_prober_impl::sample_index(uint64_t t_ns) const
{
    return d_time_base_offset +
           (double)(int64_t)(t_ns - d_time_base_timestamp) * d_sample_rate / 1e9;
}

uint64_t tch_c_slot_prober_impl::time_ns(double sample_index) const
{
    return d_time_base_timestamp +
           (int64_t)std::floor((sample_index - d_time_base_offset) * 1e9 / d_sample_rate);
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

    d_probes[id] = { type, t_first_uw_ns, freq_hz, rho2, track, position, slot, false };
    d_unacked_times.insert(t_first_uw_ns);
    d_outstanding++;
    d_newest = std::max(d_newest, t_first_uw_ns);
    d_n_probes++;

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
    if (d_debug) {
        fprintf(stderr,
                "tch_c_slot_prober: IRA t_ns=%" PRIu64 " sample=%.1f f=%.1f\n",
                t,
                sample_index(t),
                f);
    }
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
        }
    }

    on_hit(hit);
}

/*
 * HOOK: a ring alert at t_ns (time of its first UW symbol) on f_hz. The track table
 * takes it from here.
 */
void tch_c_slot_prober_impl::on_ira(uint64_t t_ns, double f_hz) {}

/*
 * HOOK: a TCH_C frame. CRC ok downlink frames update the track table.
 */
void tch_c_slot_prober_impl::on_hit(const tch_c_hit& hit) {}

int tch_c_slot_prober_impl::work(int noutput_items,
                                 gr_vector_const_void_star& input_items,
                                 gr_vector_void_star& output_items)
{
    // Reading from a file: wait for tch_c_demod. The acknowledgements arrive as
    // messages, handled between the calls of work().
    if (!d_drop_overflow && blocked()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        return 0;
    }

    std::vector<tag_t> tags;
    get_tags_in_window(tags, 0, 0, noutput_items, pmt::mp("time_base"));
    if (!tags.empty()) {
        std::sort(tags.begin(), tags.end(), tag_t::offset_compare);
        d_time_base_offset = pmt::to_uint64(pmt::tuple_ref(tags.back().value, 0));
        d_time_base_timestamp = pmt::to_uint64(pmt::tuple_ref(tags.back().value, 1));
    }

    // Append the input to the ring buffer
    const gr_complex* in = (const gr_complex*)input_items[0];
    assert(nitems_read(0) == d_end);
    const size_t size = d_ring.size();
    size_t n = noutput_items;
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

    while (!d_iras.empty()) {
        on_ira(d_iras.front().first, d_iras.front().second);
        d_iras.pop_front();
    }

    if (!d_selftest.empty()) {
        selftest_run();
    }

    forget_old_probes();

    return noutput_items;
}

/*
 * TEMPORARY self-test of the probe path, remove with the real probing.
 *
 * IRIDIUM_PROBE_SELFTEST names a file of known TCH_C frames, one per line:
 * "t_ns f_hz type" with the time of the first UW symbol (as in the timestamps of the
 * frames) and the carrier. Once its samples are in the buffer, every frame of the
 * types of this block is cut out with a plain rotation, RRC matched filter and
 * decimation to one sample per symbol, and published as a probe. With
 * IRIDIUM_PROBE_DEBUG the UW correlation is evaluated around the expected time and the
 * timing error printed ("selftest: ... dt_us=").
 */
void tch_c_slot_prober_impl::selftest_load(const char* filename)
{
    std::ifstream file(filename);
    if (!file) {
        fprintf(stderr, "tch_c_slot_prober: cannot read %s\n", filename);
        return;
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream s(line);
        selftest_frame frame;
        if (s >> frame.t_ns >> frame.f_hz >> frame.type) {
            if (std::find(d_types.begin(), d_types.end(), frame.type) != d_types.end()) {
                d_selftest.push_back(frame);
            }
        }
    }
    std::sort(
        d_selftest.begin(),
        d_selftest.end(),
        [](const selftest_frame& a, const selftest_frame& b) { return a.t_ns < b.t_ns; });
    fprintf(stderr,
            "tch_c_slot_prober: self-test with %zu known frames from %s\n",
            d_selftest.size(),
            filename);
}

// Root raised cosine pulse, t in symbols
static float rrc(double t, double alpha)
{
    if (std::abs(t) < 1e-9) {
        return 1 - alpha + 4 * alpha / M_PI;
    }
    if (std::abs(std::abs(t) - 1 / (4 * alpha)) < 1e-9) {
        return alpha / M_SQRT2 *
               ((1 + 2 / M_PI) * sin(M_PI / (4 * alpha)) +
                (1 - 2 / M_PI) * cos(M_PI / (4 * alpha)));
    }
    return (sin(M_PI * t * (1 - alpha)) + 4 * alpha * t * cos(M_PI * t * (1 + alpha))) /
           (M_PI * t * (1 - (4 * alpha * t) * (4 * alpha * t)));
}

static const double SELFTEST_RRC_ALPHA = 0.2;
static const int SELFTEST_RRC_SPAN = 6; // symbols on each side

void tch_c_slot_prober_impl::selftest_run()
{
    while (d_selftest_next < d_selftest.size()) {
        const selftest_frame& frame = d_selftest[d_selftest_next];
        const tch_c_format& f = tch_c_format_by_type(frame.type);
        const double sps = d_sample_rate / f.symbols_per_second;
        const double n0 = sample_index(frame.t_ns);
        // Two symbols more on each side for the timing search of the debug output
        const double reach = (SELFTEST_RRC_SPAN + 2) * sps + 2;
        const double first = n0 - reach;
        const double last = n0 + f.n_symbols * sps + reach;
        if (last >= d_end) {
            return;
        }
        d_selftest_next++;
        if (first < 0 || (uint64_t)first < ring_begin()) {
            d_selftest_missed++;
            continue;
        }
        selftest_extract(frame);
    }
}

void tch_c_slot_prober_impl::selftest_extract(const selftest_frame& frame)
{
    const tch_c_format& f = tch_c_format_by_type(frame.type);
    const double sps = d_sample_rate / f.symbols_per_second;
    const double n0 = sample_index(frame.t_ns);
    const double reach = (SELFTEST_RRC_SPAN + 2) * sps + 2;
    const uint64_t start = (uint64_t)(n0 - reach);
    const size_t len = (size_t)(n0 + f.n_symbols * sps + reach) - start;

    // Rotate the carrier to 0 Hz
    const gr_complex* x = samples(start, len);
    std::vector<gr_complex> rotated(len);
    const double w = -2 * M_PI * (frame.f_hz - d_center_frequency) / d_sample_rate;
    std::complex<double> phasor(1, 0);
    const std::complex<double> step = std::polar(1.0, w);
    for (size_t i = 0; i < len; i++) {
        rotated[i] = x[i] * gr_complex(phasor);
        phasor *= step;
        if ((i & 1023) == 1023) {
            phasor /= std::abs(phasor);
        }
    }

    // Matched filter output at a (fractional) sample position
    const int half = std::ceil(SELFTEST_RRC_SPAN * sps);
    std::vector<float> taps(2 * half + 1);
    double taps_frac = -1;
    auto filter_at = [&](double position) {
        const double rel = position - start;
        const int base = std::floor(rel);
        const double frac = rel - base;
        if (std::abs(frac - taps_frac) > 1e-6) {
            for (int j = -half; j <= half; j++) {
                taps[j + half] = rrc((j - frac) / sps, SELFTEST_RRC_ALPHA);
            }
            taps_frac = frac;
        }
        gr_complex y;
        volk_32fc_32f_dot_prod_32fc(&y, &rotated[base - half], &taps[0], taps.size());
        return y;
    };

    // Timing check: UW correlation around the expected time
    float dt_us = 0, rho2 = 0;
    {
        float best = -1, best_dt = 0;
        float energy_at_best = 1;
        const int steps = 20;
        for (int s = -2 * steps; s <= 2 * steps; s++) {
            if (!d_debug && s != 0) {
                continue;
            }
            const double dt = s * sps / steps;
            gr_complex acc_dl = 0, acc_ul = 0;
            float energy = 0;
            for (int k = 0; k < TCH_C_UW_LENGTH; k++) {
                const gr_complex y = filter_at(n0 + dt + k * sps);
                acc_dl += y * (TCH_C_UW_DL[k] ? -1.f : 1.f);
                acc_ul += y * (TCH_C_UW_UL[k] ? -1.f : 1.f);
                energy += std::norm(y);
            }
            const float c = std::max(std::norm(acc_dl), std::norm(acc_ul));
            if (c > best) {
                best = c;
                best_dt = dt;
                energy_at_best = energy;
            }
            if (s == 0) {
                rho2 = c / (TCH_C_UW_LENGTH * energy);
            }
        }
        dt_us = best_dt / d_sample_rate * 1e6;
        if (d_debug) {
            fprintf(stderr,
                    "tch_c_slot_prober: selftest t_ns=%" PRIu64
                    " f=%.0f type=%d rho2=%.3f dt_us=%.2f rho2_best=%.3f\n",
                    frame.t_ns,
                    frame.f_hz,
                    frame.type,
                    rho2,
                    dt_us,
                    best / (TCH_C_UW_LENGTH * energy_at_best));
        }
    }

    std::vector<gr_complex> symbols(f.n_symbols + 1);
    for (int k = 0; k <= f.n_symbols; k++) {
        symbols[k] = filter_at(n0 + k * sps);
    }
    publish_probe(
        frame.type, frame.t_ns, frame.f_hz, &symbols[0], symbols.size(), 0, 0, rho2);
}

} /* namespace iridium */
} /* namespace gr */

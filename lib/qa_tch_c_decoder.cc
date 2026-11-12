/* -*- c++ -*- */
/*
 * Copyright 2026 gr-iridium author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "qa_tch_c_encoder.h"
#include "tch_c_format.h"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace gr::iridium;

namespace {

const tch_c_format* const FORMATS[] = { &TCH_C1, &TCH_C2, &TCH_C8Q, &TCH_C8A };

const char* name_of(const tch_c_format& f)
{
    if (&f == &TCH_C8A) {
        return "TCH_C8 16APSK";
    }
    return &f == &TCH_C8Q ? "TCH_C8 QPSK" : &f == &TCH_C2 ? "TCH_C2" : "TCH_C1";
}

// Synthetic payload; std::mt19937 output is the same on every platform
std::vector<uint8_t> random_payload(const tch_c_format& f, uint32_t seed)
{
    std::mt19937 rng(seed);
    std::vector<uint8_t> payload(f.n_codewords * f.code.info_bits - 10);
    for (auto& b : payload) {
        b = rng() >> 31;
    }
    return payload;
}

/*
 * Soft bits of a transmitted burst: antipodal (0 -> +1, 1 -> -1) plus Gaussian noise of
 * standard deviation sigma, as LLRs 2 y / sigma^2. sigma = 0: +-4.
 */
std::vector<float> channel(const std::vector<uint8_t>& bits, float sigma, uint32_t seed)
{
    std::mt19937 rng(seed);
    auto uniform = [&]() { return (rng() + 0.5) / 4294967296.0; };
    std::vector<float> llr;
    for (uint8_t b : bits) {
        const float x = b ? -1.f : 1.f;
        if (sigma == 0) {
            llr.push_back(4 * x);
            continue;
        }
        // Box-Muller
        const double n =
            std::sqrt(-2 * std::log(uniform())) * std::cos(2 * M_PI * uniform());
        llr.push_back(2 * (x + sigma * n) / (sigma * sigma));
    }
    return llr;
}

// Like tch_c_demod: decode all codewords, descramble across them, check the CRC-10
bool decode_burst(const tch_c_format& f,
                  const std::vector<float>& llr,
                  std::vector<uint8_t>& msg)
{
    tch_c_decoder decoder(f.code);
    const int K = f.code.info_bits;
    const int L = f.code.coded_bits;
    msg.assign(f.n_codewords * K, 0);
    for (int j = 0; j < f.n_codewords; j++) {
        decoder.decode_raw(&llr[j * L], &msg[j * K]);
    }
    tch_c_decoder::descramble(&msg[0], f.n_codewords * K, 0);
    return tch_c_decoder::check_crc(&msg[0], f.n_codewords * K);
}

bool payload_matches(const std::vector<uint8_t>& payload, const std::vector<uint8_t>& msg)
{
    return std::equal(payload.begin(), payload.end(), msg.begin());
}

std::vector<uint8_t> hex_to_bits(const char* hex, int n_bits)
{
    std::vector<uint8_t> bits;
    for (const char* p = hex; *p; p++) {
        int v = std::stoi(std::string(1, *p), nullptr, 16);
        for (int m = 3; m >= 0; m--) {
            bits.push_back((v >> m) & 1);
        }
    }
    bits.resize(n_bits);
    return bits;
}

} // namespace

BOOST_AUTO_TEST_CASE(test_tch_c_format_tables)
{
    for (const tch_c_format* f : FORMATS) {
        BOOST_TEST_CONTEXT(name_of(*f))
        {
            // Every symbol is UW, pilot or data, exactly once
            std::vector<int> use(f->n_symbols, 0);
            for (int k = 0; k < TCH_C_UW_LENGTH; k++) {
                use[k]++;
            }
            for (int p = 0; p < f->n_pilots; p++) {
                for (int k = 0; k < TCH_C_PILOT_LENGTH; k++) {
                    use.at(f->pilot_start[p] + k)++;
                }
            }
            int n_data = 0;
            for (int s = 0; s < f->n_segments; s++) {
                for (int k = f->segment[s][0]; k < f->segment[s][0] + f->segment[s][1];
                     k++) {
                    use.at(k)++;
                    n_data++;
                }
            }
            for (int k = 0; k < f->n_symbols; k++) {
                BOOST_CHECK_EQUAL(use[k], 1);
            }

            // The data symbols carry the codewords and the padding
            BOOST_CHECK_EQUAL(n_data * f->bits_per_symbol,
                              f->n_codewords * f->code.coded_bits + f->padding_bits);
            BOOST_CHECK(n_data * f->bits_per_symbol <= TCH_C_MAX_DATA_BITS);
            BOOST_CHECK(f->code.coded_bits <= tch_c_decoder::MAX_CODED);
            BOOST_CHECK(f->code.info_bits <= tch_c_decoder::MAX_INFO);
            BOOST_CHECK(f->n_codewords * f->code.info_bits <=
                        tch_c_decoder::MAX_SCRAMBLED);

            // All formats last 8.27 ms
            BOOST_CHECK_CLOSE(f->n_symbols / (double)f->symbols_per_second, 0.00827, 0.1);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_tch_c_decoder_clean)
{
    // Synthetic bursts without noise decode to their payload
    for (const tch_c_format* f : FORMATS) {
        for (uint32_t seed = 1; seed <= 3; seed++) {
            BOOST_TEST_CONTEXT(name_of(*f) << " seed " << seed)
            {
                const std::vector<uint8_t> payload = random_payload(*f, seed);
                const std::vector<uint8_t> tx =
                    qa_tch_c_encoder(*f).encode_burst(payload);
                std::vector<uint8_t> msg;
                BOOST_CHECK(decode_burst(*f, channel(tx, 0, seed), msg));
                BOOST_CHECK(payload_matches(payload, msg));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_tch_c_decoder_noise)
{
    // Gaussian noise decodes: sigma 0.5 for rate 4/5 (TCH_C1), 0.6 for rate 2/3. At these
    // levels 50 of 50 random bursts decode; at 0.6 / 0.7 it is only 18 / 37 to 43 of 50.
    for (const tch_c_format* f : FORMATS) {
        BOOST_TEST_CONTEXT(name_of(*f))
        {
            const float sigma = f == &TCH_C1 ? 0.5f : 0.6f;
            const std::vector<uint8_t> payload = random_payload(*f, 10);
            const std::vector<uint8_t> tx = qa_tch_c_encoder(*f).encode_burst(payload);
            const std::vector<float> llr = channel(tx, sigma, 11);

            // The hard decisions do contain errors
            int errors = 0;
            for (size_t i = 0; i < tx.size(); i++) {
                errors += (llr[i] < 0) != (tx[i] == 1);
            }
            BOOST_CHECK(errors > (int)tx.size() / 50);

            std::vector<uint8_t> msg;
            BOOST_CHECK(decode_burst(*f, llr, msg));
            BOOST_CHECK(payload_matches(payload, msg));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_tch_c_decoder_errors_and_erasures)
{
    // TCH_C1 (rate 4/5): every 40th coded bit wrong but weak, every 40th erased:
    // 10 + 10 of 400 bits. With 4% each it no longer decodes.
    const tch_c_format& f = TCH_C1;
    const std::vector<uint8_t> payload = random_payload(f, 20);
    std::vector<float> llr = channel(qa_tch_c_encoder(f).encode_burst(payload), 0, 0);
    for (size_t i = 0; i < llr.size(); i += 40) {
        llr[i] = -0.5f * llr[i];
    }
    for (size_t i = 20; i < llr.size(); i += 40) {
        llr[i] = 0;
    }

    std::vector<uint8_t> msg;
    BOOST_CHECK(decode_burst(f, llr, msg));
    BOOST_CHECK(payload_matches(payload, msg));
}

BOOST_AUTO_TEST_CASE(test_tch_c_decoder_rejects_noise)
{
    // Far too much noise: the CRC fails
    for (const tch_c_format* f : FORMATS) {
        BOOST_TEST_CONTEXT(name_of(*f))
        {
            const std::vector<uint8_t> payload = random_payload(*f, 30);
            const std::vector<uint8_t> tx = qa_tch_c_encoder(*f).encode_burst(payload);
            std::vector<uint8_t> msg;
            BOOST_CHECK(!decode_burst(*f, channel(tx, 3.f, 31), msg));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_tch_c_crc10)
{
    std::vector<uint8_t> msg = random_payload(TCH_C1, 40);
    const uint16_t crc = tch_c_decoder::crc10(&msg[0], msg.size());
    BOOST_CHECK_EQUAL(crc, qa_tch_c_encoder::crc10(msg));
    for (int i = 9; i >= 0; i--) {
        msg.push_back((crc >> i) & 1);
    }
    BOOST_CHECK(tch_c_decoder::check_crc(&msg[0], msg.size()));

    // Any single bit error is detected
    for (size_t i = 0; i < msg.size(); i++) {
        msg[i] ^= 1;
        BOOST_CHECK(!tch_c_decoder::check_crc(&msg[0], msg.size()));
        msg[i] ^= 1;
    }

    // The register starts at 0 and has no final XOR: all zeros check out
    std::vector<uint8_t> zeros(100, 0);
    BOOST_CHECK(tch_c_decoder::check_crc(&zeros[0], zeros.size()));
}

BOOST_AUTO_TEST_CASE(test_tch_c_descramble)
{
    // DVB energy dispersal PRBS 1 + X^14 + X^15, seed 100101010000000: 03 F6 08 34 ...
    std::vector<uint8_t> bits(32, 0);
    tch_c_decoder::descramble(&bits[0], bits.size(), 0);
    BOOST_CHECK(bits == hex_to_bits("03f60834", 32));

    // The decoder uses the same sequence as the encoder over all codewords
    const int n = tch_c_decoder::MAX_SCRAMBLED;
    std::vector<uint8_t> prbs(n, 0);
    tch_c_decoder::descramble(&prbs[0], n, 0);
    BOOST_CHECK(prbs == qa_tch_c_encoder::dvb_prbs(n));

    // Descrambling twice restores the input, also at an offset
    std::vector<uint8_t> x = hex_to_bits("0123456789abcdef", 64), y = x;
    tch_c_decoder::descramble(&y[0], y.size(), 1616);
    BOOST_CHECK(!std::equal(x.begin(), x.end(), y.begin()));
    tch_c_decoder::descramble(&y[0], y.size(), 1616);
    BOOST_CHECK(std::equal(x.begin(), x.end(), y.begin()));
}

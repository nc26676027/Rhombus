#pragma once

/**
 * SealCompat.h
 *
 * Compatibility layer: maps RBOOT's modified-SEAL (3.6.6-modified) private API calls
 * to standard SEAL 4.1.1 public API.
 *
 * KEY INSIGHT (corrected): seal-modified's "lazy scaling" does NOT mean rescale
 * or multiply_plain skip updating the scale field. In seal-modified, rescale
 * does divide_and_round_q_last AND sets scale /= q_last; multiply_plain sets
 * scale = ct.scale * plain.scale. The "lazy scaling" is implemented in the
 * Bootstrapper's coefficient preprocessing (fftcoeff*_scale), not at the
 * Evaluator level.
 *
 * Therefore, our CompatEvaluator must let standard SEAL update scale fields
 * normally. The ONLY adaptation needed is temporarily lowering ct.scale
 * before certain operations so that SEAL's internal "is_scale_within_bounds"
 * check passes. After the operation, we compute the correct scale value
 * (what SEAL would have computed with the original ct.scale) and set it.
 *
 * This file must be included AFTER "seal/seal.h".
 */

#include "seal/seal.h"
#include <NTL/RR.h>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

using namespace seal;

// ─── iter() macro ───
// In seal-modified, iter() returns an iterator adapter.
#define iter(x) (x)

/**
 * CompatEvaluator
 *
 * Wraps seal::Evaluator and provides the missing "reduced error" methods.
 * Strategy: temporarily lower ct.scale to pass SEAL's bounds checks,
 * then compute the correct scale after the operation.
 */
class CompatEvaluator {
public:
    Evaluator &evaluator;
    CKKSEncoder &encoder;
    const SEALContext &context;
    double default_scale;

    // Max exponent safe for pow(2.0, n) without overflow to infinity.
    // IEEE 754 double: max_exponent=1023, pow(2,1024)=inf.
    static constexpr int MAX_SAFE_POW2_EXP = 1023;

    CompatEvaluator(Evaluator &ev, CKKSEncoder &enc, const SEALContext &ctx, double scale)
        : evaluator(ev), encoder(enc), context(ctx), default_scale(scale) {}

    // ─── Safe scale helpers ───

    // Safe encode scale: ensures log2(enc_scale) < total_bits so CKKSEncoder::encode
    // doesn't throw "scale out of bounds". Returns the largest power-of-2 scale
    // that passes the encode bounds check. Falls back to default_scale if needed.
    double safe_encode_scale(const Ciphertext &ct) const {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return default_scale;
        int total_bits = ctx_data->total_coeff_modulus_bit_count();
        // CKKSEncoder::encode checks: log2(scale) < total_bits
        int safe_bits = std::min(total_bits - 2, MAX_SAFE_POW2_EXP);
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    double safe_encode_scale_for_parms(parms_id_type pid) const {
        auto ctx_data = context.get_context_data(pid);
        if (!ctx_data) return default_scale;
        int total_bits = ctx_data->total_coeff_modulus_bit_count();
        int safe_bits = std::min(total_bits - 2, MAX_SAFE_POW2_EXP);
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    // For multiply_plain: after lowering ct_scale and plain_scale, the result
    // ct.scale * plain.scale must also fit within total_bits.
    double safe_scale_for_multiply_plain(const Ciphertext &ct, double plain_scale) const {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return 1.0;
        int total_bits = ctx_data->total_coeff_modulus_bit_count();
        int plain_bits = (plain_scale > 1) ? (int)ceil(log2(plain_scale)) : 1;
        int max_result_bits = std::min(total_bits - 1, MAX_SAFE_POW2_EXP);
        int safe_bits = max_result_bits - plain_bits;
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    double safe_scale_for_square(const Ciphertext &ct) const {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return 1.0;
        int total_bits = ctx_data->total_coeff_modulus_bit_count();
        int max_result_bits = std::min(total_bits - 1, MAX_SAFE_POW2_EXP);
        int safe_bits = max_result_bits / 2;
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    double safe_scale_for_rescale(const Ciphertext &ct) const {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return 1.0;
        auto next_data = ctx_data->next_context_data();
        if (!next_data) return ct.scale();
        int next_bits = next_data->total_coeff_modulus_bit_count();
        int safe_bits = std::min(next_bits - 2, MAX_SAFE_POW2_EXP);
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    // ─── encode helper: ensures encode scale passes CKKSEncoder bounds ───
    // CKKSEncoder::encode checks: log2(scale) < total_coeff_modulus_bit_count()
    // If ct_scale is too large for the current level, we lower it.
    // Returns (enc_scale, needs_compensation) where:
    //   enc_scale = the scale we actually pass to encode
    //   needs_compensation = true if enc_scale < ct_scale (scale factor absorbed into plaintext value)
    double compute_safe_encode_scale(double desired_scale, parms_id_type pid) const {
        auto ctx_data = context.get_context_data(pid);
        if (!ctx_data) return desired_scale;
        int total_bits = ctx_data->total_coeff_modulus_bit_count();
        int desired_bits = (desired_scale > 1) ? (int)ceil(log2(desired_scale)) : 1;
        if (desired_bits < total_bits) return desired_scale;  // desired_scale fits
        // Too large: lower to safe value
        int safe_bits = std::min(total_bits - 2, MAX_SAFE_POW2_EXP);
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    // ─── multiply_vector_reduced_error ───
    // In seal-modified: encode at encrypted.scale(), multiply_plain.
    // Result scale = ct_scale * enc_scale.
    // If ct_scale is too large for encode, we lower the encode scale
    // and absorb the ratio into the plaintext values.

    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const std::vector<std::complex<double>> &coeff_vec,
        Ciphertext &destination)
    {
        if (coeff_vec.empty()) { destination = encrypted; return; }
        double ct_scale = encrypted.scale();
        double enc_scale = compute_safe_encode_scale(ct_scale, encrypted.parms_id());
        double scale_ratio = ct_scale / enc_scale;  // compensation factor
        std::vector<std::complex<double>> compensated_vec(coeff_vec.size());
        for (size_t i = 0; i < coeff_vec.size(); ++i) {
            compensated_vec[i] = coeff_vec[i] * scale_ratio;
        }
        Plaintext plain;
        encoder.encode(compensated_vec, encrypted.parms_id(), enc_scale, plain);
        Ciphertext ct_tmp = encrypted;
        ct_tmp.scale() = safe_scale_for_multiply_plain(encrypted, plain.scale());
        evaluator.multiply_plain(ct_tmp, plain, destination);
        // Result scale = ct_scale * enc_scale (since we encode at enc_scale and ct has ct_scale)
        // But the values are compensated by scale_ratio, so the effective scale is ct_scale * ct_scale
        destination.scale() = ct_scale * ct_scale;
    }

    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const Plaintext &plain,
        Ciphertext &destination)
    {
        double ct_scale = encrypted.scale();
        Ciphertext ct_tmp = encrypted;
        ct_tmp.scale() = safe_scale_for_multiply_plain(encrypted, plain.scale());
        evaluator.multiply_plain(ct_tmp, plain, destination);
        destination.scale() = ct_scale * plain.scale();
    }

    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const std::vector<Plaintext> &plain_vec,
        Ciphertext &destination)
    {
        if (plain_vec.empty()) return;
        double ct_scale = encrypted.scale();
        Ciphertext ct_tmp = encrypted;
        ct_tmp.scale() = safe_scale_for_multiply_plain(encrypted, plain_vec[0].scale());
        Ciphertext tmp;
        evaluator.multiply_plain(ct_tmp, plain_vec[0], destination);
        destination.scale() = ct_scale * plain_vec[0].scale();
        for (size_t i = 1; i < plain_vec.size(); ++i) {
            ct_tmp.scale() = safe_scale_for_multiply_plain(encrypted, plain_vec[i].scale());
            evaluator.multiply_plain(ct_tmp, plain_vec[i], tmp);
            tmp.scale() = ct_scale * plain_vec[i].scale();
            evaluator.add_inplace(destination, tmp);
        }
    }

    // ─── align_levels ───
    void align_levels(Ciphertext &ct1, Ciphertext &ct2) {
        if (ct1.parms_id() == ct2.parms_id()) return;
        auto lv1 = context.get_context_data(ct1.parms_id())->chain_index();
        auto lv2 = context.get_context_data(ct2.parms_id())->chain_index();
        if (lv1 > lv2) {
            evaluator.mod_switch_to_inplace(ct1, ct2.parms_id());
        } else {
            evaluator.mod_switch_to_inplace(ct2, ct1.parms_id());
        }
    }

    // ─── add_reduced_error variants ───
    void add_inplace_reduced_error(Ciphertext &ct1, const Ciphertext &ct2) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct2_copy = ct2;
            align_levels(ct1, ct2_copy);
        }
        evaluator.add_inplace(ct1, ct2);
    }
    void add_inplace_reduced_error(Ciphertext &ct1, Ciphertext &ct2) {
        align_levels(ct1, ct2);
        evaluator.add_inplace(ct1, ct2);
    }

    void add_reduced_error(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &result) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct1_copy = ct1, ct2_copy = ct2;
            align_levels(ct1_copy, ct2_copy);
            evaluator.add(ct1_copy, ct2_copy, result);
        } else {
            evaluator.add(ct1, ct2, result);
        }
    }
    void add_reduced_error(const Ciphertext &ct1, Ciphertext &ct2, Ciphertext &result) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct1_copy = ct1, ct2_copy = ct2;
            align_levels(ct1_copy, ct2_copy);
            evaluator.add(ct1_copy, ct2_copy, result);
        } else {
            evaluator.add(ct1, ct2, result);
        }
    }

    // ─── multiply_const ───
    void multiply_const(Ciphertext &ct, const NTL::RR &constant, Ciphertext &dest) {
        double ct_scale = ct.scale();
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), ct_scale, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain(ct_tmp, plain_const, dest);
        dest.scale() = ct_scale * ct_scale;
    }
    void multiply_const(Ciphertext &ct, double constant, Ciphertext &dest) {
        double ct_scale = ct.scale();
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), ct_scale, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain(ct_tmp, plain_const, dest);
        dest.scale() = ct_scale * ct_scale;
    }

    // ─── multiply_reduced_error: ct1 × ct2 with relinearization ───
    void multiply_reduced_error(
        const Ciphertext &ct1, const Ciphertext &ct2,
        const RelinKeys &rk, Ciphertext &dest)
    {
        double safe_scale = safe_scale_for_square(ct1);
        Ciphertext ct1_tmp = ct1, ct2_tmp = ct2;
        ct1_tmp.scale() = safe_scale;
        ct2_tmp.scale() = safe_scale;
        evaluator.multiply(ct1_tmp, ct2_tmp, dest);
        dest.scale() = ct1.scale() * ct2.scale();
        evaluator.relinearize_inplace(dest, rk);
    }
    void multiply_reduced_error(
        Ciphertext &ct1, Ciphertext &ct2,
        const RelinKeys &rk, Ciphertext &dest)
    {
        double safe_scale = safe_scale_for_square(ct1);
        Ciphertext ct1_tmp = ct1, ct2_tmp = ct2;
        ct1_tmp.scale() = safe_scale;
        ct2_tmp.scale() = safe_scale;
        evaluator.multiply(ct1_tmp, ct2_tmp, dest);
        dest.scale() = ct1.scale() * ct2.scale();
        evaluator.relinearize_inplace(dest, rk);
    }

    // ─── sub_inplace_reduced_error ───
    void sub_inplace_reduced_error(Ciphertext &ct1, const Ciphertext &ct2) {
        evaluator.sub_inplace(ct1, ct2);
    }

    // ─── Forward common Evaluator methods ───
    void rotate_vector(const Ciphertext &ct, int steps, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.rotate_vector(ct, steps, gk, dest);
    }
    void complex_conjugate(const Ciphertext &ct, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.complex_conjugate(ct, gk, dest);
    }

    void multiply_plain(const Ciphertext &ct, const Plaintext &pt, Ciphertext &dest) {
        double ct_scale = ct.scale();
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, pt.scale());
        evaluator.multiply_plain(ct_tmp, pt, dest);
        dest.scale() = ct_scale * pt.scale();
    }

    void add_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct2_copy = ct2;
            align_levels(ct1, ct2_copy);
        }
        evaluator.add_inplace(ct1, ct2);
    }
    void add(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &dest) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct1_copy = ct1, ct2_copy = ct2;
            align_levels(ct1_copy, ct2_copy);
            evaluator.add(ct1_copy, ct2_copy, dest);
        } else {
            evaluator.add(ct1, ct2, dest);
        }
    }
    void sub_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct2_copy = ct2;
            align_levels(ct1, ct2_copy);
        }
        evaluator.sub_inplace(ct1, ct2);
    }

    // ─── rescale_to_next_inplace ───
    // SEAL: divide_and_round_q_last + scale /= q_last.
    // We temporarily lower scale for bounds check, then compute correct scale.
    void rescale_to_next_inplace(Ciphertext &ct) {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return;
        auto next_data = ctx_data->next_context_data();
        if (!next_data) return;
        double orig_scale = ct.scale();
        ct.scale() = safe_scale_for_rescale(ct);
        evaluator.rescale_to_next_inplace(ct);
        // Compute q_last to get the correct result scale = orig_scale / q_last
        double q_last = static_cast<double>(ctx_data->parms().coeff_modulus().back().value());
        ct.scale() = orig_scale / q_last;
    }

    void relinearize_inplace(Ciphertext &ct, const RelinKeys &rk) {
        evaluator.relinearize_inplace(ct, rk);
    }

    // ─── mod_switch: preserves scale, only needs bounds bypass ───
    void mod_switch_to_next_inplace(Ciphertext &ct) {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data || !ctx_data->next_context_data()) return;
        double orig_scale = ct.scale();
        int next_bits = ctx_data->next_context_data()->total_coeff_modulus_bit_count();
        int safe_bits = std::min(next_bits - 2, MAX_SAFE_POW2_EXP);
        if (safe_bits < 10) safe_bits = 10;
        ct.scale() = pow(2.0, safe_bits);
        evaluator.mod_switch_to_next_inplace(ct);
        ct.scale() = orig_scale;
    }

    void mod_switch_to_inplace(Ciphertext &ct, parms_id_type pid) {
        double orig_scale = ct.scale();
        auto target_data = context.get_context_data(pid);
        if (!target_data) return;
        int target_bits = target_data->total_coeff_modulus_bit_count();
        int safe_bits = std::min(target_bits - 2, MAX_SAFE_POW2_EXP);
        if (safe_bits < 10) safe_bits = 10;
        ct.scale() = pow(2.0, safe_bits);
        evaluator.mod_switch_to_inplace(ct, pid);
        ct.scale() = orig_scale;
    }

    void mod_switch_to_inplace(Plaintext &pt, parms_id_type pid) {
        evaluator.mod_switch_to_inplace(pt, pid);
    }

    void rotate_rows(Ciphertext &ct, int steps, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.rotate_rows(ct, steps, gk, dest);
    }
    void rotate_columns(Ciphertext &ct, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.rotate_columns(ct, gk, dest);
    }

    // ─── square ───
    void square(const Ciphertext &ct, Ciphertext &dest) {
        double ct_scale = ct.scale();
        double safe_scale = safe_scale_for_square(ct);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale;
        evaluator.square(ct_tmp, dest);
        dest.scale() = ct_scale * ct_scale;
    }
    void square_inplace(Ciphertext &ct) {
        double ct_scale = ct.scale();
        double safe_scale = safe_scale_for_square(ct);
        ct.scale() = safe_scale;
        evaluator.square_inplace(ct);
        ct.scale() = ct_scale * ct_scale;
    }

    // ─── multiply_plain_inplace ───
    void multiply_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        double ct_scale = ct.scale();
        ct.scale() = safe_scale_for_multiply_plain(ct, pt.scale());
        evaluator.multiply_plain_inplace(ct, pt);
        ct.scale() = ct_scale * pt.scale();
    }

    void add_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        evaluator.add_plain_inplace(ct, pt);
    }
    void sub_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        evaluator.sub_plain_inplace(ct, pt);
    }

    // ─── multiply ct×ct ───
    void multiply(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &dest) {
        double safe_scale = safe_scale_for_square(ct1);
        Ciphertext ct1_tmp = ct1, ct2_tmp = ct2;
        ct1_tmp.scale() = safe_scale;
        ct2_tmp.scale() = safe_scale;
        evaluator.multiply(ct1_tmp, ct2_tmp, dest);
        dest.scale() = ct1.scale() * ct2.scale();
    }
    void multiply_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        double ct1_scale = ct1.scale();
        double safe_scale = safe_scale_for_square(ct1);
        ct1.scale() = safe_scale;
        Ciphertext ct2_tmp = ct2;
        ct2_tmp.scale() = safe_scale;
        evaluator.multiply_inplace(ct1, ct2_tmp);
        ct1.scale() = ct1_scale * ct2.scale();
    }

    void transform_to_ntt_inplace(Ciphertext &ct) {
        evaluator.transform_to_ntt_inplace(ct);
    }
    void transform_from_ntt_inplace(Ciphertext &ct) {
        evaluator.transform_from_ntt_inplace(ct);
    }
    void negate_inplace(Ciphertext &ct) {
        evaluator.negate_inplace(ct);
    }

    // ─── seal-modified specific methods ───

    // double_inplace: multiply by 2
    void double_inplace(Ciphertext &ct) {
        double ct_scale = ct.scale();
        Plaintext plain_two;
        encoder.encode(2.0, ct.parms_id(), ct_scale, plain_two);
        ct.scale() = safe_scale_for_multiply_plain(ct, plain_two.scale());
        evaluator.multiply_plain_inplace(ct, plain_two);
        ct.scale() = ct_scale * ct_scale;
    }

    // add_const
    void add_const(Ciphertext &ct, double constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain(ct, plain_const, dest);
    }
    void add_const(Ciphertext &ct, const NTL::RR &constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain(ct, plain_const, dest);
    }
    void add_const(Ciphertext &ct, double constant) {
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain_inplace(ct, plain_const);
    }

    void add_const_inplace(Ciphertext &ct, double constant) {
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain_inplace(ct, plain_const);
    }
    void add_const_inplace(Ciphertext &ct, const NTL::RR &constant) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain_inplace(ct, plain_const);
    }

    // ─── multiply_const_inplace variants ───
    void multiply_const_inplace(Ciphertext &ct, double constant) {
        double ct_scale = ct.scale();
        double enc_scale = ct_scale;
        auto ctx_data = context.get_context_data(ct.parms_id());
        int total_bits = ctx_data ? ctx_data->total_coeff_modulus_bit_count() : 0;
        int ct_bits = (ct_scale > 1) ? (int)ceil(log2(ct_scale)) : 1;
        if (2 * ct_bits >= total_bits) {
            enc_scale = safe_encode_scale(ct);
        }
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), enc_scale, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain_inplace(ct_tmp, plain_const);
        ct = ct_tmp;
        ct.scale() = ct_scale * enc_scale;
    }
    void multiply_const_inplace(Ciphertext &ct, const NTL::RR &constant) {
        double ct_scale = ct.scale();
        double enc_scale = ct_scale;
        auto ctx_data = context.get_context_data(ct.parms_id());
        int total_bits = ctx_data ? ctx_data->total_coeff_modulus_bit_count() : 0;
        int ct_bits = (ct_scale > 1) ? (int)ceil(log2(ct_scale)) : 1;
        if (2 * ct_bits >= total_bits) {
            enc_scale = safe_encode_scale(ct);
        }
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), enc_scale, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain_inplace(ct_tmp, plain_const);
        ct = ct_tmp;
        ct.scale() = ct_scale * enc_scale;
    }
    void multiply_const_inplace(Ciphertext &ct, double constant, Ciphertext &dest) {
        double ct_scale = ct.scale();
        double enc_scale = ct_scale;
        auto ctx_data = context.get_context_data(ct.parms_id());
        int total_bits = ctx_data ? ctx_data->total_coeff_modulus_bit_count() : 0;
        int ct_bits = (ct_scale > 1) ? (int)ceil(log2(ct_scale)) : 1;
        if (2 * ct_bits >= total_bits) {
            enc_scale = safe_encode_scale(ct);
        }
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), enc_scale, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain(ct_tmp, plain_const, dest);
        dest.scale() = ct_scale * enc_scale;
    }

    // multiply_const_inplace_fit_to_scale
    void multiply_const_inplace_fit_to_scale(Ciphertext &ct, const NTL::RR &constant, double scale_to_align) {
        double ct_scale = ct.scale();
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), scale_to_align, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain_inplace(ct_tmp, plain_const);
        ct = ct_tmp;
        ct.scale() = ct_scale * scale_to_align;
    }
    void multiply_const_inplace_fit_to_scale(Ciphertext &ct, const NTL::RR &constant, const NTL::RR &scale_to_align) {
        multiply_const_inplace_fit_to_scale(ct, constant, to_double(scale_to_align));
    }

    // ─── multiply_inplace_reduced_error ───
    void multiply_inplace_reduced_error(Ciphertext &ct1, Ciphertext &ct2, const RelinKeys &rk) {
        double ct1_scale = ct1.scale();
        double ct2_scale = ct2.scale();
        double safe_scale = safe_scale_for_square(ct1);
        ct1.scale() = safe_scale;
        ct2.scale() = safe_scale;
        evaluator.multiply_inplace(ct1, ct2);
        ct1.scale() = ct1_scale * ct2_scale;
        evaluator.relinearize_inplace(ct1, rk);
    }

    // ─── sub_reduced_error ───
    void sub_reduced_error(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &result) {
        evaluator.sub(ct1, ct2, result);
    }
    void sub_reduced_error(Ciphertext &ct1, Ciphertext &ct2, Ciphertext &result) {
        evaluator.sub(ct1, ct2, result);
    }
};
